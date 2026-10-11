"""Read feature extraction using the RepSeq.IO partitioning coordinates."""
import re


def targets(record):
    return record.get('targets') or [{
        'sequence': record.get('oriented_sequence', record.get('sequence', '')),
        'quality': record.get('oriented_quality', record.get('quality', ''))}]


def target_alignment(hit, index):
    alignments = hit.get('alignments')
    if alignments is not None:
        return alignments[index] if 0 <= index < len(alignments) else None
    return hit if index == 0 else None


def boundary_map(hit):
    """Restore maps from serialized JSON, including earlier format v1 files."""
    if 'boundary_map' in hit:
        return {int(k): v for k, v in hit['boundary_map'].items()}
    i, j = hit['ref_start'], hit['query_start']
    result = {i: j}
    for op in hit['operations']:
        if op in 'MX':
            i += 1
            j += 1
        elif op == 'I':
            j += 1
        elif op == 'D':
            i += 1
        result[i] = j
    return result


def point_position(record, reference, point, target_index=0):
    match = re.fullmatch(r'(.+?)([+-]\d+)?', point)
    base, offset = match.group(1), int(match.group(2) or 0)
    descriptor = reference.point_descriptors.get(base)
    if descriptor is None:
        raise ValueError(f'Unknown reference point: {point}')
    hits = record.get(descriptor['type'].lower()+'_hits', [])
    if not hits:
        return None
    best = hits[0]
    hit = target_alignment(best, target_index)
    if hit is None:
        return None
    gene = reference.by_name[best['gene']]
    coordinates = gene.get('feature_points', gene['points'])
    if descriptor['attached']:
        left = descriptor['left']
        position = hit['ref_start' if left else 'ref_end']
        activation = descriptor.get('activation')
        if activation is not None:
            act_match = re.fullmatch(r'(.+?)([+-]\d+)?', activation)
            active = coordinates.get(act_match.group(1), -1)
            if active >= 0:
                active += int(act_match.group(2) or 0)
            if active < 0 or (position > active if left else position < active):
                return None
    else:
        position = coordinates.get(base, -1)
        if position < 0:
            return None
    mapped = boundary_map(hit).get(position+offset)
    if mapped is not None and mapped < -1:
        mapped = -2-mapped
    return None if mapped == -1 else mapped


def feature_range(record, reference, name, target_index=0):
    descriptor = reference.features.get(name)
    if descriptor is None:
        raise ValueError(f'Unsupported gene feature: {name}')
    start = point_position(record, reference, descriptor['begin'], target_index)
    end = point_position(record, reference, descriptor['end'], target_index)
    if start is None or end is None or end < start:
        return None
    sequence = targets(record)[target_index]['sequence']
    return (start, end) if 0 <= start <= end <= len(sequence) else None


def extract_feature(record, reference, name):
    descriptor = reference.features.get(name)
    if descriptor is None:
        raise ValueError(f'Unsupported gene feature: {name}')
    selected = None
    for index, target in enumerate(targets(record)):
        region = feature_range(record, reference, name, index)
        if region is None:
            continue
        start, end = region
        quality = target['quality'][start:end]
        score = min(map(ord, quality), default=127)
        if selected is not None and selected[0] >= score:
            continue
        value = {'sequence': target['sequence'][start:end], 'quality': quality,
                 'target_index': index, 'target_range': [start,end],
                 'contains_junction': descriptor['contains_junction'],
                 'contains_d': descriptor['contains_d'],
                 'reference_points': [descriptor['begin'], descriptor['end']]}
        if descriptor['contains_cdr3']:
            cdr3 = feature_range(record, reference, 'CDR3', index)
            if cdr3 is not None and start <= cdr3[0] <= cdr3[1] <= end:
                value['cdr3_range'] = [cdr3[0]-start,cdr3[1]-start]
        selected = score, value
    if selected is None:
        return None
    value = selected[1]
    # MiXCR 3.0.13 getRelativeRange uses getTargetContainingFeature, whose
    # quality accumulator stays -1: it selects the last available target.
    # getFeature independently chooses the first target with maximum quality.
    n_regions = [region for junction in ('VDJunction', 'DJJunction')
                 if (region := relative_range(record, reference, name, junction)) is not None]
    if not n_regions:
        region = relative_range(record, reference, name, 'VJJunction')
        if region is not None:
            n_regions.append(region)
    value['n_regions'] = n_regions
    return value


def relative_range(record, reference, feature, subfeature):
    """Port VDJCObject.getRelativeRange, including its target-selection behavior."""
    selected = None
    for index, target in enumerate(targets(record)):
        region = feature_range(record, reference, feature, index)
        if region is not None:
            quality = target['quality'][region[0]:region[1]]
            if min((ord(q)-33 for q in quality), default=127) > -1:
                selected = index, region
    if selected is None:
        return None
    index, (start, end) = selected
    subregion = feature_range(record, reference, subfeature, index)
    if subregion is None or not start <= subregion[0] <= subregion[1] <= end:
        return None
    return [subregion[0]-start, subregion[1]-start]


def extract_features(record, reference, names):
    """Attach available requested features; absent features remain absent."""
    record['features'] = {name: value for name in names
                          if (value := extract_feature(record, reference, name)) is not None}
    return record['features']
