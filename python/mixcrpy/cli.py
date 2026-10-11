"""Portable CLI; every runtime computation stays inside the Python package."""
from collections import Counter
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

from . import __version__
from .alignment import align_read
from .assembly import assemble, annotate_clonal_hits
from .fastq import read_fastq, paired_fastq
from .reference import load_reference
from .features import extract_features, boundary_map


def save_json(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')


def hit_text(hits):
    return ','.join(f'{h["gene"]}({h["score"]:g})' for h in hits)


def write_tsv(path, fields, records):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)


def export_alignments(alignments, path):
    fields = ['readId', 'allVHitsWithScore', 'allDHitsWithScore', 'allJHitsWithScore',
              'allCHitsWithScore', 'nSeqCDR3', 'aaSeqCDR3', 'qualCDR3', 'status']
    rows = []
    for a in alignments:
        rows.append(dict(readId=','.join(map(str,a.get('read_ids',[a['read_id']]))), nSeqCDR3=a['cdr3_nt'], aaSeqCDR3=a['cdr3_aa'],
                         qualCDR3=a['cdr3_quality'], status=a['status'],
                         **{f'all{t}HitsWithScore':hit_text(a[t.lower()+'_hits']) for t in 'VDJC'}))
    write_tsv(path, fields, rows)


def export_clones(result, prefix):
    fields = ['cloneId','cloneCount','cloneFraction','allVHitsWithScore','allDHitsWithScore',
              'allJHitsWithScore','allCHitsWithScore','nSeqCDR3','aaSeqCDR3','qualCDR3']
    rows = []
    feature_names = sorted({t['feature'] for c in result['clones'] for t in c.get('targets',[])})
    for feature in feature_names:
        if feature != 'CDR3':
            fields += ['nSeq'+feature, 'qual'+feature]
    for c in result['clones']:
        row = dict(cloneId=c['clone_id'], cloneCount=float(c['count']), cloneFraction=c['fraction'],
                         nSeqCDR3=c['cdr3_nt'], aaSeqCDR3=c['cdr3_aa'], qualCDR3=c['quality'],
                         **{f'all{t}HitsWithScore':hit_text(c[t.lower()+'_hits']) for t in 'VDJC'})
        for target in c.get('targets',[]):
            if target['feature'] != 'CDR3':
                row['nSeq'+target['feature']] = target['sequence']
                row['qual'+target['feature']] = target['quality']
        rows.append(row)
    write_tsv(str(prefix)+'.clones.tsv',fields,rows)
    assignments = []
    for a in result['read_assignments']:
        status = a['status']
        if status == 'core':
            if a.get('clustered'):
                status = 'clustered'
            elif a.get('preclustered'):
                status = 'preclustered'
        assignments.append(dict(readId=a['read_id'], cloneId='' if a['clone_id'] is None else a['clone_id'],
                                mappingType=status, reason=a.get('reason','')))
    write_tsv(str(prefix)+'.read-assignments.tsv',['readId','cloneId','mappingType','reason'],assignments)


def assemble_records(alignments, reference, chains, features=None, separate_by=()):
    options = {'clone_annotator':lambda c: annotate_clonal_hits(c,reference),
               'assemblingFeatures':features or ['CDR3']}
    options.update({'separateBy'+kind:True for kind in separate_by})
    return assemble(alignments, options)


def aligned_records(args, reference):
    options = {'allowPartialAlignments':args.allow_partial_alignments,
               'allowNoCDR3PartAlignments':args.allow_no_cdr3_part_alignments,
               'no_merge':args.no_merge}
    if args.read2:
        from .paired import align_pair
        for a,b in paired_fastq(args.fastq,args.read2,check_names=not args.pair_by_order):
            yield align_pair(a['read_id'],a['name'],a['sequence'],a['quality'],
                             b['name'],b['sequence'],b['quality'],reference,args.chains,options)
    else:
        for r in read_fastq(args.fastq):
            yield align_read(**r,reference=reference,chains=args.chains,
                             options={k:v for k,v in options.items() if k != 'no_merge'})


def read_alignment_document(path):
    data = json.loads(path.read_text(encoding='utf-8'))
    if data.get('format') != 'mixcrpy.alignments.v1':
        raise ValueError('Expected a mixcrpy.alignments.v1 document')
    for record in data['alignments']:
        for kind in 'vdjc':
            for hit in record.get(kind+'_hits',[]):
                for aligned in hit.get('alignments',[hit]):
                    if aligned is not None and 'operations' in aligned:
                        aligned['boundary_map'] = boundary_map(aligned)
                if 'operations' in hit:
                    hit['boundary_map'] = boundary_map(hit)
    return data


def assemble_contigs_document(data, clone_document, reference, prefix):
    from .contig import assemble_contig
    if clone_document.get('format') != 'mixcrpy.clones.v1':
        raise ValueError('Expected a mixcrpy.clones.v1 document')
    for key,default in (('chains','ALL'),('v_feature','VRegionWithP')):
        if clone_document.get(key,data.get(key,default)) != data.get(key,default):
            raise ValueError(f'Alignment and clone documents use different {key}')
    assignments = {a['read_id']:a['clone_id'] for a in clone_document['read_assignments']}
    evidence_ids = [i for r in data['alignments'] for i in r.get('read_ids',[r['read_id']])]
    if len(set(evidence_ids)) != len(evidence_ids) or set(evidence_ids) != set(assignments):
        raise ValueError('Alignment read identifiers must match clone assignments exactly')
    groups = {}
    for record in data['alignments']:
        ids = record.get('read_ids',[record['read_id']])
        clone_ids = {assignments.get(i) for i in ids}
        if len(clone_ids) > 1:
            raise ValueError('A merged alignment has inconsistent clone assignments')
        clone_id = next(iter(clone_ids))
        if clone_id is not None:
            groups.setdefault(clone_id,[]).append(record)
    rows, branches, fasta = [], [], []
    for clone in clone_document['clones']:
        if clone.get('assembling_features',['CDR3']) != ['CDR3']:
            raise NotImplementedError('Contig assembly currently requires the CDR3 clone feature')
        evidence = groups.get(clone['clone_id'],[])
        if not evidence:
            raise ValueError(f'No alignment evidence for clone {clone["clone_id"]}')
        for branch_id, branch in enumerate(assemble_contig(evidence,clone,reference)):
            branch.update(source_clone_id=clone['clone_id'],branch_id=branch_id,cdr3_nt=clone['cdr3_nt'])
            branches.append(branch)
            for target_id,target in enumerate(branch['targets']):
                name = f'clone{clone["clone_id"]}_branch{branch_id}_target{target_id}'
                fasta.append(f'>{name} count={branch["count"]:g}\n{target["sequence"]}\n')
                rows.append(dict(cloneId=clone['clone_id'],branchId=branch_id,targetId=target_id,
                                 count=branch['count'],sequence=target['sequence'],quality=target['quality'],
                                 containsCDR3=target_id==branch['cdr3_target'],
                                 **{f'all{kind}HitsWithScore':hit_text(branch[kind.lower()+'_hits']) for kind in 'VDJC'}))
    report = {'version':__version__,'input_clones':len(clone_document['clones']),
              'output_branches':len(branches),'output_targets':len(rows),
              'branch_count_sum':sum(b['count'] for b in branches)}
    save_json(str(prefix)+'.contigs.json',{'format':'mixcrpy.contigs.v1','branches':branches,'report':report})
    write_tsv(str(prefix)+'.contigs.tsv',['cloneId','branchId','targetId','count','sequence','quality','containsCDR3']+
              [f'all{kind}HitsWithScore' for kind in 'VDJC'],rows)
    Path(str(prefix)+'.contigs.fasta').write_bytes(''.join(fasta).encode('ascii'))
    save_json(str(prefix)+'.report.json',report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(prog='mixcrpy',description='V(D)J alignment, clone assembly and contig reconstruction')
    parser.add_argument('--version', action='version', version=__version__)
    sub = parser.add_subparsers(dest='command', required=True)
    for command in ('align','analyze'):
        p = sub.add_parser(command)
        p.add_argument('fastq', type=Path)
        p.add_argument('output_prefix',type=Path)
        p.add_argument('--read2',type=Path)
        p.add_argument('--allow-partial-alignments',action='store_true')
        p.add_argument('--allow-no-cdr3-part-alignments',action='store_true')
        p.add_argument('--no-merge',action='store_true')
        p.add_argument('--pair-by-order',action='store_true',help='Pair equal-length FASTQ files by record order, preserving both original headers')
        p.add_argument('--chains', default='ALL', choices=['ALL','TCR','BCR','TRA','TRB','TRG','TRD','IGH','IGK','IGL'])
        p.add_argument('--reference',type=Path)
        p.add_argument('--v-feature',choices=['VRegionWithP','VGeneWithP'],default='VRegionWithP')
        p.add_argument('--features',nargs='+',default=['CDR3'])
        p.add_argument('--separate-by',nargs='+',choices=['V','J','C'],default=[])
        p.add_argument('-f','--force',action='store_true')
    p = sub.add_parser('assemble')
    p.add_argument('alignments',type=Path)
    p.add_argument('output_prefix',type=Path)
    p.add_argument('--reference',type=Path)
    p.add_argument('--features',nargs='+',default=['CDR3'])
    p.add_argument('--separate-by',nargs='+',choices=['V','J','C'],default=[])
    p.add_argument('-f','--force',action='store_true')
    p = sub.add_parser('assemble-partial',help='Assemble overlapping partial alignments and ancillary targets')
    p.add_argument('alignments',type=Path)
    p.add_argument('output_prefix',type=Path)
    p.add_argument('--reference',type=Path)
    p.add_argument('--overlapped-only',action='store_true')
    p.add_argument('-f','--force',action='store_true')
    p = sub.add_parser('assemble-contigs',help='Reconstruct contig sequences from assigned CDR3 reads')
    p.add_argument('alignments',type=Path)
    p.add_argument('clones',type=Path)
    p.add_argument('output_prefix',type=Path)
    p.add_argument('--reference',type=Path)
    p.add_argument('-f','--force',action='store_true')
    args = parser.parse_args(argv)
    prefix = args.output_prefix
    suffixes = ['.clones.json','.clones.tsv','.read-assignments.tsv','.report.json']
    if args.command != 'assemble': suffixes += ['.alignments.json','.alignments.tsv']
    if args.command == 'assemble-contigs': suffixes = ['.contigs.json','.contigs.tsv','.contigs.fasta','.report.json']
    existing = [str(prefix)+s for s in suffixes if Path(str(prefix)+s).exists()]
    if existing and not args.force:
        parser.error(f'Output exists; use --force to replace: {existing[0]}')
    prefix.parent.mkdir(parents=True,exist_ok=True)
    try:
        data = read_alignment_document(args.alignments) if args.command in ('assemble','assemble-partial','assemble-contigs') else None
        v_feature = data.get('v_feature','VRegionWithP') if data else args.v_feature
        reference = load_reference(args.reference,v_feature=v_feature)
        if args.command == 'assemble-contigs':
            reference.select_chains(data.get('chains','ALL'))
            clones = json.loads(args.clones.read_text(encoding='utf-8'))
            print(json.dumps(assemble_contigs_document(data,clones,reference,prefix)))
            return
        if args.command == 'assemble-partial':
            from .partial import assemble_partial
            reference.select_chains(data.get('chains','ALL'))
            result = assemble_partial([a for a in data['alignments'] if a['status']=='aligned'],reference,
                                      {'chains':reference.selected_chains,'overlappedOnly':args.overlapped_only})
            save_json(str(prefix)+'.alignments.json',{'format':'mixcrpy.alignments.v1',
                      'chains':reference.selected_chains,'v_feature':v_feature,'alignments':result['alignments']})
            export_alignments(result['alignments'],str(prefix)+'.alignments.tsv')
            save_json(str(prefix)+'.report.json',result['report'])
            print(json.dumps(result['report']))
            return
        for name in args.features:
            if name not in reference.features:
                raise ValueError(f'Unsupported gene feature: {name}')
        if args.command == 'assemble':
            alignments = data['alignments']
            reference.select_chains(data.get('chains','ALL'))
        else:
            reference.select_chains(args.chains)
            alignments = list(aligned_records(args,reference))
            for a in alignments:
                extract_features(a,reference,args.features)
            save_json(str(prefix)+'.alignments.json',{'format':'mixcrpy.alignments.v1','chains':args.chains,
                      'v_feature':v_feature,'alignments':alignments})
            export_alignments(alignments,str(prefix)+'.alignments.tsv')
        report = {'version':__version__,'reference_version':'MiXCR 3.0.13 / repseqio.v1.6',
                  'chains':reference.selected_chains,
                  'v_feature':v_feature,'assembling_features':args.features,'separate_by':args.separate_by,
                  'alignment_status':dict(Counter(a['status'] for a in alignments)),
                  'total_reads':sum(a.get('number_of_reads',len(a.get('read_ids',[a['read_id']]))) for a in alignments),
                  'compatibility_notes':['Validated human alignment and assembly profiles, dataset coverage, and compatibility boundaries: https://github.com/thekingofall/Biosoft-Project-MiXCR-Python/blob/main/docs/validation.md']}
        if args.command != 'align':
            for a in alignments:
                extract_features(a,reference,args.features)
            result = assemble_records(alignments,reference,reference.selected_chains,args.features,args.separate_by)
            save_json(str(prefix)+'.clones.json',{'format':'mixcrpy.clones.v1',
                      'chains':reference.selected_chains,'v_feature':v_feature,**result})
            export_clones(result,prefix)
            report['assembly'] = result['report']
        save_json(str(prefix)+'.report.json',report)
        print(json.dumps(report,ensure_ascii=False))
    except (ValueError,NotImplementedError) as error:
        parser.exit(2,f'mixcrpy: {error}\n')


if __name__ == '__main__': main()
