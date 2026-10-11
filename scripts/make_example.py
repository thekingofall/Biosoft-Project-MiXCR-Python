"""Build a small deterministic human TRB example from the bundled reference."""
import json
from pathlib import Path
from importlib.resources import files


def main():
    genes = {g['name']: g for g in json.loads(
        files('mixcrpy').joinpath('data/human_reference.json').read_text())['genes']}
    output = Path(__file__).resolve().parents[1] / 'examples'
    output.mkdir(exist_ok=True)
    records = []
    for index, (v_name, j_name, count, junction) in enumerate([
        ('TRBV12-3*00', 'TRBJ2-7*00', 12, 'GACAGGGG'),
        ('TRBV7-9*00', 'TRBJ1-2*00', 8, 'GCGTACCTAGG'),
    ], 1):
        v, j = genes[v_name], genes[j_name]
        start, end = v['points']['CDR3Begin'] - 87, v['points']['VEnd'] - 2
        j_start = j['points']['JBegin'] + 2
        length = end - v['points']['CDR3Begin'] + len(junction) + j['points']['CDR3End'] - j_start
        junction += 'A' * (-length % 3)
        sequence = v['sequence'][start:end] + junction + j['sequence'][j_start:]
        records.extend((f'template{index}_{k:02d}', sequence) for k in range(count))
    complement = str.maketrans('ACGTN', 'TGCAN')
    for name, mate in [('trb.fastq', 0), ('trb_R1.fastq', 1), ('trb_R2.fastq', 2)]:
        contents = []
        for identifier, sequence in records:
            if mate == 1:
                sequence = sequence[:110]
            elif mate == 2:
                sequence = sequence[-110:].translate(complement)[::-1]
            suffix = f'/{mate}' if mate else ''
            contents.append(f'@{identifier}{suffix}\n{sequence}\n+\n{"I" * len(sequence)}\n')
        (output / name).write_bytes(''.join(contents).encode('ascii'))
    print('Created 20 synthetic TRB reads and 20 matching read pairs.')


if __name__ == '__main__':
    main()
