"""Check the published FASTQ example through the installed command-line package."""
import json
from pathlib import Path
import sys


def main():
    folder = Path(sys.argv[1] if len(sys.argv) > 1 else 'results')
    for sample in ('single', 'paired'):
        document = json.loads((folder / f'{sample}.clones.json').read_text())
        counts = sorted(c['count'] for c in document['clones'])
        if counts != [8, 12]:
            raise SystemExit(f'{sample}: expected clone counts [8, 12], got {counts}')
        assigned = [r for r in document['read_assignments'] if r['clone_id'] is not None]
        if len(assigned) != 20:
            raise SystemExit(f'{sample}: expected 20 assigned reads, got {len(assigned)}')
    contigs = json.loads((folder / 'contigs.contigs.json').read_text())
    counts = sorted(b['count'] for b in contigs['branches'])
    if counts != [8, 12]:
        raise SystemExit(f'contigs: expected branch counts [8, 12], got {counts}')
    print('Example passed: single reads, paired reads, two clones, and two contig branches.')


if __name__ == '__main__':
    main()
