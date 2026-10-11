"""Strict streaming FASTQ input, with gzip and paired-read validation."""
import gzip
from itertools import zip_longest
from pathlib import Path


def read_fastq(path):
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == '.gz' else open
    with opener(path, 'rt', encoding='ascii', newline=None) as handle:
        index = 0
        while True:
            header = handle.readline()
            if not header:
                return
            sequence, plus, quality = (handle.readline() for _ in range(3))
            if not sequence or not plus or not quality:
                raise ValueError(f'{path}: incomplete FASTQ record {index+1}')
            if not header.startswith('@') or not plus.startswith('+'):
                raise ValueError(f'{path}: invalid FASTQ header/separator at record {index+1}')
            if not header[1:].strip():
                raise ValueError(f'{path}: empty FASTQ read identifier at record {index+1}')
            sequence, quality = sequence.rstrip('\r\n').upper(), quality.rstrip('\r\n')
            if not sequence or len(sequence) != len(quality):
                raise ValueError(f'{path}: sequence/quality length mismatch at record {index+1}')
            if not set(sequence) <= set('ACGTRYKMSWBDHVN'):
                raise ValueError(f'{path}: invalid nucleotide at record {index+1}')
            if any(ord(q) < 33 or ord(q) > 126 for q in quality):
                raise ValueError(f'{path}: invalid Phred+33 quality at record {index+1}')
            yield {'read_id':index, 'name':header[1:].rstrip('\r\n'),
                   'sequence':sequence, 'quality':quality}
            index += 1


def paired_fastq(first, second, *, check_names=True):
    for a, b in zip_longest(read_fastq(first), read_fastq(second)):
        if a is None or b is None:
            raise ValueError('Paired FASTQ files contain different record counts')
        # CASAVA names end at whitespace; older headers use /1 and /2.
        def key(name):
            name = name.split()[0]
            return name[:-2] if name.endswith(('/1','/2')) else name
        if check_names and key(a['name']) != key(b['name']):
            raise ValueError(f'Paired FASTQ read names differ at record {a["read_id"]+1}')
        yield a, b
