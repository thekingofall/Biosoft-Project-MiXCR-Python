import gzip
import pytest
from mixcrpy.fastq import read_fastq, paired_fastq


def test_gzip_and_crlf(tmp_path):
    p = tmp_path/'r.fastq.gz'
    with gzip.open(p,'wb') as f: f.write(b'@one\r\nAcGtN\r\n+\r\nIIIII\r\n')
    assert list(read_fastq(p))[0]['sequence'] == 'ACGTN'


@pytest.mark.parametrize('data', [b'@x\nAC\n+\n', b'x\nAC\n+\nII\n',b'@x\nAC\n+\nI\n',b'@x\nAZ\n+\nII\n',b'@\nAC\n+\nII\n'])
def test_invalid_fastq(tmp_path,data):
    p = tmp_path/'r.fastq'; p.write_bytes(data)
    with pytest.raises(ValueError): list(read_fastq(p))


def test_pair_validation(tmp_path):
    a,b = tmp_path/'a',tmp_path/'b'
    a.write_bytes(b'@read/1\nAC\n+\nII\n'); b.write_bytes(b'@read/2\nGT\n+\nII\n')
    assert len(list(paired_fastq(a,b))) == 1
    b.write_bytes(b'@other/2\nGT\n+\nII\n')
    with pytest.raises(ValueError): list(paired_fastq(a,b))
    b.write_bytes(b'')
    with pytest.raises(ValueError): list(paired_fastq(a,b))
