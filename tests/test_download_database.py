from pathlib import Path
import hashlib
import tempfile
from viralrecall import download_database


def test_check_file_hash(tmp_path):
    f = tmp_path / 'tmp.txt'
    data = b'TestData123'
    f.write_bytes(data)
    m = hashlib.md5(data).hexdigest()
    assert download_database.check_file_hash(f, m) is True
    assert download_database.check_file_hash(f, 'deadbeef') is False
