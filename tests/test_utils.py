import os
import tempfile
from pathlib import Path
import pytest
import hashlib
import stat

from viralrecall import utils


class FakeGenome:
    def __init__(self, seqs: dict):
        self._seqs = seqs
    def keys(self):
        return list(self._seqs.keys())
    def __getitem__(self, key):
        return self._seqs[key]


def test_check_directory_permissions_true(tmp_path):
    assert utils.check_directory_permissions(tmp_path) is True


def test_check_directory_permissions_false(tmp_path):
    d = tmp_path / 'noaccess'
    d.mkdir()
    # remove read/write for owner
    d.chmod(0)
    try:
        assert utils.check_directory_permissions(d) is False
    finally:
        # restore permissions so pytest can clean up
        d.chmod(stat.S_IRWXU)


def test_mp_cpu_none_and_bounds(monkeypatch):
    # patch cpu_count to deterministic value
    monkeypatch.setattr(utils, 'cpu_count', lambda: 4)
    assert utils.mp_cpu(None) == 4
    assert utils.mp_cpu(2) == 2
    assert utils.mp_cpu(8) == 4


def test_find_db_not_dir(tmp_path):
    with pytest.raises(FileNotFoundError):
        utils.find_db(tmp_path / 'nope')


def test_find_db_missing_files(tmp_path):
    d = tmp_path / 'hmm'
    d.mkdir()
    # missing files
    with pytest.raises(FileNotFoundError):
        utils.find_db(d)


def test_find_db_success(tmp_path):
    d = tmp_path / 'hmm'
    d.mkdir()
    (d / 'gvog_mirus_cat.hmm').write_text('a')
    (d / 'NCLDV_markers.hmm').write_text('b')
    (d / 'gvog_annotation.tsv').write_text('c')
    assert utils.find_db(d) == d

