"""Tests for the aslscan fixture registry, geometry helpers and phantoms.

None of these need a generated fixture or the aslscan binary unless marked otherwise.
"""

import json

import pytest

from aslprep.tests import aslscan_fixtures as af


def test_spec_file_matches():
    """The committed spec file (the CI cache key) must match the registry."""
    committed = af.SPEC_FILE.read_text()
    assert committed == af.spec_text(), (
        f'{af.SPEC_FILE.name} is out of date with the fixture registry. Regenerate it with:\n'
        f'    {af.SPEC_COMMAND}'
    )


def test_hashed_modules_exist():
    for mod in af.HASHED_MODULES:
        assert (af.TESTS_DIR / mod).is_file(), mod


def test_digest_dir_line_endings(tmp_path):
    """Line endings do not change a digest; names and lengths do."""
    lf, crlf = tmp_path / 'lf', tmp_path / 'crlf'
    for d, sep in ((lf, b'\n'), (crlf, b'\r\n')):
        (d / 'sub').mkdir(parents=True)
        (d / 'a.json').write_bytes(b'{' + sep + b'"x": 1' + sep + b'}' + sep)
        (d / 'sub' / 'b.tsv').write_bytes(b'volume_type' + sep + b'control' + sep)
    assert af._digest_dir(lf) == af._digest_dir(crlf)

    before = af._digest_dir(lf)
    (lf / 'sub' / 'b.tsv').rename(lf / 'sub' / 'c.tsv')
    assert af._digest_dir(lf) != before

    before = af._digest_dir(lf)
    (lf / 'a.json').write_bytes(b'{\n"x": 10\n}\n')
    assert af._digest_dir(lf) != before


def test_check_aslscan_requires_matching_stamp(tmp_path):
    binary = tmp_path / 'aslscan'
    binary.write_bytes(b'not really a binary')
    with pytest.raises(af.AslscanUnavailable, match=r'no aslscan\.build\.json'):
        af.check_aslscan(binary)

    stamp = {**af._expected_stamp(), 'sha256': af._sha256_file(binary)}
    (tmp_path / af.STAMP_NAME).write_text(json.dumps(stamp))
    assert af.check_aslscan(binary) == binary

    binary.write_bytes(b'rebuilt from something else')
    with pytest.raises(af.AslscanUnavailable, match='does not match the hash'):
        af.check_aslscan(binary)

    stamp = {**af._expected_stamp(), 'aslscan': '0' * 40, 'sha256': af._sha256_file(binary)}
    (tmp_path / af.STAMP_NAME).write_text(json.dumps(stamp))
    with pytest.raises(af.AslscanUnavailable, match='was built from'):
        af.check_aslscan(binary)


def test_find_aslscan_reports_every_candidate(tmp_path, monkeypatch):
    monkeypatch.setenv('ASLSCAN', str(tmp_path / 'missing'))
    monkeypatch.setenv('PATH', str(tmp_path))
    with pytest.raises(af.AslscanUnavailable, match='missing does not exist'):
        af.find_aslscan()


def test_spec_digest_tracks_inputs(monkeypatch):
    before = af.spec_digest()
    monkeypatch.setattr(af, 'CACHE_EPOCH', af.CACHE_EPOCH + 1)
    assert af.spec_digest() != before
