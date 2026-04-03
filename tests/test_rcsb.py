#!/usr/bin/env python3
# coding: utf-8

import os

from pdb_cpp import Coor, rcsb
import pdb_cpp._pyprops as _pyprops

from .datafiles import MMCIF_1Y0M


class _DummyResponse:
    def __init__(self, payload):
        self.payload = payload

    def read(self):
        return self.payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


def test_build_download_url_asymmetric_unit():
    assert (
        rcsb.build_download_url("4hhb")
        == "https://files.rcsb.org/download/4hhb.cif"
    )
    assert (
        rcsb.build_download_url("4hhb", file_format="pdb")
        == "https://files.rcsb.org/download/4hhb.pdb"
    )


def test_build_download_url_biological_assembly():
    assert (
        rcsb.build_download_url(
            "5a9z",
            structure="biological_assembly",
            assembly_id=1,
        )
        == "https://files.rcsb.org/download/5a9z-assembly1.cif"
    )
    assert (
        rcsb.build_download_url(
            "1hh3",
            structure="assembly",
            file_format="pdb",
            assembly_id=2,
        )
        == "https://files.rcsb.org/download/1hh3.pdb2"
    )


def test_download_structure_uses_cache(monkeypatch, tmp_path):
    calls = []
    payload = b"data_test\n#\n"

    def fake_urlopen(url):
        calls.append(url)
        return _DummyResponse(payload)

    monkeypatch.setattr(rcsb.urllib.request, "urlopen", fake_urlopen)

    first_path = rcsb.download_structure("1abc", cache_dir=tmp_path)
    second_path = rcsb.download_structure("1abc", cache_dir=tmp_path)

    assert first_path == second_path
    assert os.path.exists(first_path)
    assert calls == ["https://files.rcsb.org/download/1abc.cif"]


def test_load_structure_returns_coor(monkeypatch, tmp_path):
    with open(MMCIF_1Y0M, "rb") as handle:
        payload = handle.read()

    def fake_urlopen(url):
        return _DummyResponse(payload)

    monkeypatch.setattr(rcsb.urllib.request, "urlopen", fake_urlopen)

    coor = rcsb.load("1y0m", cache_dir=tmp_path)

    assert isinstance(coor, Coor)
    assert coor.len == 648
    assert coor.model_num == 1


def test_coor_init_accepts_rcsb_options(monkeypatch, tmp_path):
    captured = {}

    def fake_download_structure(
        pdb_id,
        structure="asymmetric_unit",
        file_format="cif",
        assembly_id=1,
        cache_dir=None,
        force_download=False,
    ):
        captured["pdb_id"] = pdb_id
        captured["structure"] = structure
        captured["file_format"] = file_format
        captured["assembly_id"] = assembly_id
        captured["cache_dir"] = cache_dir
        captured["force_download"] = force_download
        return MMCIF_1Y0M

    monkeypatch.setattr(_pyprops.rcsb, "download_structure", fake_download_structure)

    coor = Coor(
        pdb_id="1y0m",
        rcsb_structure="biological_assembly",
        assembly_id=2,
        cache_dir=str(tmp_path),
        force_download=True,
    )

    assert coor.len == 648
    assert captured == {
        "pdb_id": "1y0m",
        "structure": "biological_assembly",
        "file_format": "cif",
        "assembly_id": 2,
        "cache_dir": str(tmp_path),
        "force_download": True,
    }