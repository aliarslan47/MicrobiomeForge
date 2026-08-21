"""M0 birim testleri: referans kayıt defteri."""

import re

import pytest

from microbiomeforge.references import CATALOG, ToolRegistry

# Basit ama gerçek DOI biçim doğrulaması (10.xxxx/...).
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")


def test_all_catalog_dois_wellformed():
    for key, ref in CATALOG.items():
        assert DOI_RE.match(ref.doi), f"{key}: geçersiz DOI biçimi {ref.doi!r}"
        assert ref.year >= 2015
        assert ref.journal and ref.name and ref.role


def test_bibliography_only_contains_used_tools():
    reg = ToolRegistry()
    reg.mark_used("kraken2", version="2.1.3")
    reg.mark_used("sylph", version="0.6.1")
    biblio = reg.bibliography()
    keys = {r.key for r in biblio}
    assert keys == {"kraken2", "sylph"}
    # Kullanılmayan araç kaynakçaya girmez.
    assert "metaphlan4" not in keys


def test_bibliography_sorted_and_numbered():
    reg = ToolRegistry()
    reg.mark_used("sylph", "0.6.1")       # 2024
    reg.mark_used("bracken", "2.9")       # 2017
    reg.mark_used("kraken2", "2.1.3")     # 2019
    lines = reg.formatted_bibliography()
    assert lines[0].startswith("[1] ") and "Bracken" in lines[0]  # en eski önce
    assert lines[-1].startswith("[3] ") and "sylph" in lines[-1]
    assert "v0.6.1" in lines[-1]


def test_unknown_tool_rejected():
    reg = ToolRegistry()
    with pytest.raises(KeyError):
        reg.mark_used("madeup_tool_9000")
