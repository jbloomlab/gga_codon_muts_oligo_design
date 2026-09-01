"""Tests for gga_codon_muts_oligo_design."""

from pathlib import Path

import pandas as pd
import pytest

from gga_codon_muts_oligo_design import gga_codon_muts_oligo_design

CODON_FREQS_CSV = Path(__file__).parent / "human_codon_freq.csv"

# Minimal one-tile sequence encoding Met-Ala-Lys (no stop codons, no gaps).
# With initial_sequential_site=35: site 35=M, 36=A, 37=K.
# With initial_sequential_site=1:  site  1=M,  2=A,  3=K.
_NT_MAK = "ATGGCTAAA"


def _tiles_csv(tmp_path, nt_seq=_NT_MAK):
    p = tmp_path / "tiles.csv"
    pd.DataFrame(
        {
            "fragment": ["f1"],
            "fragment_sequence": [nt_seq],
            "inframe_mutated_region": [nt_seq],
        }
    ).to_csv(p, index=False)
    return p


def _mutations_csv(tmp_path, rows):
    p = tmp_path / "mutations.csv"
    pd.DataFrame(
        rows, columns=["sequential_site", "wildtype_aa", "mutant_aa", "representation"]
    ).to_csv(p, index=False)
    return p


def _run(tmp_path, mutations_rows, initial_sequential_site=1):
    out = tmp_path / "oligos.fa"
    gga_codon_muts_oligo_design(
        tiles_csv=_tiles_csv(tmp_path),
        mutations_to_make_csv=_mutations_csv(tmp_path, mutations_rows),
        output_oligos_fasta=out,
        max_representation=1,
        wildtype_frac=0,
        avoid_motifs=[],
        codon_freqs_csv=CODON_FREQS_CSV,
        initial_sequential_site=initial_sequential_site,
    )
    return out.read_text()


def test_default_site_numbering(tmp_path):
    """With default initial_sequential_site=1, site 2 (Ala) is named A2G."""
    content = _run(tmp_path, [{"sequential_site": 2, "wildtype_aa": "A", "mutant_aa": "G", "representation": 1}])
    assert "A2G" in content


def test_offset_site_numbering(tmp_path):
    """With initial_sequential_site=35, site 36 (Ala) is named A36G."""
    content = _run(
        tmp_path,
        [{"sequential_site": 36, "wildtype_aa": "A", "mutant_aa": "G", "representation": 1}],
        initial_sequential_site=35,
    )
    assert "A36G" in content


def test_site_below_range_raises(tmp_path):
    """A sequential_site below initial_sequential_site raises ValueError."""
    with pytest.raises(ValueError, match="outside range"):
        _run(
            tmp_path,
            [{"sequential_site": 34, "wildtype_aa": "M", "mutant_aa": "A", "representation": 1}],
            initial_sequential_site=35,
        )


def test_site_above_range_raises(tmp_path):
    """A sequential_site above the last tile site raises ValueError."""
    with pytest.raises(ValueError, match="outside range"):
        _run(
            tmp_path,
            [{"sequential_site": 38, "wildtype_aa": "K", "mutant_aa": "A", "representation": 1}],
            initial_sequential_site=35,
        )
