"""M5 birim testleri: istatistik motoru."""

import numpy as np
import pandas as pd

from microbiomeforge.stats import (
    alpha_diversity,
    benjamini_hochberg,
    bray_curtis_matrix,
    compare_alpha,
    differential_abundance,
    pcoa,
    permanova,
)


def _matrix():
    # 4 örnek, 2 grup (pre/post). pre: taxonA baskın; post: taxonB baskın.
    return pd.DataFrame({
        "taxon": ["A", "B", "C"],
        "pre1": [0.8, 0.1, 0.1],
        "pre2": [0.75, 0.15, 0.1],
        "post1": [0.1, 0.8, 0.1],
        "post2": [0.15, 0.75, 0.1],
    })


GROUPS = {"pre1": "pre", "pre2": "pre", "post1": "post", "post2": "post"}


def test_alpha_diversity_values():
    a = alpha_diversity(_matrix())
    assert set(a["sample"]) == {"pre1", "pre2", "post1", "post2"}
    assert (a["richness"] == 3).all()
    # eşit dağılım daha yüksek Shannon verir; buradaki örnekler tek-baskın → orta.
    assert (a["shannon"] > 0).all()


def test_shannon_uniform_is_max():
    uniform = pd.DataFrame({"taxon": ["A", "B"], "s": [0.5, 0.5]})
    a = alpha_diversity(uniform)
    assert abs(a.iloc[0]["shannon"] - np.log(2)) < 1e-3  # çıktı 4 haneye yuvarlı


def test_bray_curtis_symmetric_zero_diag():
    D = bray_curtis_matrix(_matrix())
    assert np.allclose(np.diag(D.to_numpy()), 0)
    assert np.allclose(D.to_numpy(), D.to_numpy().T)
    # pre1-pre2 mesafesi pre1-post1'den küçük olmalı
    assert D.loc["pre1", "pre2"] < D.loc["pre1", "post1"]


def test_pcoa_shape_and_explained():
    D = bray_curtis_matrix(_matrix())
    coords = pcoa(D, n_axes=2)
    assert list(coords.columns) == ["sample", "PCo1", "PCo2"]
    assert len(coords.attrs["explained"]) == 2


def test_permanova_separates_groups():
    D = bray_curtis_matrix(_matrix())
    res = permanova(D, GROUPS, permutations=199)
    assert res.n_groups == 2
    assert res.pseudo_f > 1  # gruplar ayrışıyor
    assert 0 < res.p_value <= 1


def test_compare_alpha_runs():
    a = alpha_diversity(_matrix())
    gt = compare_alpha(a, GROUPS, metric="shannon")
    assert gt.n_groups == 2
    assert gt.test == "Mann-Whitney U"


def test_benjamini_hochberg_monotone():
    p = np.array([0.001, 0.5, 0.02, 0.8])
    q = benjamini_hochberg(p)
    assert (q >= p).all() or np.allclose(q[0], p[0] * 4)
    assert (q <= 1).all()


def test_differential_abundance_detects_shift():
    df = differential_abundance(_matrix(), GROUPS)
    assert set(df.columns) >= {"taxon", "effect_clr", "p_value", "q_value", "significant"}
    # A ve B gruplar arası kayar; C sabit → C en yüksek q (en az anlamlı) civarı
    cval = df[df["taxon"] == "C"].iloc[0]
    assert cval["q_value"] >= df[df["taxon"] == "A"].iloc[0]["q_value"]
