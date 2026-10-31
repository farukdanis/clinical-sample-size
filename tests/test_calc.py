import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import calc  # noqa: E402


# Two proportions: textbook example p1=0.40, p2=0.60, alpha=0.05, power=0.80
# -> 97 per group without, 107 with Fleiss continuity correction.
def test_two_proportions_textbook():
    assert calc.two_proportions(0.4, 0.6)["n1"] == 97
    assert calc.two_proportions(0.4, 0.6, continuity=True)["n1"] == 107


def test_two_proportions_unequal_ratio_power_check():
    # Simulated power of a pooled z-test at the computed n should be ~0.80+
    rng = np.random.default_rng(1)
    r = calc.two_proportions(0.3, 0.5, ratio=2)
    n1, n2 = r["n1"], r["n2"]
    x1 = rng.binomial(n1, 0.3, 20000); x2 = rng.binomial(n2, 0.5, 20000)
    p = (x1 + x2) / (n1 + n2)
    zstat = (x2 / n2 - x1 / n1) / np.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    assert (np.abs(zstat) > 1.959964).mean() > 0.78


# Two means: G*Power, d = 0.5 -> 64 per group; d = 0.8 -> 26 per group
def test_two_means_gpower():
    assert calc.two_means(d=0.5)["n1"] == 64
    assert calc.two_means(d=0.8)["n1"] == 26
    assert calc.two_means(diff=5, sd=10)["n1"] == 64


# Buderer: Se 0.90, precision 0.05, prevalence 0.20 -> 139 diseased, 692 total
def test_buderer():
    r = calc.diagnostic_accuracy(0.90, 0.85, 0.20, 0.05)
    assert r["diseased_needed"] == 139
    assert r["total_for_sens"] == 692
    assert r["total_for_spec"] == 245
    assert r["total"] == 692 and r["driver"] == "sensitivity"


# Correlation: Hulley et al. table, r = 0.30 -> 85
def test_correlation_hulley():
    assert calc.correlation(0.30)["total"] == 85


# Schoenfeld: HR 0.70, 1:1 -> 247 events; Freedman -> 253
def test_survival():
    assert calc.survival_logrank(0.7, 0.5)["events"] == 247
    assert calc.survival_logrank(0.7, 0.5)["total"] == 494
    assert calc.survival_logrank(0.7, 0.5, method="freedman")["events"] == 253


def test_roc_minimal_n():
    r = calc.roc_auc(0.70)
    assert r["achieved_power"] >= 0.80
    n = r["n_pos"]
    assert calc.auc_power(0.70, n - 1, n - 1) < 0.80
    # Larger AUC needs fewer cases
    assert calc.roc_auc(0.80)["n_pos"] < n


def test_epv_and_riley():
    assert calc.logistic_epv(10, 0.20)["total"] == 500
    r = calc.riley_criteria(10, 0.2, 0.1)
    expected = 10 / ((0.9 - 1) * math.log(1 - 0.1 / 0.9))
    assert r["n_shrinkage"] == math.ceil(expected)
    assert r["n_intercept"] == math.ceil((1.959964 / 0.05) ** 2 * 0.16)


def test_noninferiority_and_equivalence():
    assert calc.ni_means(1.0, 0.5)["n1"] == 63
    assert calc.equivalence_means(1.0, 0.5)["n1"] == 69
    r = calc.ni_proportions(0.8, 0.8, 0.1)
    assert r["n1"] == math.ceil((1.959964 + 0.841621) ** 2 * 0.32 / 0.01)


def test_dropout():
    assert calc.adjust_dropout(100, 0.10) == 112
    assert calc.adjust_dropout(90, 0.10) == 100
    with pytest.raises(ValueError):
        calc.adjust_dropout(10, 1.0)


def test_input_validation():
    with pytest.raises(ValueError):
        calc.two_proportions(0.5, 0.5)
    with pytest.raises(ValueError):
        calc.diagnostic_accuracy(1.2, 0.9, 0.1, 0.05)
