"""Sample size formulas for clinical research.

Every function returns a dict with the raw (unrounded) result and the
rounded-up integer sample sizes. Dropout adjustment is applied separately by
``adjust_dropout`` so that the unadjusted figure is always reported too.
"""
from __future__ import annotations

import math

from scipy.optimize import brentq
from scipy.stats import nct, norm
from scipy.stats import t as t_dist


def z(p: float) -> float:
    return float(norm.ppf(p))


def z_alpha(alpha: float, sides: int = 2) -> float:
    return z(1 - alpha / sides)


def _check_prob(name: str, v: float, open_interval: bool = True) -> None:
    if open_interval and not (0 < v < 1):
        raise ValueError(f"{name} must be between 0 and 1 (exclusive).")


def adjust_dropout(n: int, dropout: float) -> int:
    """N_adjusted = N / (1 - dropout), rounded up."""
    if not (0 <= dropout < 1):
        raise ValueError("Dropout must be in [0, 1).")
    return math.ceil(n / (1 - dropout) - 1e-9)


# 1. Two independent proportions --------------------------------------------
def two_proportions(p1, p2, alpha=0.05, power=0.80, ratio=1.0, sides=2,
                    continuity=False):
    """Normal approximation (pooled variance under H0), Fleiss et al. 2003.

    ratio = n2 / n1. Optional Fleiss continuity correction.
    """
    _check_prob("p1", p1); _check_prob("p2", p2); _check_prob("power", power)
    if p1 == p2:
        raise ValueError("p1 and p2 must differ.")
    k = ratio
    za, zb = z_alpha(alpha, sides), z(power)
    pbar = (p1 + k * p2) / (1 + k)
    d = abs(p1 - p2)
    n1 = (za * math.sqrt((1 + 1 / k) * pbar * (1 - pbar))
          + zb * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2) / k)) ** 2 / d ** 2
    if continuity:
        n1 = n1 / 4 * (1 + math.sqrt(1 + 2 * (k + 1) / (k * n1 * d))) ** 2
    n1c = math.ceil(n1 - 1e-9)
    n2c = math.ceil(n1c * k - 1e-9)
    return {"n1_raw": n1, "n1": n1c, "n2": n2c, "total": n1c + n2c}


# 2. Two independent means (t-test) ------------------------------------------
def t_power(d, n1, ratio=1.0, alpha=0.05, sides=2):
    """Power of the two-sample t-test (noncentral t), n1 may be fractional."""
    n2 = n1 * ratio
    df = n1 + n2 - 2
    ncp = d * math.sqrt(n1 * n2 / (n1 + n2))
    crit = t_dist.ppf(1 - alpha / sides, df)
    pw = nct.sf(crit, df, ncp)
    if sides == 2:
        pw += nct.cdf(-crit, df, ncp)
    return float(pw)


def two_means(diff=None, sd=None, d=None, alpha=0.05, power=0.80, ratio=1.0,
              sides=2):
    """Exact noncentral-t solution (same as G*Power / statsmodels)."""
    if d is None:
        if diff is None or sd is None or sd <= 0:
            raise ValueError("Give Cohen's d, or mean difference and SD.")
        d = abs(diff) / sd
    if d <= 0:
        raise ValueError("Effect size must be > 0.")
    lo = max(1.5, 1.0 + 1.0 / ratio)  # keeps df > 0
    hi = 4.0
    while t_power(d, hi, ratio, alpha, sides) < power:
        hi *= 2
        if hi > 1e8:
            raise ValueError("No solution below 100,000,000 per group.")
    n1 = brentq(lambda n: t_power(d, n, ratio, alpha, sides) - power, lo, hi,
                xtol=1e-10)
    n1c = math.ceil(n1 - 1e-9)
    n2c = math.ceil(n1c * ratio - 1e-9)
    return {"d": d, "n1_raw": float(n1), "n1": n1c, "n2": n2c,
            "total": n1c + n2c}


# 3. Diagnostic accuracy (Buderer 1996) --------------------------------------
def diagnostic_accuracy(sens, spec, prevalence, precision, conf=0.95):
    """Precision-based sample size for sensitivity and specificity."""
    for n_, v in (("sensitivity", sens), ("specificity", spec),
                  ("prevalence", prevalence), ("precision", precision)):
        _check_prob(n_, v)
    zz = z(1 - (1 - conf) / 2)
    n_dis = zz ** 2 * sens * (1 - sens) / precision ** 2
    n_nondis = zz ** 2 * spec * (1 - spec) / precision ** 2
    total_sens = n_dis / prevalence
    total_spec = n_nondis / (1 - prevalence)
    total = max(total_sens, total_spec)
    return {
        "diseased_needed": math.ceil(n_dis - 1e-9),
        "nondiseased_needed": math.ceil(n_nondis - 1e-9),
        "total_for_sens": math.ceil(total_sens - 1e-9),
        "total_for_spec": math.ceil(total_spec - 1e-9),
        "total": math.ceil(total - 1e-9),
        "driver": "sensitivity" if total_sens >= total_spec else "specificity",
    }


# 4. ROC AUC vs 0.5 (Hanley & McNeil 1982) ------------------------------------
def _auc_var(a: float, n_pos: float, n_neg: float) -> float:
    q1 = a / (2 - a)
    q2 = 2 * a ** 2 / (1 + a)
    return (a * (1 - a) + (n_pos - 1) * (q1 - a ** 2)
            + (n_neg - 1) * (q2 - a ** 2)) / (n_pos * n_neg)


def auc_power(auc, n_pos, n_neg, auc0=0.5, alpha=0.05, sides=2):
    se0 = math.sqrt(_auc_var(auc0, n_pos, n_neg))
    se1 = math.sqrt(_auc_var(auc, n_pos, n_neg))
    za = z_alpha(alpha, sides)
    return float(norm.cdf((abs(auc - auc0) - za * se0) / se1))


def roc_auc(auc, auc0=0.5, alpha=0.05, power=0.80, ratio=1.0, sides=2):
    """Smallest number of positives (cases) so that the test of
    H0: AUC = auc0 reaches the target power. ratio = negatives / positives."""
    if not (0.5 <= auc0 < auc < 1):
        raise ValueError("Need 0.5 <= AUC0 < AUC < 1.")
    n = 2
    while auc_power(auc, n, max(2, n * ratio), auc0, alpha, sides) < power:
        n += 1
        if n > 10 ** 6:
            raise ValueError("No solution below 1,000,000.")
    n_neg = math.ceil(n * ratio - 1e-9)
    return {"n_pos": n, "n_neg": n_neg, "total": n + n_neg,
            "achieved_power": auc_power(auc, n, n_neg, auc0, alpha, sides)}


# 5. Logistic regression / prediction models ----------------------------------
def logistic_epv(n_predictors, outcome_prop, epv=10):
    """Peduzzi 1996 rule: events = EPV x number of candidate parameters."""
    _check_prob("outcome proportion", outcome_prop)
    events = epv * n_predictors
    n = events / min(outcome_prop, 1 - outcome_prop)
    return {"events": events, "total": math.ceil(n - 1e-9)}


def riley_criteria(n_predictors, outcome_prop, r2_cs, shrinkage=0.9,
                   intercept_margin=0.05):
    """Riley et al. 2020 (BMJ) criteria (i) and (iii) for binary outcomes.

    (i)  n = p / ((S - 1) * ln(1 - R2_CS / S))
    (iii) n = (1.96 / margin)^2 * phi * (1 - phi)
    """
    _check_prob("outcome proportion", outcome_prop)
    if not (0 < r2_cs < shrinkage):
        raise ValueError("R2_CS must be > 0 and smaller than S.")
    n1 = n_predictors / ((shrinkage - 1) * math.log(1 - r2_cs / shrinkage))
    n3 = (1.959964 / intercept_margin) ** 2 * outcome_prop * (1 - outcome_prop)
    total = max(n1, n3)
    n_tot = math.ceil(total - 1e-9)
    return {"n_shrinkage": math.ceil(n1 - 1e-9),
            "n_intercept": math.ceil(n3 - 1e-9),
            "total": n_tot, "events": math.ceil(n_tot * outcome_prop - 1e-9),
            "epv": n_tot * outcome_prop / n_predictors}


# 6. Correlation (Fisher z) ----------------------------------------------------
def correlation(r, alpha=0.05, power=0.80, sides=2):
    if not (0 < abs(r) < 1):
        raise ValueError("r must be between -1 and 1, not 0.")
    c = 0.5 * math.log((1 + abs(r)) / (1 - abs(r)))
    n = ((z_alpha(alpha, sides) + z(power)) / c) ** 2 + 3
    return {"n_raw": n, "total": math.ceil(n - 1e-9)}


# 7. Survival: log-rank -------------------------------------------------------
def survival_logrank(hr, event_prob, alpha=0.05, power=0.80, alloc=0.5,
                     sides=2, method="schoenfeld"):
    """Required number of events (Schoenfeld 1983 or Freedman 1982) and total
    N = events / overall probability of observing an event."""
    if hr <= 0 or hr == 1:
        raise ValueError("HR must be > 0 and != 1.")
    _check_prob("event probability", event_prob); _check_prob("allocation", alloc)
    zz = (z_alpha(alpha, sides) + z(power)) ** 2
    if method == "schoenfeld":
        d = zz / (alloc * (1 - alloc) * math.log(hr) ** 2)
    elif method == "freedman":
        if alloc != 0.5:
            raise ValueError("Freedman formula here assumes 1:1 allocation.")
        d = zz * (1 + hr) ** 2 / (1 - hr) ** 2
    else:
        raise ValueError("method must be schoenfeld or freedman")
    events = math.ceil(d - 1e-9)
    total = math.ceil(events / event_prob - 1e-9)
    return {"events_raw": d, "events": events, "total": total}


# 8. Non-inferiority / equivalence -------------------------------------------
def ni_means(sd, margin, true_diff=0.0, alpha=0.025, power=0.80, ratio=1.0):
    """Non-inferiority (or superiority-by-margin), one-sided alpha.
    margin > 0 is the NI margin; true_diff = expected (test - control),
    with higher = better. n1 = (1 + 1/k) sd^2 (za + zb)^2 / (true_diff + margin)^2.
    """
    eff = true_diff + margin
    if sd <= 0 or eff <= 0:
        raise ValueError("Need SD > 0 and true difference + margin > 0.")
    n1 = (1 + 1 / ratio) * sd ** 2 * (z(1 - alpha) + z(power)) ** 2 / eff ** 2
    n1c = math.ceil(n1 - 1e-9); n2c = math.ceil(n1c * ratio - 1e-9)
    return {"n1_raw": n1, "n1": n1c, "n2": n2c, "total": n1c + n2c}


def equivalence_means(sd, margin, true_diff=0.0, alpha=0.05, power=0.80,
                      ratio=1.0):
    """TOST equivalence, normal approximation (Chow, Shao & Wang).
    With true_diff = 0: n1 = (1 + 1/k) sd^2 (z_{1-a} + z_{1-b/2})^2 / margin^2.
    """
    if sd <= 0 or margin <= abs(true_diff):
        raise ValueError("Need SD > 0 and margin > |true difference|.")
    zb = z(1 - (1 - power) / 2) if true_diff == 0 else z(power)
    n1 = ((1 + 1 / ratio) * sd ** 2 * (z(1 - alpha) + zb) ** 2
          / (margin - abs(true_diff)) ** 2)
    n1c = math.ceil(n1 - 1e-9); n2c = math.ceil(n1c * ratio - 1e-9)
    return {"n1_raw": n1, "n1": n1c, "n2": n2c, "total": n1c + n2c}


def ni_proportions(p_control, p_test, margin, alpha=0.025, power=0.80,
                   ratio=1.0):
    """Non-inferiority for two proportions (higher = better), one-sided alpha.
    n1 = (za + zb)^2 [p1(1-p1)/k + p2(1-p2)] / (p_test - p_control + margin)^2
    where n1 is the control group and k = n_test / n_control."""
    _check_prob("p_control", p_control); _check_prob("p_test", p_test)
    eff = p_test - p_control + margin
    if margin <= 0 or eff <= 0:
        raise ValueError("Need margin > 0 and (p_test - p_control + margin) > 0.")
    var = p_control * (1 - p_control) + p_test * (1 - p_test) / ratio
    n1 = (z(1 - alpha) + z(power)) ** 2 * var / eff ** 2
    n1c = math.ceil(n1 - 1e-9); n2c = math.ceil(n1c * ratio - 1e-9)
    return {"n1_raw": n1, "n1": n1c, "n2": n2c, "total": n1c + n2c}
