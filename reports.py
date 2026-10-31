"""Report builders shared by the Gradio app and the static (Pyodide) app.

Author: Faruk Danış, MD (ORCID 0000-0001-8722-5402)
"""
import calc

REFS = {
    "fleiss": "Fleiss JL, Levin B, Paik MC. *Statistical Methods for Rates and Proportions*. 3rd ed. Hoboken (NJ): Wiley; 2003.",
    "cohen": "Cohen J. *Statistical Power Analysis for the Behavioral Sciences*. 2nd ed. Hillsdale (NJ): Lawrence Erlbaum; 1988.",
    "buderer": "Buderer NM. Statistical methodology: I. Incorporating the prevalence of disease into the sample size calculation for sensitivity and specificity. *Acad Emerg Med*. 1996;3(9):895-900.",
    "hanley": "Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) curve. *Radiology*. 1982;143(1):29-36.",
    "peduzzi": "Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR. A simulation study of the number of events per variable in logistic regression analysis. *J Clin Epidemiol*. 1996;49(12):1373-9.",
    "riley": "Riley RD, Ensor J, Snell KIE, et al. Calculating the sample size required for developing a clinical prediction model. *BMJ*. 2020;368:m441.",
    "hulley": "Hulley SB, Cummings SR, Browner WS, Grady DG, Newman TB. *Designing Clinical Research*. 4th ed. Philadelphia: Lippincott Williams & Wilkins; 2013.",
    "schoenfeld": "Schoenfeld DA. Sample-size formula for the proportional-hazards regression model. *Biometrics*. 1983;39(2):499-503.",
    "freedman": "Freedman LS. Tables of the number of patients required in clinical trials using the logrank test. *Stat Med*. 1982;1(2):121-9.",
    "chow": "Chow SC, Shao J, Wang H. *Sample Size Calculations in Clinical Research*. 2nd ed. Boca Raton (FL): Chapman & Hall/CRC; 2008.",
}

LATEX = [{"left": "$$", "right": "$$", "display": True},
         {"left": "$", "right": "$", "display": False}]


def trn(x, fmt="g"):
    """Turkish number formatting (decimal comma)."""
    return format(x, fmt).replace(".", ",")


def pct(x):
    return f"{x * 100:g}%"


def report(result_rows, formula, assumptions, method_en, method_tr, refs):
    rows = "\n".join(f"| {a} | **{b}** |" for a, b in result_rows)
    refs_md = "\n".join(f"{i}. {REFS[r]}" for i, r in enumerate(refs, 1))
    return f"""### Sonuç / Result
| | |
|---|---|
{rows}

### Formül / Formula
{formula}

### Varsayımlar / Assumptions
{assumptions}

### Yöntem cümlesi / Methods sentence
**EN:** {method_en}

**TR:** {method_tr}

### Kaynak / Reference
{refs_md}
"""


def safe(fn):
    def wrapper(*args):
        try:
            return fn(*args)
        except ValueError as e:
            return f"⚠️ **Girdi hatası / Input error:** {e}"
    return wrapper


def drop_rows(n_total, dropout):
    adj = calc.adjust_dropout(n_total, dropout)
    return adj, [(f"Kayıp düzeltmeli toplam / Total after {pct(dropout)} dropout", adj)]


# ---------------------------------------------------------------- handlers
@safe
def h_props(p1, p2, alpha, power, ratio, sides, cc, dropout):
    r = calc.two_proportions(p1, p2, alpha, power, ratio, int(sides), cc)
    adj, extra = drop_rows(r["total"], dropout)
    cc_en = " with Fleiss continuity correction" if cc else ""
    cc_tr = " (Fleiss süreklilik düzeltmesiyle)" if cc else ""
    return report(
        [("Grup 1 / Group 1", r["n1"]), ("Grup 2 / Group 2", r["n2"]),
         ("Toplam / Total", r["total"])] + extra,
        r"$$n_1=\frac{\left[z_{1-\alpha/2}\sqrt{(1+\tfrac1k)\bar p\bar q}+z_{1-\beta}\sqrt{p_1q_1+\tfrac{p_2q_2}{k}}\right]^2}{(p_1-p_2)^2},\quad k=\frac{n_2}{n_1}$$"
        + ("\n\nFleiss CC: $n'=\\frac{n}{4}\\left(1+\\sqrt{1+\\frac{2(k+1)}{k\\,n\\,|p_1-p_2|}}\\right)^2$" if cc else ""),
        "- Bağımsız iki grup, ikili sonuç / independent groups, binary outcome\n- Normal yaklaşım; H0 altında havuzlanmış varyans / normal approximation, pooled variance under H0",
        f"Assuming an outcome proportion of {pct(p1)} in group 1 and {pct(p2)} in group 2, a {int(sides)}-sided α of {alpha}, {pct(power)} power and an allocation ratio of 1:{ratio:g}, {r['n1']} and {r['n2']} participants are required{cc_en} ({r['total']} in total); allowing for {pct(dropout)} attrition, we will enrol {adj} participants.",
        f"Grup 1'de %{trn(p1*100)}, grup 2'de %{trn(p2*100)} sonuç oranı varsayılarak, {int(sides)} yönlü α={trn(alpha)}, %{trn(power*100)} güç ve 1:{trn(ratio)} dağılım oranıyla{cc_tr} gruplarda sırasıyla {r['n1']} ve {r['n2']} (toplam {r['total']}) katılımcı gerekmektedir; %{trn(dropout*100)} kayıp öngörülerek {adj} katılımcı alınacaktır.",
        ["fleiss"])


@safe
def h_means(mode, diff, sd, d, alpha, power, ratio, sides, dropout):
    if mode.startswith("Cohen"):
        r = calc.two_means(d=d, alpha=alpha, power=power, ratio=ratio, sides=int(sides))
        eff_en = f"a standardized effect size (Cohen's d) of {d:g}"
        eff_tr = f"standart etki büyüklüğü (Cohen d) {trn(d)}"
    else:
        r = calc.two_means(diff=diff, sd=sd, alpha=alpha, power=power, ratio=ratio, sides=int(sides))
        eff_en = f"a mean difference of {diff:g} with a common SD of {sd:g} (d = {r['d']:.2f})"
        eff_tr = f"{trn(diff)} birimlik ortalama farkı ve {trn(sd)} ortak standart sapma (d = {trn(r['d'], '.2f')})"
    adj, extra = drop_rows(r["total"], dropout)
    return report(
        [("Cohen d", f"{r['d']:.3f}"), ("Grup 1 / Group 1", r["n1"]),
         ("Grup 2 / Group 2", r["n2"]), ("Toplam / Total", r["total"])] + extra,
        r"$$n_1\approx\left(1+\tfrac1k\right)\frac{(z_{1-\alpha/2}+z_{1-\beta})^2}{d^2},\quad d=\frac{|\mu_1-\mu_2|}{\sigma}$$" "\n\nKesin çözüm: merkezî olmayan t dağılımı. / Exact solution: noncentral t distribution (as in G*Power).",
        "- Bağımsız t-testi, normal dağılım, eşit varyans / independent t-test, normality, equal variances",
        f"To detect {eff_en} with a {int(sides)}-sided α of {alpha} and {pct(power)} power (allocation 1:{ratio:g}), {r['n1']} and {r['n2']} participants are required ({r['total']} in total; independent-samples t-test). Allowing for {pct(dropout)} attrition, {adj} participants will be enrolled.",
        f"{eff_tr} saptamak için {int(sides)} yönlü α={trn(alpha)}, %{trn(power*100)} güç ve 1:{trn(ratio)} dağılımla gruplarda {r['n1']} ve {r['n2']} (toplam {r['total']}) katılımcı gerekmektedir (bağımsız örneklem t-testi). %{trn(dropout*100)} kayıp öngörülerek {adj} katılımcı alınacaktır.",
        ["cohen"])


@safe
def h_diag(sens, spec, prev, prec, conf, dropout):
    r = calc.diagnostic_accuracy(sens, spec, prev, prec, conf)
    adj, extra = drop_rows(r["total"], dropout)
    return report(
        [("Hasta sayısı (duyarlılık) / Diseased needed", r["diseased_needed"]),
         ("Hasta olmayan (özgüllük) / Non-diseased needed", r["nondiseased_needed"]),
         ("Toplam (duyarlılık için) / Total for sensitivity", r["total_for_sens"]),
         ("Toplam (özgüllük için) / Total for specificity", r["total_for_spec"]),
         ("Gereken toplam / Required total", r["total"])] + extra,
        r"$$n_{Se}=\frac{z^2_{1-\alpha/2}\,Se(1-Se)}{W^2\cdot P},\qquad n_{Sp}=\frac{z^2_{1-\alpha/2}\,Sp(1-Sp)}{W^2\cdot(1-P)}$$" "\n\n$W$: güven aralığı yarı genişliği / CI half-width; $P$: prevalans / prevalence.",
        "- Ardışık (kesitsel) hasta alımı; prevalans örneklemde de geçerli / consecutive enrolment, sample prevalence equals P\n- Wald güven aralığı / Wald confidence interval",
        f"Assuming an expected sensitivity of {pct(sens)}, specificity of {pct(spec)} and a disease prevalence of {pct(prev)}, {r['total']} consecutive patients are required to estimate {r['driver']} with a {pct(conf)} confidence interval half-width of {prec:g} (Buderer method). Allowing for {pct(dropout)} exclusions, {adj} patients will be enrolled.",
        f"Beklenen duyarlılık %{trn(sens*100)}, özgüllük %{trn(spec*100)} ve hastalık prevalansı %{trn(prev*100)} varsayılarak, {'duyarlılığın' if r['driver']=='sensitivity' else 'özgüllüğün'} ±{trn(prec)} yarı genişlikli %{trn(conf*100)} güven aralığıyla kestirilmesi için {r['total']} ardışık hasta gerekmektedir (Buderer yöntemi). %{trn(dropout*100)} dışlanma öngörülerek {adj} hasta alınacaktır.",
        ["buderer"])


@safe
def h_auc(auc, auc0, alpha, power, ratio, sides, dropout):
    r = calc.roc_auc(auc, auc0, alpha, power, ratio, int(sides))
    adj, extra = drop_rows(r["total"], dropout)
    return report(
        [("Pozitif (olgu) / Positives (cases)", r["n_pos"]),
         ("Negatif (kontrol) / Negatives", r["n_neg"]),
         ("Toplam / Total", r["total"]),
         ("Ulaşılan güç / Achieved power", f"{r['achieved_power']:.3f}")] + extra,
        r"$$SE^2(A)=\frac{A(1-A)+(n_A-1)(Q_1-A^2)+(n_N-1)(Q_2-A^2)}{n_A n_N},\; Q_1=\frac{A}{2-A},\; Q_2=\frac{2A^2}{1+A}$$"
        "\n\n$$|A-A_0|\\ge z_{1-\\alpha/2}\\,SE(A_0)+z_{1-\\beta}\\,SE(A)$$\n\nEn küçük $n_A$ sayısal olarak bulunur. / The smallest $n_A$ is found numerically.",
        "- Tek test, H0: AUC = AUC0 / single test against a fixed AUC0\n- Hanley–McNeil üstel varyans yaklaşımı / Hanley–McNeil exponential variance approximation",
        f"To show that an AUC of {auc:g} differs from {auc0:g} with a {int(sides)}-sided α of {alpha} and {pct(power)} power, {r['n_pos']} patients with and {r['n_neg']} patients without the outcome are required ({r['total']} in total; Hanley–McNeil variance). Allowing for {pct(dropout)} attrition, {adj} patients will be enrolled.",
        f"{trn(auc)} düzeyindeki bir AUC'nin {trn(auc0)} değerinden farkını {int(sides)} yönlü α={trn(alpha)} ve %{trn(power*100)} güçle göstermek için sonucu olan {r['n_pos']} ve olmayan {r['n_neg']} (toplam {r['total']}) hasta gerekmektedir (Hanley–McNeil varyansı). %{trn(dropout*100)} kayıp öngörülerek {adj} hasta alınacaktır.",
        ["hanley"])


@safe
def h_logit(k, prop, epv, r2, dropout):
    e = calc.logistic_epv(int(k), prop, epv)
    rows = [(f"EPV={epv:g} → olay / events", e["events"]),
            (f"EPV={epv:g} → toplam / total", e["total"])]
    refs = ["peduzzi"]
    riley_en = riley_tr = ""
    total = e["total"]
    if r2 and r2 > 0:
        rl = calc.riley_criteria(int(k), prop, r2)
        rows += [("Riley (i) büzülme ≤ %10 / shrinkage", rl["n_shrinkage"]),
                 ("Riley (iii) kesişim ±0.05 / intercept", rl["n_intercept"]),
                 ("Riley toplam / total", rl["total"]),
                 ("Riley olay / events (EPP)", f"{rl['events']} ({rl['epv']:.1f})")]
        refs.append("riley")
        total = max(total, rl["total"])
        riley_en = f" Following Riley et al., with an anticipated Cox–Snell R² of {r2:g}, at least {rl['total']} participants ({rl['events']} events) are needed to limit overfitting (expected shrinkage ≥ 0.90) and to estimate the overall risk precisely."
        riley_tr = f" Riley ve ark. yöntemine göre, öngörülen Cox–Snell R² {trn(r2)} iken aşırı uyumu sınırlamak (beklenen büzülme ≥ 0,90) ve genel riski kesin kestirmek için en az {rl['total']} katılımcı ({rl['events']} olay) gerekir."
    adj, extra = drop_rows(total, dropout)
    rows += [("Gereken toplam / Required total", total)] + extra
    return report(
        rows,
        r"$$N_{EPV}=\frac{EPV\times k}{\min(\phi,1-\phi)}$$" "\n\n"
        r"Riley (i): $$N=\frac{k}{(S-1)\ln\left(1-\frac{R^2_{CS}}{S}\right)},\ S=0.9\qquad$$ Riley (iii): $$N=\left(\frac{1.96}{0.05}\right)^2\phi(1-\phi)$$",
        "- $k$: aday parametre sayısı (kategorik değişkenlerde gösterge sayısı) / candidate parameters (dummy variables count separately)\n- EPV ≥ 10 kaba bir kuraldır; tahmin modelleri için Riley ölçütleri önerilir / EPV ≥ 10 is a rule of thumb; Riley criteria are preferred for prediction models",
        f"With {int(k)} candidate predictor parameters and an expected outcome proportion of {pct(prop)}, at least {e['events']} events ({e['total']} participants) are required to keep {epv:g} events per variable.{riley_en} Allowing for {pct(dropout)} missing data, {adj} participants will be included.",
        f"{int(k)} aday parametre ve %{trn(prop*100)} beklenen sonuç oranıyla değişken başına {trn(epv)} olay sağlamak için en az {e['events']} olay ({e['total']} katılımcı) gerekmektedir.{riley_tr} %{trn(dropout*100)} eksik veri öngörülerek {adj} katılımcı alınacaktır.",
        refs)


@safe
def h_corr(r, alpha, power, sides, dropout):
    res = calc.correlation(r, alpha, power, int(sides))
    adj, extra = drop_rows(res["total"], dropout)
    return report(
        [("Toplam / Total", res["total"])] + extra,
        r"$$n=\left(\frac{z_{1-\alpha/2}+z_{1-\beta}}{C}\right)^2+3,\qquad C=\tfrac12\ln\frac{1+r}{1-r}$$",
        "- İki değişkenli normal dağılım; H0: ρ = 0 / bivariate normality; H0: ρ = 0\n- Fisher z dönüşümü / Fisher z transformation",
        f"To detect a correlation coefficient of {r:g} with a {int(sides)}-sided α of {alpha} and {pct(power)} power, {res['total']} participants are required (Fisher z approximation); allowing for {pct(dropout)} attrition, {adj} participants will be enrolled.",
        f"{trn(r)} düzeyinde bir korelasyon katsayısını {int(sides)} yönlü α={trn(alpha)} ve %{trn(power*100)} güçle saptamak için {res['total']} katılımcı gerekmektedir (Fisher z yaklaşımı); %{trn(dropout*100)} kayıp öngörülerek {adj} katılımcı alınacaktır.",
        ["hulley"])


@safe
def h_surv(hr, pev, alpha, power, alloc, sides, method, dropout):
    m = "schoenfeld" if method.startswith("Schoenfeld") else "freedman"
    r = calc.survival_logrank(hr, pev, alpha, power, alloc, int(sides), m)
    adj, extra = drop_rows(r["total"], dropout)
    f = (r"$$D=\frac{(z_{1-\alpha/2}+z_{1-\beta})^2}{\pi(1-\pi)(\ln HR)^2}$$" if m == "schoenfeld"
         else r"$$D=(z_{1-\alpha/2}+z_{1-\beta})^2\frac{(1+HR)^2}{(1-HR)^2}$$")
    f += "\n\n$$N=D/P(\\text{olay/event})$$"
    return report(
        [("Gereken olay / Required events", r["events"]),
         ("Toplam hasta / Total patients", r["total"])] + extra,
        f,
        "- Orantılı hazardlar / proportional hazards\n- $\\pi$: grup 1'e ayrılma oranı / allocation proportion\n- Olay olasılığı izlem süresi ve kayıp yapısına bağlıdır / event probability depends on follow-up",
        f"To detect a hazard ratio of {hr:g} with a {int(sides)}-sided α of {alpha} and {pct(power)} power, {r['events']} events are required ({'Schoenfeld' if m=='schoenfeld' else 'Freedman'} method). With an anticipated overall event probability of {pct(pev)}, {r['total']} patients are needed; allowing for {pct(dropout)} loss to follow-up, {adj} patients will be enrolled.",
        f"{trn(hr)} düzeyinde bir hazard oranını {int(sides)} yönlü α={trn(alpha)} ve %{trn(power*100)} güçle saptamak için {r['events']} olay gerekmektedir ({'Schoenfeld' if m=='schoenfeld' else 'Freedman'} yöntemi). Öngörülen genel olay olasılığı %{trn(pev*100)} olduğundan {r['total']} hasta gerekir; %{trn(dropout*100)} izlem kaybı öngörülerek {adj} hasta alınacaktır.",
        [m])


@safe
def h_ni(kind, sd, margin, tdiff, pc, pt, alpha, power, ratio, dropout):
    if kind.startswith("NI — ortalama"):
        r = calc.ni_means(sd, margin, tdiff, alpha, power, ratio)
        f = r"$$n_1=\left(1+\tfrac1k\right)\frac{\sigma^2(z_{1-\alpha}+z_{1-\beta})^2}{(\Delta+\delta)^2}$$"
        en = f"Using a non-inferiority margin of {margin:g} (SD {sd:g}, expected difference {tdiff:g}), a one-sided α of {alpha} and {pct(power)} power"
        tr = f"{trn(margin)} daha kötü olmama (non-inferiority) sınırı (SS {trn(sd)}, beklenen fark {trn(tdiff)}), tek yönlü α={trn(alpha)} ve %{trn(power*100)} güçle"
    elif kind.startswith("Eşdeğerlik"):
        r = calc.equivalence_means(sd, margin, tdiff, alpha, power, ratio)
        f = r"$$n_1=\left(1+\tfrac1k\right)\frac{\sigma^2(z_{1-\alpha}+z_{1-\beta/2})^2}{\delta^2}\quad(\Delta=0)$$" "\n\n$\\Delta\\neq0$ ise $z_{1-\\beta}$ ve $(\\delta-|\\Delta|)^2$ kullanılır (yaklaşık) / if $\\Delta\\neq0$, $z_{1-\\beta}$ and $(\\delta-|\\Delta|)^2$ are used (approximate)."
        en = f"Using equivalence margins of ±{margin:g} (SD {sd:g}, expected difference {tdiff:g}), two one-sided tests at α = {alpha} and {pct(power)} power"
        tr = f"±{trn(margin)} eşdeğerlik sınırları (SS {trn(sd)}, beklenen fark {trn(tdiff)}), α={trn(alpha)} düzeyinde iki tek yönlü test (TOST) ve %{trn(power*100)} güçle"
    else:
        r = calc.ni_proportions(pc, pt, margin, alpha, power, ratio)
        f = r"$$n_1=\frac{(z_{1-\alpha}+z_{1-\beta})^2\left[p_c(1-p_c)+\tfrac{p_t(1-p_t)}{k}\right]}{(p_t-p_c+\delta)^2}$$"
        en = f"Assuming success rates of {pct(pc)} (control) and {pct(pt)} (test), a non-inferiority margin of {margin*100:g} percentage points, a one-sided α of {alpha} and {pct(power)} power"
        tr = f"Kontrolde %{trn(pc*100)}, test grubunda %{trn(pt*100)} başarı oranı, {trn(margin*100)} yüzde puanlık daha kötü olmama (non-inferiority) sınırı, tek yönlü α={trn(alpha)} ve %{trn(power*100)} güçle"
    adj, extra = drop_rows(r["total"], dropout)
    return report(
        [("Grup 1 (kontrol) / Group 1 (control)", r["n1"]),
         ("Grup 2 (test) / Group 2 (test)", r["n2"]), ("Toplam / Total", r["total"])] + extra,
        f,
        "- Yüksek değer = daha iyi / higher = better\n- Normal yaklaşım / normal approximation\n- Sınır (δ) klinik gerekçeyle önceden belirlenmeli / the margin must be pre-specified on clinical grounds",
        f"{en}, {r['n1']} and {r['n2']} participants are required ({r['total']} in total); allowing for {pct(dropout)} attrition, {adj} participants will be enrolled.",
        f"{tr} gruplarda {r['n1']} ve {r['n2']} (toplam {r['total']}) katılımcı gerekmektedir; %{trn(dropout*100)} kayıp öngörülerek {adj} katılımcı alınacaktır.",
        ["chow"])
