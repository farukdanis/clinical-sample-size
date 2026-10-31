"""Clinical Sample Size Calculator — Gradio app (TR/EN).

Author: Faruk Danış, MD (ORCID 0000-0001-8722-5402)
"""
import gradio as gr

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
        f"Grup 1'de %{p1*100:g}, grup 2'de %{p2*100:g} sonuç oranı varsayılarak, {int(sides)} yönlü α={alpha}, %{power*100:g} güç ve 1:{ratio:g} dağılım oranıyla{cc_tr} gruplarda sırasıyla {r['n1']} ve {r['n2']} (toplam {r['total']}) katılımcı gerekmektedir; %{dropout*100:g} kayıp öngörülerek {adj} katılımcı alınacaktır.",
        ["fleiss"])


@safe
def h_means(mode, diff, sd, d, alpha, power, ratio, sides, dropout):
    if mode.startswith("Cohen"):
        r = calc.two_means(d=d, alpha=alpha, power=power, ratio=ratio, sides=int(sides))
        eff_en = f"a standardized effect size (Cohen's d) of {d:g}"
        eff_tr = f"standart etki büyüklüğü (Cohen d) {d:g}"
    else:
        r = calc.two_means(diff=diff, sd=sd, alpha=alpha, power=power, ratio=ratio, sides=int(sides))
        eff_en = f"a mean difference of {diff:g} with a common SD of {sd:g} (d = {r['d']:.2f})"
        eff_tr = f"{diff:g} birimlik ortalama farkı ve {sd:g} ortak standart sapma (d = {r['d']:.2f})"
    adj, extra = drop_rows(r["total"], dropout)
    return report(
        [("Cohen d", f"{r['d']:.3f}"), ("Grup 1 / Group 1", r["n1"]),
         ("Grup 2 / Group 2", r["n2"]), ("Toplam / Total", r["total"])] + extra,
        r"$$n_1\approx\left(1+\tfrac1k\right)\frac{(z_{1-\alpha/2}+z_{1-\beta})^2}{d^2},\quad d=\frac{|\mu_1-\mu_2|}{\sigma}$$" "\n\nKesin çözüm: merkezî olmayan t dağılımı (statsmodels). / Exact solution: noncentral t distribution (statsmodels).",
        "- Bağımsız t-testi, normal dağılım, eşit varyans / independent t-test, normality, equal variances",
        f"To detect {eff_en} with a {int(sides)}-sided α of {alpha} and {pct(power)} power (allocation 1:{ratio:g}), {r['n1']} and {r['n2']} participants are required ({r['total']} in total; independent-samples t-test). Allowing for {pct(dropout)} attrition, {adj} participants will be enrolled.",
        f"{eff_tr} saptamak için {int(sides)} yönlü α={alpha}, %{power*100:g} güç ve 1:{ratio:g} dağılımla gruplarda {r['n1']} ve {r['n2']} (toplam {r['total']}) katılımcı gerekmektedir (bağımsız örneklem t-testi). %{dropout*100:g} kayıp öngörülerek {adj} katılımcı alınacaktır.",
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
        f"Beklenen duyarlılık %{sens*100:g}, özgüllük %{spec*100:g} ve hastalık prevalansı %{prev*100:g} varsayılarak, {'duyarlılığın' if r['driver']=='sensitivity' else 'özgüllüğün'} ±{prec:g} yarı genişlikli %{conf*100:g} güven aralığıyla kestirilmesi için {r['total']} ardışık hasta gerekmektedir (Buderer yöntemi). %{dropout*100:g} dışlanma öngörülerek {adj} hasta alınacaktır.",
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
        f"{auc:g} düzeyindeki bir AUC'nin {auc0:g} değerinden farkını {int(sides)} yönlü α={alpha} ve %{power*100:g} güçle göstermek için sonucu olan {r['n_pos']} ve olmayan {r['n_neg']} (toplam {r['total']}) hasta gerekmektedir (Hanley–McNeil varyansı). %{dropout*100:g} kayıp öngörülerek {adj} hasta alınacaktır.",
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
        riley_tr = f" Riley ve ark. yöntemine göre, öngörülen Cox–Snell R² {r2:g} iken aşırı uyumu sınırlamak (beklenen büzülme ≥ 0,90) ve genel riski kesin kestirmek için en az {rl['total']} katılımcı ({rl['events']} olay) gerekir."
    adj, extra = drop_rows(total, dropout)
    rows += [("Gereken toplam / Required total", total)] + extra
    return report(
        rows,
        r"$$N_{EPV}=\frac{EPV\times k}{\min(\phi,1-\phi)}$$" "\n\n"
        r"Riley (i): $$N=\frac{k}{(S-1)\ln\left(1-\frac{R^2_{CS}}{S}\right)},\ S=0.9\qquad$$ Riley (iii): $$N=\left(\frac{1.96}{0.05}\right)^2\phi(1-\phi)$$",
        "- $k$: aday parametre sayısı (kategorik değişkenlerde gösterge sayısı) / candidate parameters (dummy variables count separately)\n- EPV ≥ 10 kaba bir kuraldır; tahmin modelleri için Riley ölçütleri önerilir / EPV ≥ 10 is a rule of thumb; Riley criteria are preferred for prediction models",
        f"With {int(k)} candidate predictor parameters and an expected outcome proportion of {pct(prop)}, at least {e['events']} events ({e['total']} participants) are required to keep {epv:g} events per variable.{riley_en} Allowing for {pct(dropout)} missing data, {adj} participants will be included.",
        f"{int(k)} aday parametre ve %{prop*100:g} beklenen sonuç oranıyla değişken başına {epv:g} olay sağlamak için en az {e['events']} olay ({e['total']} katılımcı) gerekmektedir.{riley_tr} %{dropout*100:g} eksik veri öngörülerek {adj} katılımcı alınacaktır.",
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
        f"{r:g} düzeyinde bir korelasyon katsayısını {int(sides)} yönlü α={alpha} ve %{power*100:g} güçle saptamak için {res['total']} katılımcı gerekmektedir (Fisher z yaklaşımı); %{dropout*100:g} kayıp öngörülerek {adj} katılımcı alınacaktır.",
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
        f"{hr:g} düzeyinde bir hazard oranını {int(sides)} yönlü α={alpha} ve %{power*100:g} güçle saptamak için {r['events']} olay gerekmektedir ({'Schoenfeld' if m=='schoenfeld' else 'Freedman'} yöntemi). Öngörülen genel olay olasılığı %{pev*100:g} olduğundan {r['total']} hasta gerekir; %{dropout*100:g} izlem kaybı öngörülerek {adj} hasta alınacaktır.",
        [m])


@safe
def h_ni(kind, sd, margin, tdiff, pc, pt, alpha, power, ratio, dropout):
    if kind.startswith("NI — ortalama"):
        r = calc.ni_means(sd, margin, tdiff, alpha, power, ratio)
        f = r"$$n_1=\left(1+\tfrac1k\right)\frac{\sigma^2(z_{1-\alpha}+z_{1-\beta})^2}{(\Delta+\delta)^2}$$"
        en = f"Using a non-inferiority margin of {margin:g} (SD {sd:g}, expected difference {tdiff:g}), a one-sided α of {alpha} and {pct(power)} power"
        tr = f"{margin:g} eşdeğer-olmama (non-inferiority) sınırı (SS {sd:g}, beklenen fark {tdiff:g}), tek yönlü α={alpha} ve %{power*100:g} güçle"
    elif kind.startswith("Eşdeğerlik"):
        r = calc.equivalence_means(sd, margin, tdiff, alpha, power, ratio)
        f = r"$$n_1=\left(1+\tfrac1k\right)\frac{\sigma^2(z_{1-\alpha}+z_{1-\beta/2})^2}{\delta^2}\quad(\Delta=0)$$" "\n\n$\\Delta\\neq0$ ise $z_{1-\\beta}$ ve $(\\delta-|\\Delta|)^2$ kullanılır (yaklaşık) / if $\\Delta\\neq0$, $z_{1-\\beta}$ and $(\\delta-|\\Delta|)^2$ are used (approximate)."
        en = f"Using equivalence margins of ±{margin:g} (SD {sd:g}, expected difference {tdiff:g}), two one-sided tests at α = {alpha} and {pct(power)} power"
        tr = f"±{margin:g} eşdeğerlik sınırları (SS {sd:g}, beklenen fark {tdiff:g}), α={alpha} düzeyinde iki tek yönlü test (TOST) ve %{power*100:g} güçle"
    else:
        r = calc.ni_proportions(pc, pt, margin, alpha, power, ratio)
        f = r"$$n_1=\frac{(z_{1-\alpha}+z_{1-\beta})^2\left[p_c(1-p_c)+\tfrac{p_t(1-p_t)}{k}\right]}{(p_t-p_c+\delta)^2}$$"
        en = f"Assuming success rates of {pct(pc)} (control) and {pct(pt)} (test), a non-inferiority margin of {margin*100:g} percentage points, a one-sided α of {alpha} and {pct(power)} power"
        tr = f"Kontrolde %{pc*100:g}, test grubunda %{pt*100:g} başarı oranı, {margin*100:g} yüzde puanlık non-inferiority sınırı, tek yönlü α={alpha} ve %{power*100:g} güçle"
    adj, extra = drop_rows(r["total"], dropout)
    return report(
        [("Grup 1 (kontrol) / Group 1 (control)", r["n1"]),
         ("Grup 2 (test) / Group 2 (test)", r["n2"]), ("Toplam / Total", r["total"])] + extra,
        f,
        "- Yüksek değer = daha iyi / higher = better\n- Normal yaklaşım / normal approximation\n- Sınır (δ) klinik gerekçeyle önceden belirlenmeli / the margin must be pre-specified on clinical grounds",
        f"{en}, {r['n1']} and {r['n2']} participants are required ({r['total']} in total); allowing for {pct(dropout)} attrition, {adj} participants will be enrolled.",
        f"{tr} gruplarda {r['n1']} ve {r['n2']} (toplam {r['total']}) katılımcı gerekmektedir; %{dropout*100:g} kayıp öngörülerek {adj} katılımcı alınacaktır.",
        ["chow"])


# ---------------------------------------------------------------- UI
def common(two_sided=True, alpha=0.05):
    a = gr.Number(alpha, label="α (anlamlılık / significance)", minimum=0.0001, maximum=0.5)
    p = gr.Slider(0.5, 0.99, 0.80, step=0.01, label="Güç / Power (1−β)")
    d = gr.Slider(0, 0.6, 0.10, step=0.01, label="Kayıp / Dropout")
    return a, p, d


def sides_radio():
    return gr.Radio(["2", "1"], value="2", label="Yön / Sides")


def out():
    return gr.Markdown(latex_delimiters=LATEX)


INTRO = """# 🧮 Klinik Örneklem Büyüklüğü · Clinical Sample Size
Etik kurul başvuruları ve makale yöntem bölümleri için a priori örneklem hesabı. Her sonuç **formül, varsayımlar, kopyalanabilir yöntem cümlesi (EN/TR) ve kaynakla** birlikte verilir.
A priori sample size for ethics applications and methods sections — each result comes with the **formula, assumptions, a ready-to-paste methods sentence (EN/TR) and the reference**.

<small>Faruk Danış, MD · Emergency Medicine · [ORCID 0000-0001-8722-5402](https://orcid.org/0000-0001-8722-5402) · Eğitim ve planlama amaçlıdır; nihai hesabı bir biyoistatistikçiyle doğrulayın. / For education and planning; confirm final calculations with a biostatistician.</small>
"""

with gr.Blocks(title="Clinical Sample Size Calculator") as demo:
    gr.Markdown(INTRO)
    with gr.Tabs():
        with gr.Tab("İki oran / Two proportions"):
            with gr.Row():
                p1 = gr.Number(0.40, label="p1 (grup 1 oranı / group 1)")
                p2 = gr.Number(0.60, label="p2 (grup 2 oranı / group 2)")
                k = gr.Number(1.0, label="Oran n2/n1 / Ratio n2/n1", minimum=0.1)
            with gr.Row():
                a, pw, dr = common(); s = sides_radio()
                cc = gr.Checkbox(False, label="Fleiss süreklilik düzeltmesi / continuity correction")
            o = out()
            gr.Button("Hesapla / Calculate", variant="primary").click(h_props, [p1, p2, a, pw, k, s, cc, dr], o)

        with gr.Tab("İki ortalama / Two means"):
            mode = gr.Radio(["Fark + SS / Difference + SD", "Cohen d"], value="Fark + SS / Difference + SD", label="Girdi / Input")
            with gr.Row():
                diff = gr.Number(5, label="Ortalama farkı / Mean difference")
                sd = gr.Number(10, label="Ortak SS / Common SD")
                d = gr.Number(0.5, label="Cohen d")
                k = gr.Number(1.0, label="Oran n2/n1 / Ratio", minimum=0.1)
            with gr.Row():
                a, pw, dr = common(); s = sides_radio()
            o = out()
            gr.Button("Hesapla / Calculate", variant="primary").click(h_means, [mode, diff, sd, d, a, pw, k, s, dr], o)

        with gr.Tab("Tanısal doğruluk / Diagnostic accuracy"):
            with gr.Row():
                se = gr.Number(0.90, label="Beklenen duyarlılık / Expected sensitivity")
                sp = gr.Number(0.85, label="Beklenen özgüllük / Expected specificity")
                pv = gr.Number(0.20, label="Prevalans / Prevalence")
                pr = gr.Number(0.05, label="Kesinlik (GA yarı genişliği) / Precision (CI half-width)")
            with gr.Row():
                cf = gr.Radio([0.90, 0.95, 0.99], value=0.95, label="Güven düzeyi / Confidence level")
                dr = gr.Slider(0, 0.6, 0.10, step=0.01, label="Dışlanma / Exclusions")
            o = out()
            gr.Button("Hesapla / Calculate", variant="primary").click(h_diag, [se, sp, pv, pr, cf, dr], o)

        with gr.Tab("ROC AUC"):
            with gr.Row():
                au = gr.Number(0.70, label="Beklenen AUC / Expected AUC")
                au0 = gr.Number(0.50, label="H0 AUC")
                k = gr.Number(1.0, label="Negatif/pozitif oranı / Negatives per positive", minimum=0.1)
            with gr.Row():
                a, pw, dr = common(); s = sides_radio()
            o = out()
            gr.Button("Hesapla / Calculate", variant="primary").click(h_auc, [au, au0, a, pw, k, s, dr], o)

        with gr.Tab("Lojistik regresyon / Logistic regression"):
            with gr.Row():
                kk = gr.Number(10, label="Aday parametre sayısı / Candidate parameters", precision=0)
                ph = gr.Number(0.20, label="Sonuç oranı / Outcome proportion")
                ep = gr.Number(10, label="EPV")
                r2 = gr.Number(0.0, label="Cox–Snell R² (Riley; 0 = kullanma / skip)")
            dr = gr.Slider(0, 0.6, 0.10, step=0.01, label="Eksik veri / Missing data")
            o = out()
            gr.Button("Hesapla / Calculate", variant="primary").click(h_logit, [kk, ph, ep, r2, dr], o)

        with gr.Tab("Korelasyon / Correlation"):
            with gr.Row():
                rr = gr.Number(0.30, label="Beklenen r / Expected r")
                a, pw, dr = common(); s = sides_radio()
            o = out()
            gr.Button("Hesapla / Calculate", variant="primary").click(h_corr, [rr, a, pw, s, dr], o)

        with gr.Tab("Sağkalım / Survival (log-rank)"):
            with gr.Row():
                hr = gr.Number(0.70, label="Hazard oranı / Hazard ratio")
                pe = gr.Number(0.50, label="Genel olay olasılığı / Overall event probability")
                al = gr.Number(0.5, label="Grup 1 payı / Allocation to group 1")
                me = gr.Radio(["Schoenfeld (1983)", "Freedman (1982)"], value="Schoenfeld (1983)", label="Yöntem / Method")
            with gr.Row():
                a, pw, dr = common(); s = sides_radio()
            o = out()
            gr.Button("Hesapla / Calculate", variant="primary").click(h_surv, [hr, pe, a, pw, al, s, me, dr], o)

        with gr.Tab("Eşdeğerlik · NI / Equivalence · NI"):
            kind = gr.Radio(["NI — ortalama / means", "Eşdeğerlik — ortalama / Equivalence — means", "NI — oran / proportions"],
                            value="NI — ortalama / means", label="Tasarım / Design")
            with gr.Row():
                sd = gr.Number(1.0, label="SS / SD (ortalama / means)")
                mg = gr.Number(0.5, label="Sınır δ / Margin δ")
                td = gr.Number(0.0, label="Beklenen fark Δ / Expected difference Δ")
            with gr.Row():
                pc = gr.Number(0.80, label="Kontrol oranı / Control proportion")
                pt = gr.Number(0.80, label="Test oranı / Test proportion")
                k = gr.Number(1.0, label="Oran n2/n1 / Ratio", minimum=0.1)
            with gr.Row():
                a, pw, dr = common(alpha=0.025)
            gr.Markdown("<small>NI için α tek yönlüdür (genellikle 0,025); eşdeğerlikte TOST α'sı (genellikle 0,05). / One-sided α for NI (usually 0.025); TOST α for equivalence (usually 0.05).</small>")
            o = out()
            gr.Button("Hesapla / Calculate", variant="primary").click(h_ni, [kind, sd, mg, td, pc, pt, a, pw, k, dr], o)

    gr.Markdown("<small>Kaynak kod / Source: [GitHub](https://github.com/farukdanis/clinical-sample-size) · MIT License</small>")

if __name__ == "__main__":
    demo.launch()
