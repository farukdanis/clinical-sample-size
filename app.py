"""Clinical Sample Size Calculator — Gradio app (TR/EN).

Author: Faruk Danış, MD (ORCID 0000-0001-8722-5402)
"""
import gradio as gr

from reports import (LATEX, h_auc, h_corr, h_diag, h_logit, h_means, h_ni,
                     h_props, h_surv)


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
