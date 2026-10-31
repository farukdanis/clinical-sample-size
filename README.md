---
title: Clinical Sample Size Calculator
emoji: 🧮
colorFrom: blue
colorTo: indigo
sdk: static
app_file: index.html
pinned: true
license: mit
short_description: A priori sample size for clinical studies (TR/EN)
tags:
  - medical
  - statistics
  - sample-size
  - clinical-research
  - emergency-medicine
---

# 🧮 Clinical Sample Size Calculator · Klinik Örneklem Büyüklüğü

A bilingual (Turkish / English) a priori sample size calculator for clinical researchers, ethics committee applications and manuscript methods sections.
Every result is returned with the **formula, assumptions, a ready-to-paste methods sentence in English and Turkish, and the original reference**.

Etik kurul başvuruları ve makale yöntem bölümleri için iki dilli a priori örneklem hesaplayıcı. Her sonuç **formül, varsayımlar, kopyalanabilir İngilizce/Türkçe yöntem cümlesi ve özgün kaynakla** verilir.

## Modules / Modüller

| Module | Method | Reference |
|---|---|---|
| Two proportions | Normal approximation, pooled variance under H0, optional Fleiss continuity correction, unequal allocation | Fleiss, Levin & Paik 2003 |
| Two means | Exact noncentral-t solution (matches G\*Power and statsmodels) | Cohen 1988 |
| Diagnostic accuracy | Precision-based sensitivity/specificity, prevalence-adjusted | Buderer 1996 |
| ROC AUC | Test of AUC vs AUC0, Hanley–McNeil variance, numeric search | Hanley & McNeil 1982 |
| Logistic regression | Events-per-variable rule + Riley criteria (i) shrinkage and (iii) intercept precision | Peduzzi 1996; Riley 2020 |
| Correlation | Fisher z transformation | Hulley et al. 2013 |
| Survival (log-rank) | Required events by Schoenfeld or Freedman, converted to N by event probability | Schoenfeld 1983; Freedman 1982 |
| Non-inferiority / equivalence | Means (NI, TOST) and proportions (NI), normal approximation | Chow, Shao & Wang 2008 |

All modules include a dropout / exclusion adjustment: $N_{adj} = N / (1 - d)$.

## Validation / Doğrulama

`tests/test_calc.py` checks the formulas against published or widely used reference values:

- Two proportions 0.40 vs 0.60 → 97 per group (107 with continuity correction); Monte Carlo power check for unequal allocation
- Two means, d = 0.5 → 64 per group; d = 0.8 → 26 per group (G\*Power); identical to statsmodels to 0.001 for unequal allocation and one-sided tests
- Buderer: Se 0.90, precision 0.05, prevalence 0.20 → 139 diseased, 692 total
- Correlation r = 0.30 → 85 (Hulley et al.)
- Schoenfeld HR 0.70 → 247 events; Freedman → 253 events
- ROC: minimal-n property and monotonicity; NI/TOST and Riley formulas; dropout and input validation

## Two ways to run / İki çalışma biçimi

- **In the browser (Hugging Face Space):** `index.html` loads [Pyodide](https://pyodide.org) and runs the same `calc.py` + `reports.py` locally in the browser. No server, no data leaves the device.
- **Locally with Gradio:** `app.py` wraps the same functions in a Gradio interface.

```bash
pip install -r requirements.txt -r requirements-dev.txt
pytest -q          # validation tests (statsmodels used only as a cross-check)
python app.py      # Gradio UI
python -m http.server  # then open http://localhost:8000 for the static UI
```

## References / Kaynaklar

1. Fleiss JL, Levin B, Paik MC. *Statistical Methods for Rates and Proportions*. 3rd ed. Hoboken (NJ): Wiley; 2003.
2. Cohen J. *Statistical Power Analysis for the Behavioral Sciences*. 2nd ed. Hillsdale (NJ): Lawrence Erlbaum; 1988.
3. Buderer NM. Statistical methodology: I. Incorporating the prevalence of disease into the sample size calculation for sensitivity and specificity. *Acad Emerg Med*. 1996;3(9):895-900.
4. Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) curve. *Radiology*. 1982;143(1):29-36.
5. Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR. A simulation study of the number of events per variable in logistic regression analysis. *J Clin Epidemiol*. 1996;49(12):1373-9.
6. Riley RD, Ensor J, Snell KIE, et al. Calculating the sample size required for developing a clinical prediction model. *BMJ*. 2020;368:m441.
7. Hulley SB, Cummings SR, Browner WS, Grady DG, Newman TB. *Designing Clinical Research*. 4th ed. Philadelphia: Lippincott Williams & Wilkins; 2013.
8. Schoenfeld DA. Sample-size formula for the proportional-hazards regression model. *Biometrics*. 1983;39(2):499-503.
9. Freedman LS. Tables of the number of patients required in clinical trials using the logrank test. *Stat Med*. 1982;1(2):121-9.
10. Chow SC, Shao J, Wang H. *Sample Size Calculations in Clinical Research*. 2nd ed. Boca Raton (FL): Chapman & Hall/CRC; 2008.

## Disclaimer / Uyarı

For education and study planning. Confirm final calculations for regulated trials with a biostatistician.
Eğitim ve çalışma planlaması amaçlıdır; düzenlemeye tabi çalışmalarda nihai hesabı bir biyoistatistikçiyle doğrulayın.

## Author / Yazar

**Faruk Danış, MD** — Assistant Professor of Emergency Medicine, Bolu Abant İzzet Baysal University · [ORCID 0000-0001-8722-5402](https://orcid.org/0000-0001-8722-5402) · [GitHub](https://github.com/farukdanis) · [faruk.md](https://www.faruk.md)

License: MIT
