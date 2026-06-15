---
title: "ADVANCE Network-Diffusion Study"
subtitle: "Peer exposure, threshold, susceptibility and infectiousness in adolescent e-cigarette adoption (W1–W10)"
author: "Y. Cao, T. Valente, K. Miljkovic, A. Olivera"
date: \today
geometry: "margin=2.5cm"
fontsize: 10pt
colorlinks: true
linkcolor: teal
toc: true
toc-depth: 3
header-includes:
  - \usepackage{amsmath}
  - \usepackage{amssymb}
  - \usepackage{booktabs}
  - \usepackage{array}
  - \usepackage{newunicodechar}
  - \newunicodechar{≈}{\ensuremath{\approx}}
  - \newunicodechar{≤}{\ensuremath{\le}}
  - \newunicodechar{≥}{\ensuremath{\ge}}
  - \newunicodechar{≠}{\ensuremath{\neq}}
  - \newunicodechar{×}{\ensuremath{\times}}
  - \newunicodechar{→}{\ensuremath{\rightarrow}}
  - \newunicodechar{…}{\ensuremath{\ldots}}
  - \newunicodechar{√}{\ensuremath{\surd}}
  - \newunicodechar{⊔}{\ensuremath{\sqcup}}
---

# 1. Introduction

This is the threshold/diffusion companion to the ADVANCE *disadoption* study. Where the
disadoption paper models the **1 → 0** transition (quitting), this project models the
**0 → 1** transition (adoption) on the same USC ADVANCE panel. The lead analyst is Y. Cao; the
current pipeline is [`codes/260612_netdiffuseR_past6mo_revised13.R`](../codes/260612_netdiffuseR_past6mo_revised13.R)
on `main`.

**Scope (2026-06-15).** The study is now focused on the regressions reported in
`docs/260612_threshold_report_final.pdf`, namely:

1. **Adoption** — a mixed-effects discrete-time hazard (random intercepts for student and school),
   with a school-fixed-effects logistic regression as a sensitivity check (report Tables 2–3).
2. **Threshold** — the degree-weighted grouped-binomial model (report Tables 2–3).
3. **Peer-independent subgroups** — logistic "membership" models for users (No Exposure, No
   Perceived Friends, Neither, Low Threshold) and non-users (High Exposure, High Perceived Friends,
   Both High), including a social-media exposure term (report Tables 1a–1b).

The susceptibility and infectiousness regressions of the earlier `260603` pipeline have been
**retired** from the current report; their specifications are kept in §3 for the record but no
longer carry results. The substantive question (per T. Valente, 2026-06-03) is the
**peer-independent adopter** — students who adopt with no network exposure, no perceived friends,
or a low threshold — and whether **social-media exposure** explains them.

All regression cells report **odds ratios** $\exp(\hat\beta)$ with the model p-value. Numbers in
this document were reproduced end-to-end on the 042326 release with **netdiffuseR 1.25.0**
(see §6). Two corrected, AO-labelled re-analyses accompany the report
([`codes/07-threshold-corrected-AO.R`](../codes/07-threshold-corrected-AO.R),
[`codes/08-socialmedia-power-AO.R`](../codes/08-socialmedia-power-AO.R)); their results are in §5.

# 2. Data

## 2.1 Sources and sample

- **W1–W8**: `ADVANCE_W1-W8_Data_Complete_042326.xlsx` — 4,437 students × 10,341 columns.
- **W9–W10 HS**: `ADVANCE_W9-W10_HS_Data_Complete_042326.xlsx` — 1,060 students (subset of W1–W8).
- These are the **same physical files** the disadoption study uses
  (`A - Network-Disadoption/data/advance/Data/`). The new `260603` pipeline reads the
  consolidated XLSX directly (slicing per-wave columns `w{w}_*`); the older `260427` pipeline
  read per-wave `wN_adv_data.csv`. The friendship edges are parsed in-script from the
  `w{w}_friend{k}_{school}` nomination columns.

## 2.2 Cohorts, schools and the grade-period timeline

Fourteen schools in two cohorts. School 108 is excluded; schools with $<20$ students are
dropped. Ethnic composition (`sch_type`) is assigned by school: Asian-majority
{103,105,112,113,212,213}, Hispanic-majority {102,106,107,114,214}, Other otherwise.

| Cohort | Schools | Waves used | Grade-period (GP) map |
|:--|:--|:--|:--|
| C1 early | 101–105 | W3–W8 | GP = W |
| C1 late  | 106,107,112,113,114 | W3–W8 | GP = W |
| C2       | 201,212,213,214 | W5–W10 | GP = W − 2 |

Both cohorts are aligned on **six grade-periods**, GP3–GP8 = {10th Fall, 10th Spring, 11th Fall,
11th Spring, 12th Fall, 12th Spring}. Grade 9 (GP1–GP2) is dropped so that C1-late and C2 schools,
which only enter at 10th-Fall, are comparable. Combined panel: **4,007 students, 841 adopters
(21.0% by 12th-Spring)**.

## 2.3 Diffusion network construction

For each school $s$ a directed friendship network is built per grade-period from the nomination
columns; ties are kept only between consenting panel members. The per-school diffnet is
assembled with `as_diffnet(graph, toa, …)` and the 14 are merged with `do.call(c, …)` into the
combined object $G$ on which all network statistics are computed. Past-6-month use is
**monotonised** (once 1, stays 1), and the time-of-adoption is
$\tau_i = \min\{t : y_{it} = 1\}$.

## 2.4 Covariates

Demographics are taken at each student's **baseline wave** (first wave their school appears):
`female = 1[dem_gender = 0]`, `hispanic = 1[eth = 1]`, `asian = 1[race = 2]`,
`sex_min = 1[dem_sexuality ≠ 1]`. Time-varying covariates merged at the relevant grade-period:
`gad = rcads_gad_mean`, `mdd = rcads_mdd_mean`, `friends_ecig = friends_use_ecig` (perceived
friend use, "6" recoded NA). Social-media exposure (`sm_post_1..6`, `ecig_posted_1..3`) exists
**only W4–W8** — a fact that dominates §6.

### 2.4.1 Parent education — harmonisation, and how it differs from the disadoption study

Parent education needs cross-wave handling because the questionnaire scale changed mid-panel.
Per the codebook (042326):

- **Legacy** `w{w}_dem_high_par_edu`: `1`=≤8th grade, `2`=some HS, `3`=HS grad,
  `4`=some college, `5`=college grad, `6`=advanced degree, **`9`=Don't know**.
- **New** `w{w}_dem_high_par_edu_new` (W7–W10 only): `1..4` as above, `5`=vocational,
  `6`=associate, `7`=bachelor's, `8`=master's, `9`=Don't know.

The `260603` pipeline reads the **legacy-named** column at all waves, remaps W7+ with
`c(1,2,3,4,4,4,5,6,NA)`, and then assigns each student their **first valid value across all
waves** (a single, effectively time-invariant number).

**This is *not* equivalent to the disadoption construction**, in three ways:

| | Disadoption (`R/01-advance-panel.R`) | Diffusion (`260603`) |
|:--|:--|:--|
| Fill rule | **LOCF** per student-wave (time-varying) | **first-valid** across waves (one static value) |
| "Don't know" (9) at W7+ | $9 \to 7$ (kept as top level) | $9 \to$ **NA** (dropped) |
| "Don't know" (9) at W1–W6 | kept as **9** | kept as **9** |

Two consequences worth flagging:

1. **A latent scale bug, shared by both.** Empirically, in the 042326 release the legacy-named
   column `w{w}_dem_high_par_edu` is on the **legacy 6-level scale at every wave** (values
   $\{1,\dots,6,9\}$; no 7s or 8s, even at W7–W8 — the new scale lives in the tiny
   `…_new` column, $n=98$ at W7). So the W7+ "new→legacy" remap is applied to data that is
   *already* legacy-scaled: its `5→4, 6→4` rules **mis-collapse "college graduate" and "advanced
   degree" into "some college"** for any value drawn from W7–W10. Impact is bounded by which wave
   supplies each student's value (most come from an early legacy wave), but it is a real defect in
   both pipelines.
2. **"Don't know" leaks as a numeric level in this pipeline.** Because `9 → NA` only fires for
   $w \ge 7$, and *first-valid* almost always pulls from an early ($w<7$) wave, the value **9
   survives** and enters the GLMs as a continuous predictor ranked **above** "advanced degree" (6).
   In the threshold sample, **60 of 835 adopters (7.2%) have `par_edu = 9`**. This contaminates any
   model that treats `par_edu` as continuous — and `par_edu` is "significant" in exactly the two
   models where it should be least expected (susceptibility OR 1.18, $p=.012$; infectiousness
   OR 1.25, $p=.016$), which is a red flag rather than a finding.

**Recommended fix:** recode `9 → NA` at *all* waves; since the read column is legacy-scaled
throughout this release, drop the W7+ new-scale remap (or read `…_new` only where populated);
and treat `par_edu` as an ordered factor (or 1–6 continuous with DK = NA).

## 2.5 Temporal alignment — are the regressions lagged?

Short answer: **only the adoption model is lagged.** The four model families differ:

| Model | Predictor timing | Outcome timing | Lagged? |
|:--|:--|:--|:--|
| **Adoption / exposure** | period $t$ | adoption at $t+1$ | **Yes** (one period) |
| **Threshold** | at TOA $\tau_i$ | at TOA $\tau_i$ | No (contemporaneous) |
| **Susceptibility** | at TOA $\tau_i$ | at TOA $\tau_i$ | No (contemporaneous) |
| **Infectiousness** | at TOA $\tau_i$ | at TOA $\tau_i$ | No (contemporaneous) |

The adoption model is a discrete-time hazard: exposure and covariates measured at $t$ predict
*next-period* adoption (`adopt_next = adopt_mat[, t+1]`), so the predictor strictly precedes the
event. The threshold/susceptibility/infectiousness models put **one row per adopter at their
adoption period** and merge `gad`, `mdd`, `friends_ecig`, demographics **at $\tau_i$ itself**.
The *network* quantity is inherently ordered (threshold = friends already using *before* TOA;
susceptibility = friends who adopted in $\tau_i-1$; infectiousness = friends adopting in
$\tau_i+1$), but the **explanatory covariates are concurrent with adoption**, which raises a
simultaneity caveat: e.g. mental-health state *at* the moment of adoption may be a consequence,
not a cause. In the disadoption study predictors are taken at $w-1$; matching that convention
here (covariates at $\tau_i - 1$) would make the two studies consistent and reduce the
reverse-causation risk.

# 3. Methods — model specifications

Let $G = (V, E)$ be the combined directed diffnet, $t \in \{3,\dots,8\}$ the grade-periods,
$y_{it} \in \{0,1\}$ monotone past-6-month use, $\tau_i = \min\{t: y_{it}=1\}$ the time of
adoption, $A_i^{\text{out}}(t) = \{j : (i \to j) \in E_t\}$ ego $i$'s out-neighbours, and
$E_{it} = \frac{\sum_{j \in A_i^{\text{out}}(t)} y_{jt}}{|A_i^{\text{out}}(t)|}$ the network
exposure (share of alters already adopted). $\mathbf{x}_i$ is the covariate vector
(cohort, female, hispanic, asian, par_edu, gad, mdd, friends_ecig, sex_min), $\alpha_t$ a
grade-period fixed effect, and $\xi_{s}$ a school fixed effect (C2 schools collapsed to the
reference so `factor(school_c1fe)` yields only C1 dummies, with cohort absorbing the C1-vs-C2
mean).

## 3.1 Threshold — grouped binomial

For each non-isolate adopter ($n_i = |A_i^{\text{out}}(\tau_i)| > 0$), let
$k_i = |\{ j \in A_i^{\text{out}}(\tau_i) : \tau_j < \tau_i \}|$ be the number of out-alters who
adopted **strictly before** $i$. The empirical threshold is $\hat T_i = k_i / n_i$. We model the
counts directly:
$$
k_i \mid n_i \;\sim\; \mathrm{Binomial}(n_i,\, \pi_i),
\qquad
\operatorname{logit}(\pi_i) \;=\; \alpha_{\tau_i} + \xi_{s(i)} + \mathbf{x}_i^{\top}\boldsymbol\beta .
$$
In R: `glm(cbind(k_users, n_alters - k_users) ~ factor(toa) + cohort + female + hispanic +
asian + par_edu + factor(school_c1fe) + gad + mdd + friends_ecig + sex_min, family = binomial)`.
Isolates are excluded because $\hat T_i = 0/0$ is undefined. This **degree-weighted** binomial is
the re-specification adopted from the disadoption collaborator's 2026-05-12 diagnostic: it
replaces a quasibinomial on the raw proportion $\hat T_i$, which over-weights low-degree adopters
whose threshold lives on a coarse rational grid $\{0, \tfrac12, \tfrac13, \dots\}$.

## 3.2 Adoption — discrete-time exposure hazard

Person-period rows $(i,t)$ for every still-at-risk student ($y_{it}=0$); predictors at $t$,
outcome at $t+1$. The **current report's main model (Table 2)** is a *mixed-effects* logistic
hazard with random intercepts for student and school:
$$
\Pr\!\big(y_{i,t+1} = 1 \mid y_{it} = 0\big)
= \operatorname{logit}^{-1}\!\big(\gamma\, E_{it} + \alpha_t + \mathbf{x}_i^{\top}\boldsymbol\beta
  + u_{i} + v_{s(i)}\big),
\qquad u_i \sim \mathcal N(0,\sigma_u^2),\; v_s \sim \mathcal N(0,\sigma_v^2).
$$
`lme4::glmer(adopt_next ~ exposure + factor(grade_period) + … + gad + mdd + friends_ecig +
sex_min + (1 | id) + (1 | school), family = binomial)`. The **sensitivity model (Table 3)**
replaces the random effects with school fixed effects so individual school ORs are inspectable:
$$
\Pr\!\big(y_{i,t+1} = 1 \mid y_{it} = 0\big)
= \operatorname{logit}^{-1}\!\big(\gamma\, E_{it} + \alpha_t + \xi_{s(i)} + \mathbf{x}_i^{\top}\boldsymbol\beta\big).
$$
The mixed model returns *subject-specific* ORs (larger in magnitude, e.g. exposure 6.46) while the
fixed-effects GLM returns *population-average* ORs (exposure 2.31); both agree on sign and
significance for every predictor.

## 3.3 Susceptibility and infectiousness — quasibinomial *(retired from current scope)*

These were computed in the earlier `260603` pipeline but are **not part of the current report**;
the specification is kept here for the record only. One row per adopter; **susceptibility** is the
share of $i$'s out-alters who adopted *immediately before* $i$ (excludes $\tau_i = 3$),
**infectiousness** the share adopting *immediately after* (excludes $\tau_i = 8$), each a
proportion in $[0,1]$ modelled with an over-dispersed logit
$\operatorname{logit}\mathbb{E}[S_i] = \alpha_{\tau_i} + \xi_{s(i)} + \mathbf{x}_i^{\top}\boldsymbol\beta$,
$\operatorname{Var}(S_i) = \phi\,\mu_i(1-\mu_i)$ (`family = quasibinomial`).

## 3.4 Subgroup membership — "Tom's table"

Seven subgroups — four among **users** (report Table 1a), three among **non-users** (Table 1b) —
each modelled as binary membership within its stratum:
$$
\operatorname{logit}\Pr(m_{ic} = 1) = \alpha_{\tau_i} + \mathbf{z}_i^{\top}\boldsymbol\theta_c ,
$$
where $\mathbf z_i$ stacks demographics, GP3/GP8 degree, perceived friend use, exposure, and a
social-media term; $\alpha_{\tau_i}$ are grade-period fixed effects. The predictor that *defines*
a column is excluded from that column's model (definitional overlap). **Users:** *No Exposure*
$E_{i\tau_i}=0$; *No Perceived Friends* PFU$_{\tau_i}=0$; *Neither* both $=0$; *Low Threshold*
$\hat T_i \le 0.25$ (the sample's 70th percentile, per Valente 1996). **Non-users:** *High
Exposure* $\max_t E_{it}>0$; *High Perceived Friends* $\max_t \text{PFU}_{it}>0$; *Both High* both.

## 3.5 The social-media term

In the current report the social-media exposure indicator $\mathrm{SM}_i$ (`sm_ecig_any` at TOA =
saw e-cig content on any of six platforms, W4–W8) enters the **user subgroup model** of §3.4 as one
more covariate, and its OR $\exp(\hat\eta)$ is the `sm_ecig_any` row of Table 1a:
$$
\operatorname{logit}\Pr(m_{ic} = 1) = \eta\,\mathrm{SM}_i + \alpha_{\tau_i} + \mathbf{z}_i^{\top}\boldsymbol\theta .
$$
(An older standalone 18-model sensitivity sweep — six SM operationalisations × three subgroups —
also exists in the pipeline but is not the reported result.)

# 4. Results — what we obtained re-running the pipeline

These are **our** numbers from re-running [`codes/260612_netdiffuseR_past6mo_revised13.R`](../codes/260612_netdiffuseR_past6mo_revised13.R)
on the 042326 release (netdiffuseR 1.25.0), laid out in the **same tables** as
`260612_threshold_report_final.pdf`. Each table is **not one regression but several** — we state
the model behind each column, then interpret it in plain words. ** = $p<.01$, * = $p<.05$.

## 4.1 Tables 2–3 — Adoption and Threshold side by side

Each of Tables 2 and 3 juxtaposes **two different regressions**: an *adoption* model (does a still-
non-using student start next period?) and the *threshold* model (among those who start, how built-up
was the peer pressure?). Table 2 uses the **mixed-effects** adoption model, Table 3 the
**school-fixed-effects** version; the threshold column is the same grouped binomial in both.

**The regressions (math).** Adoption is a discrete-time hazard over person-periods $(i,t)$ with
$E_{it}$ = share of $i$'s friends already using:
$$
\textbf{Table 2 (mixed): } \operatorname{logit}\Pr(y_{i,t+1}{=}1\mid y_{it}{=}0)=\gamma E_{it}+\alpha_t+\mathbf{x}_i^\top\boldsymbol\beta+u_i+v_{s(i)},
$$
$$
\textbf{Table 3 (FE): } \operatorname{logit}\Pr(y_{i,t+1}{=}1\mid y_{it}{=}0)=\gamma E_{it}+\alpha_t+\xi_{s(i)}+\mathbf{x}_i^\top\boldsymbol\beta .
$$
Threshold is the degree-weighted grouped binomial over adopters ($n_i$ friends, $k_i$ already using):
$$
\textbf{Threshold: } k_i\mid n_i\sim\mathrm{Binomial}(n_i,\pi_i),\quad \operatorname{logit}\pi_i=\alpha_{\tau_i}+\xi_{s(i)}+\mathbf{x}_i^\top\boldsymbol\beta .
$$

**Our numbers** (Adoption ORs; Threshold ORs, n = 406 non-isolate adopters):

| Predictor | Adoption — Table 2 mixed | Adoption — Table 3 FE | Threshold (Tables 2–3) |
|:--|--:|--:|--:|
| Network exposure $\exp(\gamma)$ | 7.6 \*\* | **2.31** \*\* | — |
| Perceived friend use | 1.72 \*\* | **1.43** \*\* | **1.20** \*\* |
| MDD | (unstable) | **1.61** \*\* | 0.81 |
| GAD | (unstable) | **0.77** \* | 0.97 |
| Cohort C2 | (unstable) | 1.10 | **2.64** \* |
| Asian | (unstable) | 0.77 | 0.60 (.053) |
| Hispanic / female / par_edu / sex_min | (unstable) | ns | ns |
| School 104 / 114 (Hispanic-maj.) | — | ns | **2.99\* / 3.63\*\*** |
| Constant; GP4 … GP8 | GP8 separates | GP8 separates | **−20.6; 17.2 … 18.8 (degenerate)** |

Two reproduction notes. (1) The **FE adoption model (Table 3) reproduces the report exactly** —
exposure 2.31, MDD 1.61, GAD 0.77, PFU 1.43. (2) The **mixed model (Table 2) is numerically
unstable**: on identical data our `glmer` converged to a different point (the demographic ORs blow
up and GP8 separates), so we trust Table 3 and treat the mixed ORs as not robustly identified — a
first sign that the sparse, separated design strains these models. The Threshold column's
**Constant and all GP rows are a separation artifact** (the `−20.6 / 17.2 … 18.8` block in the
report) — see §5.1.

**For grandma.** *Whether a teen starts vaping in the next half-year depends most on how many of
their friends already vape — a kid surrounded by vapers is about **2.3× more likely** to start.
Feeling depressed nudges them toward it; feeling anxious, slightly away from it; thinking your
friends vape adds a push too. Now, among the kids who do start, the "threshold" asks how many
friends were already vaping at that moment: kids who start in 12th grade waited until far more
friends were already vaping, while early-starters jumped in with almost no one — and believing your
friends vape goes with waiting for a higher bar. (The giant numbers in the grade rows are a computer
hiccup, not a real result; we fix it in §5.)*

## 4.2 Tables 1a–1b — Peer-independent subgroups (each column is its own regression)

Every column of Tables 1a/1b is a **separate logistic regression**: among users (1a) or non-users
(1b), it asks which traits predict belonging to that "peer-independent" profile, controlling for
grade-period:
$$
\operatorname{logit}\Pr(\text{member of subgroup } c)=\alpha_{\tau_i}+\mathbf{z}_i^\top\boldsymbol\theta_c .
$$
The **social-media term** `sm_ecig_any` sits inside $\mathbf z$ for the user columns; it is undefined
for non-users (no time-of-adoption) and shown as "—".

**Table 1a — users (our ORs; ** $p{<}.01$, * $p{<}.05$):**

| predictor | No Exposure | No Perc. Friends | Neither | Low Threshold |
|:--|--:|--:|--:|--:|
| Asian | 2.49 (.05) | 2.31 | **3.92** \* | **3.53** \*\* |
| MDD | 1.07 | 0.53 | **0.46** \* | 0.79 |
| Out-degree (GP3) | **0.81** \* | 1.10 | 0.94 | 0.98 |
| Perceived friends (GP3) | **0.74** \* | 0.81 | **0.51** \*\* | 0.86 |
| Perceived friends (GP8) | 0.85 | **0.56** \*\* | **0.57** \*\* | 0.85 |
| **`sm_ecig_any` (W4–W8)** | **2.31** \* | **0.40** \* | 0.87 | 1.39 |

**Table 1b — non-users / "resisters" (our ORs):**

| predictor | High Exposure | High Perc. Friends | Both High |
|:--|--:|--:|--:|
| Asian | **0.51** \*\* | **0.34** \*\* | **0.37** \*\* |
| Female | 1.03 | **1.30** \*\* | **1.63** \*\* |
| MDD | 1.20 | **1.51** \*\* | **1.53** \*\* |
| Out-degree (GP3) | **1.39** \*\* | **1.06** \* | **1.25** \*\* |
| Perceived friends (GP8) | **1.47** \*\* | — | **3.42** \*\* |
| par_edu | 0.96 | **0.93** \*\* | 0.98 |

These reproduce the report exactly. Note `par_edu` looks significant in the "High Perceived Friends"
column (0.93\*\*) — a red flag we trace to a coding bug and fix in §5.3.

**For grandma.** *We sorted the teens into "types". Among vapers who had **no vaping friends** in
the network, Asian kids are over-represented and — surprisingly — those who'd seen vaping on social
media are **2.3× more likely** to be in this no-friends group (as if social media stood in for
friends). Among teens who **don't** vape even though vaping surrounds them (the "resisters"), Asian
kids resist strongly, while girls and kids who feel depressed resist **less**. One number about
parents' education looks meaningful, but it turns out to be an artifact of how "Don't know" answers
were coded — we clean that up next.*

# 5. Two corrected re-analyses (`_AO`)

The two issues below are reproduced from the current report and then **corrected**, with the fix
shipped as a labelled script. Each is explained intuitively first, then technically with the
relevant algebra.

## 5.1 Threshold — the grade-period block is spurious (complete separation)
*Correction: [`codes/07-threshold-corrected-AO.R`](../codes/07-threshold-corrected-AO.R) → `outputs_AO/model/threshold_corrected-AO.csv`*

**Intuitively.** A "threshold" answers: *of the friends you could be influenced by, how many were
already using when you started?* The model needs, for each adopter, a count of friends-already-using
$k_i$ out of friends-total $n_i$. But the network only starts at grade-period 3 (GP3). For someone
who adopts *at* GP3 there is no earlier period in the data, so **none** of their friends can have
adopted before them: $k_i = 0$ for **every** GP3 adopter (149 non-isolates; 44 survive into the
complete-case fit — all with $k=0$ either way). When an entire time category contains only zeros, the
model can fit it perfectly by sending that category's coefficient to $\pm\infty$. That is exactly
what happens: the report's threshold column shows a Constant of **−20.63** and grade-period rows of
**17.23, 17.88, 18.11, 18.46, 18.84** — these are not findings, they are the optimiser running off to
infinity (ORs of $10^{7}$–$10^{8}$). Yuchan already spotted the *same* pathology for School 105 and
dropped it; the grade-period block has it too, undiagnosed — **and a second copy survives in School
105 even after GP3 is removed.**

**Technically.** The model is the grouped binomial of §3.1,
$k_i \mid n_i \sim \mathrm{Binomial}(n_i,\pi_i)$,
$\operatorname{logit}\pi_i = \alpha_{\tau_i} + \xi_{s(i)} + \mathbf{x}_i^\top\boldsymbol\beta$,
with GP3 the reference level absorbed into the intercept. Because $\tau_i = 3 \Rightarrow k_i = 0$,
the reference stratum is *perfectly separated*: the likelihood increases monotonically as
$\alpha_{\text{GP3}}\to-\infty$ with the GP$\ge 4$ dummies diverging to compensate; the MLE does not
exist (printed SE $\approx 750$). The same holds for School 105 (9 rows, all $k=0$). Two fixes:

- **B. Exclude GP3 adopters** ($\tau_i>3$), as the retired susceptibility model already did. This is
  information-lossless for the slopes — algebraically, $B_{\text{intercept}}-A_{\text{intercept}}$
  equals A's GP4 coefficient and B's GP dummies are A's GP *differences* — so every substantive OR is
  bit-identical to the as-is fit (max difference $10^{-14}$). *But it leaves the School-105 separation
  in place.*
- **C. Firth-penalized GLM** (`brglm2`, guaranteed-finite estimates) on the **full** 406 adopters —
  **the recommended primary fix**: it keeps all the data, resolves *both* separations, and recovers
  the genuine steep grade-period trend. Its bias-reducing shrinkage moves a few coefficients at the
  2nd–3rd decimal, so the ORs are *substantively similar*, not bit-identical.

| parameter | A. as-is (n=406) | B. excl. GP3 (n=362) | C. Firth, full (n=406) |
|:--|--:|--:|--:|
| GP4 / GP5 / GP8 | $10^{7}$–$10^{8}$ *(sep.)* | — / 1.91\* / 5.02\*\* | 40.7\*\* / 76.4\*\* / 196.7\*\* |
| Perceived friend use | 1.20 \*\*\* | 1.20 \*\*\* | 1.20 \*\*\* |
| Asian | 0.598 (.05) | 0.598 (.05) | 0.604 (.05) |
| Cohort C2 | 2.64 * | 2.64 * | 2.43 * |
| Hispanic / GAD / MDD / sex_min | ns | ns | ns |
| separation remaining | Constant + all GP + Sch105 | **Sch105** | **none** |

The substantive ORs agree across all three (PFU $\approx 1.20$, Asian $\approx 0.60$, Cohort
$\approx 2.5$) — the separation never threatened those conclusions — but only C produces a fully
interpretable, separation-free fit. The corrected story is itself substantive: **the adoption
threshold rises steeply across grade-periods** (later adopters require far more friends already
using), which the degenerate as-is block could not express. Tables 2–3 should report **model C**; the
`−20.63 / 17.x` rows must be removed.

## 5.2 Social media — fragile and spec-dependent, on an under-powered design
*Correction: [`codes/08-socialmedia-power-AO.R`](../codes/08-socialmedia-power-AO.R) → `outputs_AO/model/socialmedia_corrected-AO.csv`*

**Intuitively.** The question is whether *social-media* e-cig exposure explains adoption among
peer-independent students. The catch: the social-media items exist only in **W4–W8**, so the 224
students who adopt at GP3 (= W3) — the single largest cohort — have **no** measurement at adoption and
drop out (only **456 of 841 adopters, 54%**, are measurable). Table 1a then shows two significant
cells — No Exposure **OR 2.31** and No Perceived Friends **OR 0.40** — pointing in *opposite*
directions, which already warns that this is not one clean "alternative channel". But the right
reading is **not** "artifact": the No-Exposure association is **fragile and specification-dependent**.
It survives a Firth penalty (OR 2.12, $p=.021$, so it is not an over-fitting blow-up), yet it
attenuates to **borderline** (OR ≈ 1.5, $p≈.10$) once we adjust for *temporally-valid baseline
confounders* instead of the model's `friends_ecig_toa` term — which is measured **at the moment of
adoption**, i.e. post-treatment relative to the outcome, so controlling for it risks collider bias.
Whether the No-Exposure effect is "real" depends on a confounder-vs-collider judgement the data
cannot settle. The No-Perceived-Friends cell, by contrast, is robust *and* protective *and*
opposite-signed — which is the opposite of what the channel hypothesis predicts.

**Technically.** With $\mathcal A$ the adopters and
$\mathcal A_{\mathrm{SM}}=\{i:\mathrm{SM}_i\text{ observed at }\tau_i\}$, collection only in
$t\in\{4,\dots,8\}$ gives $\tau_i=3\Rightarrow\mathrm{SM}_i$ missing and
$|\mathcal A_{\mathrm{SM}}|/|\mathcal A| = 456/841 = 0.54$, concentrated in the earliest cohort
(non-ignorable). The "collapse" of the No-Exposure cell decomposes into (a) covariate adjustment and
(b) a $239\!\to\!455$ listwise-deletion sample change — on a *fixed* sample the parsimonious OR is
1.52, not 1.26 — so it is not a single clean effect. The load-bearing covariate is `friends_ecig_toa`
(post-treatment); using only baseline GP3 confounders gives the fair middle estimate. As a power
floor, the post-hoc Wald minimum detectable OR,
$\mathrm{OR}_{\min}=\exp[(z_{.975}+z_{.80})\,\widehat{\operatorname{SE}}]$, is **≈ 2.0–2.2** for the
null cells: those nulls are *uninformative*, not evidence of absence. Specification ladder (each cell
OR ($p$); EPV = events-per-variable):

| subgroup | 1. full (Yuchan) | 2. Firth | 3. baseline GP3 *(fair)* | 4. parsimonious | $p_{\text{Bonf}}$ |
|:--|--:|--:|--:|--:|--:|
| **No Exposure** | 2.31 (.013), EPV 7 | 2.12 (.021) | **1.51 (.10)**, EPV 19 | 1.26 (.25) | .050 |
| **No Perceived Friends** | 0.40 (.010) | 0.43 (.017) | **0.55 (.017)** | 0.51 (.001) | .041 |
| Neither | 0.87 (.70) | 0.87 (.68) | 0.89 (.66) | 0.67 (.08) | 1.00 |
| Low Threshold | 1.39 (.31) | 1.35 (.35) | 1.11 (.67) | 0.95 (.82) | 1.00 |

(A cumulative "ever exposed W4..TOA" measure recovers 592/841 adopters but is mechanically tied to
adoption timing — more waves of opportunity for later adopters — so it is a *different* exposure
definition, not a cleaner version of the snapshot; it is reported in the CSV for completeness.)
**Conclusion.** The only social-media association stable across all specifications is the
**protective, opposite-signed** No-Perceived-Friends cell, which does *not* support
social-media-as-alternative-channel; the headline No-Exposure cell is **borderline and
specification-dependent**; and the genuinely peer-independent cells stay null with a detectable-OR
floor near 2. The honest statement is "fragile and under-powered — not yet answerable", not "social
media has (or lacks) an effect". A clean test needs a pre-registered baseline-confounder block, a
pooled W4–W8 exposure measure, and the GP3 adopters acknowledged as structurally unmeasurable.

*(Aside: the Hispanic threshold coefficient is null under school FE because individual ethnicity is
collinear with the Hispanic-majority school dummies — the signal sits in Schools 104/114. Worth
keeping in mind when reading Tables 2–3; not pursued in detail here.)*

## 5.3 The report's tables, with every fix applied
*Correction: [`codes/09-fully-corrected-tables-AO.R`](../codes/09-fully-corrected-tables-AO.R) → `outputs_AO/model/fully_corrected_*-AO.csv`*

This section re-runs the report's four tables with **all** the fixes at once: the `par_edu` recode
below, plus **Firth** for every model that suffered separation or small-sample fragility (threshold,
adoption, and the subgroup/social-media logistics).

**The `par_edu` fix.** Per the codebook the parent-education column is the **legacy 6-level scale at
every wave** (1 = ≤8th grade … 6 = advanced degree, **9 = "Don't know"**). Yuchan's pipeline (a) left
"Don't know" in as a numeric **9** — ranked *above* "advanced degree" — for the 7% of adopters who
gave it, and (b) applied a "new-scale" remap at W7+ that does not match this release. The corrected
recode homogenises the scale to a clean $\{1,\dots,6\}$ and sends "Don't know" to missing:
$$
\text{par\_edu}^{\text{corr}}_i \;=\; \text{first valid }v\in\{1,\dots,6\}\text{ across the student's waves},
\qquad 9 \mapsto \text{NA}.
$$

**Corrected Tables 2–3 (Adoption — Firth FE; Threshold — Firth grouped binomial):**

| Predictor | Adoption (Firth FE) | Threshold (Firth) |
|:--|--:|--:|
| Network exposure | **2.34** \*\* | — |
| Perceived friend use | **1.44** \*\* | **1.20** \*\* |
| MDD | **1.61** \*\* | 0.78 |
| GAD | **0.79** \* | 1.02 |
| Cohort C2 | 1.05 | 2.34 (.06) |
| Asian | 0.80 | 0.61 (.06) |
| **par_edu (recoded)** | 0.98 *(ns)* | 1.04 *(ns)* |
| Constant; GP4 / GP5 / GP8 | GP8 = **0.007\*\*** *(finite)* | **40\*\* / 70\*\* / 189\*\*** *(finite)* |

The degenerate blocks are gone: the threshold grade-period effects are now **finite, large and
increasing** (40 → 189), i.e. the genuine result the separation hid — *later adopters start only once
many more friends already use*. Every substantive OR is unchanged (exposure 2.34, PFU 1.20–1.44, MDD
1.61, GAD 0.79, Asian ≈ 0.6).

**Corrected Tables 1a–1b (subgroups — Firth + recoded `par_edu`):** the conclusions of §4.2 hold —
Asian elevated among peer-independent users (Neither **3.45\***, Low Threshold **3.25\***) and strongly
protective among resisters (**0.51 / 0.34 / 0.36**, all \*\*); MDD/female raise resistance. Two changes
matter:

| | as-is (§4.2) | corrected (§5.3) |
|:--|--:|--:|
| `par_edu` in "High Perceived Friends" | **0.93 \*\*** (looked real) | **1.01 (ns)** — spurious effect removed |
| `sm_ecig_any` in "No Exposure" | 2.31 \* | 2.40 \*\* (survives Firth) |

So the `par_edu` "effect" was manufactured by the "Don't know" coding bug and **vanishes** once it is
recoded. The social-media "No Exposure" cell, by contrast, **survives** Firth — confirming (as §5.2
already showed) that it is *not* a small-sample artifact; its fragility lives in the confounder choice
(`friends_ecig_toa` is measured at adoption), which Firth does not address. Net: Firth fully fixes the
threshold; for social media it removes one doubt (small-sample) but leaves the design under-powered
and the No-Exposure estimate confounder-dependent.

**For grandma.** *We re-ran everything with the dials set correctly. First we fixed how "parents'
education" was recorded: people who answered "I don't know" had been counted as if that were the
**highest** education level — so we moved those to "missing". The moment we did, a result that looked
real (something about parents' education among the resisters) **disappeared** — it was never real,
just a coding slip. Second, we used a sturdier statistical method (Firth) that stops the computer from
producing those impossible giant numbers. With it, the real pattern shows cleanly: kids who start
vaping later in high school only do so once **many** of their friends already vape. None of the
genuine findings changed — friends who vape, feeling depressed, and being surrounded by vaping still
drive adoption. The one stubborn loose end is social media: the sturdier method confirms the
"no-friends vapers saw more vaping online" result is not a fluke of small numbers, but we still can't
be sure it's real, because half the kids can't be measured for it and it leans on a questionable
control — so we mark it "promising but unproven", not "true".*

# 6. Reproducibility

Reproduced end-to-end on the 042326 release with **netdiffuseR 1.25.0** (CRAN). The in-development
1.26.0 aborts at `hazard_rate` because `c.diffnet` drops the new `$status` field; root cause and a
minimal reproducer are filed for the maintainer (git-ignored `playground/`).

```sh
# 1) Yuchan's pipeline (produces her Tables 1–3 and caches the analytic frames)
Rscript codes/260612_netdiffuseR_past6mo_revised13.R     # data_path → the two ADVANCE_W*.xlsx
# 2) The corrections (consume the analytic frames; write to outputs_AO/model/)
Rscript codes/07-threshold-corrected-AO.R       # threshold separation (§5.1)
Rscript codes/08-socialmedia-power-AO.R         # social-media power/robustness (§5.2)
Rscript codes/09-fully-corrected-tables-AO.R    # all fixes at once: par_edu + Firth (§5.3)
```

The `_AO` scripts read the analytic frames built by Yuchan's pipeline (set `AO_THR_SUB` / `AO_CACHE`
to point at them); script 09 additionally recomputes `par_edu` from the raw XLSX. None alter her
variable construction otherwise.

# 7. Open items

1. **Adopt the corrected threshold** (§5.1/§5.3): report the Firth fit — the `−20.63 / 17.x`
   grade-period block is a separation artifact, not a result.
2. **Re-frame the social-media result** (§5.2): the positive No-Exposure cell survives Firth but is
   confounder-dependent and the design is under-powered. Pool W4–W8 exposure, use baseline-only
   confounders, and pre-register.
3. ~~Recode `par_edu`~~ — **done (§5.3)**: "Don't know" (9) → NA, homogeneous 1–6 scale. Roll the
   recode into the upstream pipeline so all tables use it by default.
4. **Lag the threshold covariates** to $\tau_i - 1$ to match the disadoption convention and blunt the
   simultaneity concern (§2.5).
5. Decide the **Hispanic estimand** (school-FE contextual vs no-FE individual) and state it once.

---

*Companion document: `A - Network-Disadoption/docs/disadoption-study.md`. Source report:
`docs/260612_threshold_report_final.pdf`. Corrected analyses: `codes/0{7,8,9}-*-AO.R`,
`outputs_AO/model/*-AO.csv`.*
