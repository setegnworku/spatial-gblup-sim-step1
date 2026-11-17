# Steps 2 and 3: objectives, design and the questions they answer

<sub>[README](README.md) &nbsp;·&nbsp; [Reproduce](REPRODUCE.md) &nbsp;·&nbsp; [The population](docs/POPULATION.md) &nbsp;·&nbsp; **Plan for steps 2–3**</sub>

> [!IMPORTANT]
> **Plan only. This page contains no results of steps 2 and 3.** All data are SIMULATED. Step 1 and its results are
> in the [README](README.md).

## 1. Where steps 2 and 3 fit

The study is a **ladder: one change per step, on the same simulated cows**. Every step uses the same founders, bulls,
cows, true breeding values, herds, herd locations and spatial field; only the records change.

| Step | Records per cow | Permanent environment (pe) | What it adds |
|---|---|---|---|
| **1** (done, see README) | one test-day record | none | the baseline: what the models do with one record per cow |
| **2** | about 10 test-day records (the real recording pattern) | none | the effect of **repeated records** alone |
| **3** | the same records as step 2 | **added** (σ²pe = 0.04) | the effect of a **permanent environment**, the situation of the real data |

*Columns.* **Records per cow:** test-day records simulated for each cow. In steps 2 and 3, each cow copies the test
dates, calving season, age and stage of lactation of a real cow (10.2 records on average). **Permanent environment:**
an effect of the cow that is the same on all her records but is not passed to her daughters (how she was reared and
managed). **What it adds:** the single change from the step before.

Because each step changes one thing, a difference between step 1 and step 2 is due to the number of records, and a
difference between step 2 and step 3 is due to the permanent environment.

**True values** (the same in every step except pe):

| Parameter | Value |
|---|---|
| Genetic variance σ²g | 0.06 |
| Permanent environmental variance σ²pe | 0 in steps 1–2; 0.04 in step 3 |
| Herd variance σ²h | 0.15 |
| Spatial field variance σ²s (range) | 0.41 (33 km) |
| Residual variance σ²e | 0.29 |
| Heritability h² | 0.066 in steps 1–2; 0.063 in step 3 |
| Repeatability (σ²g + σ²pe) / total | 0.066 in steps 1–2; 0.105 in step 3 |

*Columns.* **Parameter:** the simulated true value, on the scale of standardised milk yield. The genetic and pe values of
step 3 split the genomic variance of the published GHS model (about 0.10–0.11, Houaga et al. 2026) into a genetic part
and a permanent-environment part. That split is an assumption of the simulation; the real split is unknown.

## 2. The questions

### Step 2: what do repeated records change?

| | Question | Why it matters |
|---|---|---|
| **Q2.1** | With about 10 records per cow, is the genetic variance estimated without bias, and does the prior still matter? | In step 1, INLA's default prior biased the genetic variance with one record per cow. |
| **Q2.2** | When a model leaves out the herd or the field, where does that variance go now? | In step 1 it went to the residual. With repeated records, a cow's herd and location are the same on all her records, so the cow effect could take it. |
| **Q2.3** | Does fitting a pe term cost anything when there is no pe? | It checks that a model with pe (GPHS) is safe to use when pe is absent. |
| **Q2.4** | Do repeated records on relatives improve the breeding values of cows without records? | Young cows and heifers are evaluated from their relatives' records. |

*Columns.* **Question:** what step 2 answers. **Why it matters:** its link to step 1 or to practice.

### Step 3: what does a permanent environment change?

| | Question | Why it matters |
|---|---|---|
| **Q3.1** | If pe exists but the model has no pe term (as in the published GHS), where does its variance go: to the genetic effect or to the herd? | If it goes to the genetic effect, the reported heritability is really a repeatability. |
| **Q3.2** | Does the model with pe, herd and field (GPHS) recover every variance component? | It tells whether the full model can be recommended. |
| **Q3.3** | Without a pe term, are the breeding values over-dispersed, for recorded cows and for cows without records? | Over-dispersed breeding values exaggerate differences between animals and the expected response to selection. |
| **Q3.4** | How much does the spatial field add to the accuracy of breeding values when pe is present, compared with the gain in predicting phenotypes? | The real-data validation measured the prediction of phenotypes, which includes the environment. |
| **Q3.5** | Does the prior matter when there are repeated records? | It tells whether the prior advice from step 1 also applies to test-day data. |
| **Q3.6** | How many cows are needed (1,000, 2,000 or 4,000) to separate the genetic from the pe variance? | A design rule for recording schemes. |

*Columns.* **Question:** what step 3 answers. **Why it matters:** its link to the real data or to practice.

### Two follow-up designs on step 3 (2,000 cows)

| | Question | Design |
|---|---|---|
| **F-a** | How many relatives are needed to separate genetic and pe variance? | 5, 2 or 1 daughters per bull instead of about 10; the same cows, herds, field and records, only the sires change. |
| **F-b** | Sensitivity analysis: does the herd term still matter when nearly every herd has one cow, and which term then takes the pe? | 30, 40, 90 and 99 % single-cow herds, beyond the main design (50–80 %). |

*Columns.* **Question:** what the follow-up answers. **Design:** what is changed compared with step 3.

## 3. Design

| Part | Records; pe | Cows | Single-cow herds | Replicates | Models | Prior |
|---|---|---|---|---|---|---|
| Step 2 | ~10; no | 1,000; 2,000 | 50, 60, 70, 80 % | 10 | G, GH, GS, GHS, GPHS | INLA default |
| Step 3 | ~10; yes | 1,000; 2,000 (4,000 at 70 %) | 50, 60, 70, 80 % | 10 (4,000 cows: 4) | all eight | INLA default |
| Step 3, prior comparison | ~10; yes | 1,000; 2,000 | 50, 60, 70, 80 % | 4 | all eight | PC prior (same data refitted) |
| F-a: daughters per bull | ~10; yes | 2,000 | 70 % | 5 | GHS, GPHS | INLA default and PC |
| F-b: herd structure | ~10; yes | 2,000 | 30, 40, 90, 99 % | 5 | GS, GHS, GPS, GPHS | INLA default |
| Cross-validation | as in each step | 1,000; 2,000 | 50, 60, 70, 80 % | 5 | G, GH, GS, GHS; in step 3 also GPH and GPHS | INLA default |

*Columns.* **Part:** the analysis. **Records; pe:** records per cow and whether a true pe effect exists. **Cows:** cows
per data set. **Single-cow herds:** the share of herds with one cow; the real data have 72 %. **Replicates:**
independently simulated populations per combination of cows and share. **Models:** see Section 4. **Prior:** the prior
on the variance components (Section 4 of the README).

**Cross-validation** is done exactly as for step 1 (README §4). Each cow is a test cow once. Each measure is computed
within each test group, averaged over the 5 groups of a data set, then averaged over the replicates of each scenario.
**Every share of single-cow herds is reported separately.** Besides the accuracy r(EBV, a), the **dispersion** of the
test cows' breeding values (the slope of the true on the estimated breeding value) is computed in the same way, and
every cow's estimated and true breeding value is saved.

## 4. The eight models

| Model | Genomic animal effect | Permanent environment | Herd | Spatial field |
|---|---|---|---|---|
| G | yes | – | – | – |
| GH | yes | – | yes | – |
| GS | yes | – | – | yes |
| GHS | yes | – | yes | yes |
| GP | yes | yes | – | – |
| GPH | yes | yes | yes | – |
| GPS | yes | yes | – | yes |
| GPHS | yes | yes | yes | yes |

*Columns.* **Model:** G, GH, GS and GHS are the models of Houaga et al. (2026); P adds a pe effect. **yes / –:** whether
the model contains the term. In step 2 the true model is GHS; in step 3 it is GPHS. Every other model deliberately
leaves out at least one true term, so that the comparison shows where the variance of a missing term goes.

## 5. What we expected, written before the results

These expectations were written before any result of the part concerned (dates in brackets) and are kept unchanged,
so that the results can be judged against them.

**Step 2 and the cross-validation** (7 November 2026):

| | Expectation |
|---|---|
| **L5** | GHS and GPHS are both unbiased for σ²g, and GPHS estimates σ²pe close to 0. The spread of σ̂²g over replicates is about half that of step 1. |
| **L6** | From step 2 to step 3, GHS's σ̂²g rises from about 0.06 to about 0.10: the whole bias seen in step 3 is the pe, not the repeated records. |
| **L9** | Steps 2 and 3 give cows without records only slightly better breeding values than step 1 (+0.02 to +0.05). |

**Step 3** (written for the step-3 design, before its results):

| | Expectation |
|---|---|
| **E1** | Without a pe term, σ̂²g ≈ σ²g + σ²pe: for GHS about 0.09–0.11 (+50 to +80 %), h² about 0.10–0.12, in every scenario. More cows do not remove this bias. |
| **E2** | GPHS recovers σ²g within ±15 % in every scenario. |
| **E3** | The share of single-cow herds decides how well pe and herd separate: from 50 to 80 %, the spread of σ̂²pe and σ̂²h at least doubles and they trade off; σ̂²g is hardly affected. |
| **E4** | Missing environmental terms move into the cow terms: G's σ̂²g ≥ 0.2; GS's σ̂²g grows with the single-cow share; GH's herd variance takes the field (about 0.45). |
| **E5** | From 1,000 to 4,000 cows, the spread of GPHS's σ̂²g over replicates roughly halves. |
| **E6** | Accuracy: GPHS beats GHS by at most +0.02. Dispersion: GHS over-dispersed (slope about 0.7–0.8), GPHS about 1. |
| **E7** | WAIC prefers GPHS over GHS, more clearly with more cows. |

**F-a, daughters per bull** (13 November 2026): at 10 and 5 daughters per bull the separation works with both priors;
at 2 and 1, σ̂²g moves into σ̂²pe with the default prior; with PC priors the estimates stay closer to the truth but
with wide intervals; σ̂²g + σ̂²pe stays right at every level.

**F-b, herd structure** (13 November 2026):

| | Expectation |
|---|---|
| **F1-E1** | With a pe term, omitting the herd matters less and less: at 99 % single-cow herds, GPS and GPHS fit almost equally well. |
| **F1-E2** | Without a pe term, GS overestimates σ²g (about 0.40) at every share. |
| **F1-E3** | At 99 %, GHS can no longer separate herd from cow, so part of the herd variance enters its genetic term and GHS becomes similar to GS. |
| **F1-E4** | GPHS's σ̂²g stays within ±15 % at every share, while its herd and pe variances trade off near 100 %. |
| **F1-E5** | GPHS accuracy falls slowly with the single-cow share (about 0.63 at 50 % to about 0.58 at 99 %). |

**Cross-validation with dispersion** (16 November 2026): in step 3 the test cows' breeding values from GHS are
over-dispersed (slope about 0.6–0.8) and those from GPHS are close to 1; in step 2 the slopes are close to 1; with one
record per cow and 1,000 cows the slope of GHS is above 1 (the default prior shrinks too much); G and GS give slopes far
below 1 in steps 2 and 3.

*Columns (all tables in this section).* First column: the label of the expectation. **Expectation:** what was predicted,
in words, before the results.

## 6. How the steps answer the study's hypotheses

| Hypothesis | Answered by |
|---|---|
| The spatial field improves the accuracy of genomic breeding values, but by less than it improves the prediction of phenotypes. | Steps 1–3, cross-validation (Q2.4, Q3.4) |
| When pe exists but is not fitted, its variance goes to the genetic effect, so the heritability is a repeatability and the breeding values are over-dispersed. | Steps 2 → 3 (Q3.1–Q3.3), F-b |
| INLA's default prior biases the genetic variance when information is weak; a PC prior scaled to the trait avoids this. | Step 1, step 3 prior comparison (Q3.5), F-a |
| The number of recorded cows matters more than the share of single-cow herds. | Steps 1–3 (Q3.6), F-a, F-b |

*Columns.* **Hypothesis:** the study's hypothesis, in words. **Answered by:** the steps and questions above that test it.

## 7. Status

- Steps 2 and 3, the prior comparison and the two follow-up designs have been simulated and fitted with the same
  pipeline as step 1.
- The cross-validation with the dispersion slope is running (started 16 November 2026).
- The results are being reviewed and will be shared in a later update.
