<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/banner-dark.svg">
    <img alt="Spatial genomic models for smallholder dairy cattle. Step 1 of 3: one record per cow. Simulated data." src="docs/assets/banner-light.svg" width="100%">
  </picture>
</p>

<p align="center">
  <img alt="Data: simulated" src="https://img.shields.io/badge/data-SIMULATED-eb6834?style=flat-square">
  <img alt="Step 1: complete" src="https://img.shields.io/badge/step%201-complete-1baf7a?style=flat-square">
  <img alt="Steps 2 and 3: plan only" src="https://img.shields.io/badge/steps%202%E2%80%933-plan%20only-8b8b8b?style=flat-square">
  <img alt="R 4.5.2" src="https://img.shields.io/badge/R-4.5.2-2a78d6?style=flat-square&logo=r&logoColor=white">
  <img alt="R-INLA 26.08.07" src="https://img.shields.io/badge/R--INLA-26.08.07-2a78d6?style=flat-square">
  <img alt="inlabru 2.15.0" src="https://img.shields.io/badge/inlabru-2.15.0-2a78d6?style=flat-square">
  <img alt="AlphaSimR 2.1.0.9014" src="https://img.shields.io/badge/AlphaSimR-2.1.0.9014-2a78d6?style=flat-square">
  <a href="REPRODUCE.md"><img alt="Reproducible from the founder genomes" src="https://img.shields.io/badge/reproducible-from%20the%20genome-52514e?style=flat-square"></a>
</p>

<p align="center">
  <a href="#at-a-glance"><b>At a glance</b></a> &nbsp;·&nbsp;
  <a href="#1-why-we-simulate">Why simulate</a> &nbsp;·&nbsp;
  <a href="#6-results">Results</a> &nbsp;·&nbsp;
  <a href="#7-what-step-1-shows">What step 1 shows</a> &nbsp;·&nbsp;
  <a href="REPRODUCE.md">Reproduce</a> &nbsp;·&nbsp;
  <a href="docs/POPULATION.md">The population</a> &nbsp;·&nbsp;
  <a href="PLAN_STEPS_2_3.md">Plan for steps 2–3</a>
</p>

> [!IMPORTANT]
> **All results in this repository come from SIMULATED data.** No real milk record, animal or farm location is
> included. The simulator borrows the *structure* of a real data set (herd locations, herd sizes, recording pattern);
> that input file is not shared (see [Data not included](#10-data-not-included)).

> [!NOTE]
> **Status (16 November 2026):** step 1 is final (85 data sets, pooled 9 November 2026). The refit of the same data
> with PC priors is final (80 data sets, 10 November 2026). The cross-validation results use the agreed method:
> accuracy computed within each test group, then averaged (§4, §6.4).
> **Steps 2 and 3:** their objectives, design and the questions they answer are in
> [PLAN_STEPS_2_3.md](PLAN_STEPS_2_3.md). That page is a plan only; it holds no results.

> [!TIP]
> **To replicate step 1 from the genome onwards, follow [REPRODUCE.md](REPRODUCE.md).** It covers every step (founder
> genomes, population, genotypes and G, herds, spatial field, records, fits, scoring, cross-validation, pooling,
> figures), with seeds and checkpoints. The founder genomes and a synthetic stand-in for the real-data backbone are
> included, so the whole pipeline runs with `bash run_step1.sh`.

## At a glance

<p align="center">
  <img src="figures/fig0_records_vs_breeding_values.png" width="100%"
       alt="Dot plot for 2,000 cows in cross-validation. Prediction of records rises from 0.30 with G to 0.67 with GHS; accuracy of breeding values rises only from 0.27 to 0.34.">
</p>

Step 1 gives each simulated cow one record and fits the four models of Houaga et al. (2026). Four findings:

<table>
<tr>
<td width="50%" valign="top">

**1 · Weak data need a sensible prior**<br>
With INLA's default prior, 15 of 40 GHS fits at 1,000 cows collapsed to σ̂²g ≈ 0. A PC prior on the same data
removed all 106 collapses and the bias ([§6.6](#66-refit-with-pc-priors-the-same-80-data-sets)).

</td>
<td width="50%" valign="top">

**2 · The spatial field earns its place**<br>
GHS had the best WAIC in 82 of 85 data sets, got the herd variance right (0.15, against 0.50 in GH) and added 0.04 to
0.05 to the accuracy of breeding values ([§6.3](#63-accuracy-and-model-fit)).

</td>
</tr>
<tr>
<td width="50%" valign="top">

**3 · Judge breeding values, not records**<br>
GHS predicts hidden records with cor(y, ŷ) 0.64–0.68, but the breeding values of the same cows have accuracy only
0.29–0.34: most of the record prediction is environment ([§6.4](#64-cross-validation-cows-predicted-without-their-own-records)).

</td>
<td width="50%" valign="top">

**4 · Cows matter more than herd structure**<br>
With the default prior σ̂²g was right only with 4,000 cows. The share of single-cow herds did not move it, but the
field helped most at 80 % single-cow herds ([§6.5](#65-by-share-of-single-cow-herds)).

</td>
</tr>
</table>

**Contents**
1. [Why we simulate](#1-why-we-simulate)
2. [The ladder, and why step 1 comes first](#2-the-ladder-and-why-step-1-comes-first)
3. [What is simulated](#3-what-is-simulated)
4. [Models, priors and validation](#4-models-priors-and-validation)
5. [Expected results, written before the results](#5-expected-results-written-before-the-results)
6. [Results](#6-results)
7. [What step 1 shows](#7-what-step-1-shows)
8. [What step 1 cannot tell](#8-what-step-1-cannot-tell)
9. [Repository layout and how to run](#9-repository-layout-and-how-to-run)
10. [Data not included](#10-data-not-included)
11. [Software and references](#11-software-and-references)

---

## 1. Why we simulate

Houaga et al. (2026, *Genetics Selection Evolution*) analysed test-day milk records of dairy cows on smallholder farms
in Tanzania with genomic models that add a **herd effect** and a **spatial field** (an SPDE Gaussian field in
R-INLA/inlabru). Their four models are G, GH, GS and GHS (G = genomic animal effect, H = herd, S = spatial field). The
spatial models fitted the data better and predicted hidden records better.

In real data, however, the **truth is unknown**. Nobody can check whether the estimated genetic variance is right,
whether herd and location effects are separated correctly, or whether better prediction of records means better
breeding values. A simulation can, because:

- we **choose the true values** (genetic, herd, spatial and residual variances; every cow's true breeding value);
- we build data that **look like the real data**: herd locations, herd sizes and recording pattern are resampled from
  the AADGG Tanzania data set of Houaga et al. (2026) (details in [Data not included](#10-data-not-included));
- we fit **the same models with the same software and settings** (Houaga's fixed effects, mesh recipe and field prior,
  and INLA's default prior on the other variances, as in his script);
- we compare every estimate with the truth, over many replicates.

Smallholder data have features that make this check necessary:
- **most herds have one cow** (about 70 % of herds in the real data), so a herd effect and a cow's own effect are hard to
  tell apart;
- herds are spread over a large area, and **neighbouring farms share** climate, feed, markets and services: this is
  what the spatial field is meant to capture;
- cows have **few close relatives** in the data, and bulls are used across the whole area;
- the variance components are estimated in a Bayesian framework, so **the prior** can matter when the data are weak.

## 2. The ladder, and why step 1 comes first

The real situation mixes several things at once: repeated records per cow, a cow **permanent environment** (pe),
crossbreds. To see what each of them does, the study is built as a **ladder: one change per step, on the same
simulated cows** (same seeds, same founders, same bulls, same breeding values, same herds and the same spatial field).

| Step | Records per cow | True permanent environment | Models fitted | In this repository |
|---|---|---|---|---|
| **1** | **one** (one random test-day) | none | G, GH, GS, GHS | **yes** |
| 2 | about ten | none | G, GH, GS, GHS + GPHS | no |
| 3 | about ten | σ²pe = 0.04 | all eight (with and without pe) | no |

*Columns.* **Step:** the rung of the ladder. **Records per cow:** test-day milk records each simulated cow has.
**True permanent environment:** the variance of the cow's lifelong non-genetic effect put into the simulation (σ²pe;
it repeats on all her records but is not shared with relatives). **Models fitted:** the models compared in that step
(letters: G genetic, H herd, S spatial field, P permanent environment). **In this repository:** whether its code and
results are here.

**Step 1 is the baseline.** With one record per cow nothing repeats, so there is no permanent environment and only
four kinds of variation are in play. Each is defined by **who shares it**:

| Term | Shared by | True variance |
|---|---|---|
| genetic (animal) | relatives, through the genomic relationship matrix G | 0.06 |
| herd | herd-mates | 0.15 |
| spatial field | neighbouring farms (smoothly, with distance) | 0.41 |
| residual | nobody | 0.29 |

*Columns.* **Term:** a source of variation in milk yield. **Shared by:** which cows have similar values of it (this is
what lets a model tell the terms apart). **True variance:** its variance in the simulation, on the scale of
standardised milk yield (total about 0.91).

If the models cannot separate these four with one record per cow, adding repeated records and a pe will not make the
picture clearer. And once step 1 is understood, each later step can be read as "what changes when one thing is added".

## 3. What is simulated

Script: [`scripts/questions/sb_simulate.R`](scripts/questions/sb_simulate.R).

| Piece | Setting |
|---|---|
| Breed | one population: MaCS cattle history, 1,000 founders, 2 generations of random mating |
| Genome | 29 chromosomes; 5,800 QTL (not on the chip); SNP chip of about 40,000 markers (MAF ≥ 1 %) |
| Bulls | one bull per 10 cows (100 bulls for 1,000 cows), used unevenly (log-normal popularity), **not genotyped**; on average 10 daughters per bull |
| Cows | N daughters, one per random dam, **all genotyped**; G* = 0.99 G + 0.01 I (VanRaden G), as in Houaga et al. |
| Herds | a chosen share of herds has **one cow**; the other herds take sizes from the real multi-cow herds (mostly 2–3 cows, mean 2.3) |
| Herd locations | resampled real herd locations, each moved by a random 1 km jitter |
| Records | each cow takes the recording pattern of a random real cow (calendar month, calving season, lactation group, age, lactation-curve covariates); **step 1 keeps one random test-day of it** |
| Fixed effects | true sizes from the real data, without the breed-class part (one breed) |
| Spatial field | Matérn (ν = 1) Gaussian field, range 33 km, drawn on a mesh of the real herd locations; the same field for every cell of a replicate |
| Bulls and environment | cows are assigned to herds **at random**: no link between genetics and environment |

*Columns.* **Piece:** one part of the simulated world. **Setting:** how it is simulated in step 1 (QTL = the DNA
positions that truly affect milk yield; SNP chip = the markers used to build the genomic relationship matrix G; MAF =
minor allele frequency; jitter = a random shift that hides the exact real location).

**True values (step 1):**

| σ²g | σ²pe | σ²h | σ²s (SPDE parameter) | Range | σ²e | h² |
|---|---|---|---|---|---|---|
| **0.06** | 0 | 0.15 | 0.41 (about 0.35 across the herds) | 33 km | 0.29 | **0.066** |

*Columns.* **σ²g:** genetic variance, how much cows differ in the genes they pass on. **σ²pe:** permanent-environment
variance (none in step 1). **σ²h:** herd variance, what cows of one farm share (management, feeding). **σ²s:** variance
of the spatial field, what neighbouring farms share (climate, feed, markets); 0.41 is the SPDE parameter, and across
the actual herd locations the field varies by about 0.35. **Range:** distance at which two farms' environments become
almost unrelated (correlation about 0.1). **σ²e:** residual variance, what is left in each record. **h²:**
heritability, σ²g ÷ the sum of all variances.

The values are calibrated to the variance components estimated on the real data. Heritability is low, as in the
real data.

**Design:** number of cows × share of single-cow herds. Replicate *r* uses the same founders and the same field in every
cell, so cells differ only in the design factor.

| Cell | Cows | Single-cow herds (of herds) | Herds (realised mean) | Cows in single-cow herds | Replicates |
|---|---|---|---|---|---|
| n1000_p50 | 1,000 | 50 % | 600 | 29.8 % | 10 |
| n1000_p60 | 1,000 | 60 % | 649 | 38.6 % | 10 |
| n1000_p70 | 1,000 | 70 % | 713 | 49.7 % | 10 |
| n1000_p80 | 1,000 | 80 % | 788 | 62.9 % | 10 |
| n2000_p50 | 2,000 | 50 % | 1,210 | 30.2 % | 10 |
| n2000_p60 | 2,000 | 60 % | 1,312 | 39.4 % | 10 |
| n2000_p70 | 2,000 | 70 % | 1,436 | 50.2 % | 10 |
| n2000_p80 | 2,000 | 80 % | 1,581 | 63.2 % | 10 |
| n4000_p70 | 4,000 | 70 % | 2,874 | 50.3 % | 5 |
| **Total** | | | | | **85 data sets** |

*Columns.* **Cell:** the design's short name, n = number of cows, p = share of single-cow herds (n1000_p70 = 1,000
cows, 70 % single-cow herds). **Cows:** cows in each data set. **Single-cow herds (of herds):** the share of herds that
have exactly one cow (the design factor). **Herds (realised mean):** the number of herds that resulted, mean over
replicates. **Cows in single-cow herds:** the share of cows that are alone in their herd (lower than the share of
herds, because the other herds hold several cows). **Replicates:** independent simulated populations per cell; one
replicate of one cell = one **data set**.

2,000 cows at 70 % is the closest to the real data set (about 1,900 cows).

**The population in detail:** [docs/POPULATION.md](docs/POPULATION.md) covers the generations, the bulls and their
daughters, the dams, how related the cows are, the herds, the map and the spatial field, and how many half-sisters live
close to each cow.

## 4. Models, priors and validation

**Models** (script [`scripts/questions/q_fit.R`](scripts/questions/q_fit.R)). All have the same fixed effects as Houaga
et al., without the breed class: test year-month, calving year-season, age and the lactation curve within lactation
group (`tyrmn + cyrsn + ageZ:lacgr + leg1:lacgr + leg2:lacgr`).

| Model | Random terms | Role in step 1 |
|---|---|---|
| G | animal | leaves out herd and field |
| GH | animal + herd | leaves out the field |
| GS | animal + field | leaves out the herd |
| **GHS** | animal + herd + field | **the true model** |

*Columns.* **Model:** the name lists the random terms (G genetic, H herd, S spatial field). **Random terms:** the
effects whose variances the model estimates ("animal" = the genetic effect, linked between relatives through G); every
model also has the fixed effects and a residual. **Role in step 1:** what the model is for here.

G, GH and GS each leave out a true term on purpose: they show **where the variance of a missing term goes**.

**Priors.** A Bayesian model needs a prior for every variance. With informative data the prior hardly matters; with weak
data it can decide the answer. Step 1 was fitted twice, on identical data sets.

- **INLA's default prior** (`--priors isidore` in the code; the flag is named after the script of Houaga et al., which
  sets no prior and so uses INLA's default): loggamma(1, 5×10⁻⁵) on each precision, which is an inverse-gamma(1, 5×10⁻⁵)
  prior on the variance. It looks vague on the precision scale, but on the variance scale it puts **61 % of its
  probability below 0.0001** and 99.9 % below 0.06.
- **PC prior** (`--priors pc`). "PC" means **penalised complexity** (Simpson et al. 2017). Each random effect is seen
  as an extension of a simpler base model without it (variance 0); the prior penalises the distance from that base model
  at a constant rate, which gives an exponential prior on the standard deviation. It is set by one statement,
  P(σ > U) = α. Here: P(σ > 0.5) = 0.5 for the animal and herd effects and P(σ > 1) = 0.5 for the residual. It puts
  **1.4 % of its probability below a variance of 0.0001**.
- The spatial field has Houaga's PC-Matérn prior in both: P(σ_field > 0.5) = 0.5 and P(range < 50 km) = 0.8.

| Prior on a variance | Probability below 0.0001 | Probability below 0.06 (true σ²g) | Median |
|---|---|---|---|
| INLA's default, inverse-gamma(1, 0.00005) | 61 % | 99.9 % | 0.00007 |
| PC prior, P(σ > 0.5) = 0.5 | 1.4 % | 29 % | 0.25 |

*Columns.* **Prior on a variance:** the prior distribution, before seeing any data. **Probability below 0.0001:** how
much of the prior's weight sits at a practically zero variance (where a fit "collapses"). **Probability below 0.06:**
the prior's weight below the true genetic variance. **Median:** the variance value the prior considers equally likely
to be exceeded or not.

**Integration:** empirical Bayes (`int.strategy = "eb"`), to make the 660 full-data fits and 800 cross-validation fits
affordable. On one data set it agreed with INLA's default integration to three decimals.

**Cross-validation** (script [`scripts/questions/q_validate.R`](scripts/questions/q_validate.R)), replicates 1–5 of
the 1,000- and 2,000-cow cells (40 data sets): the cows are split at random into 5 test groups; each group's records
are hidden in turn and predicted from a fit to the other four, so every cow is a test cow once and is predicted
**from her relatives only**.
- **Test cows** (masked cows): the cows whose records are hidden in a round. **Training cows:** the cows of the other
  four groups, whose records the model is fitted to.
- **How a result is computed:**
  1. In each test group, the measure (for example r(EBV, a)) is computed for that group's test cows.
  2. The 5 values are averaged, giving one value per data set.
  3. The data sets are averaged within each scenario (number of cows × share of single-cow herds).
  4. Standard errors come from replicate means, because the scenarios of one replicate share a simulated population.
- *Changed on 15–16 November 2026:* earlier versions of this page pooled the predictions of the five test groups into one
  correlation and split the test cows into "seen" and "unseen" herds. Both were dropped. Pooling mixes EBVs from
  five different fits, which lowered the step-1 accuracy (for example GHS with 1,000 cows: 0.249 pooled against 0.294).
- The **dispersion** of the test cows' EBVs (slope of the true on the estimated breeding value, computed in the same
  way) is being added in a rerun started on 16 November 2026; it will appear here when complete.

**What we measure** (scripts [`q_extract.R`](scripts/questions/q_extract.R) and
[`sb_aggregate.R`](scripts/questions/sb_aggregate.R)):

| Measure | Meaning | Good value |
|---|---|---|
| σ̂²g (bias) | mean estimate over data sets; bias = estimate ÷ truth − 1 | within ±15 % of 0.060 |
| Coverage | share of data sets whose 95 % interval contains the truth | about 95 % |
| Collapsed | fits with σ̂²g < 0.001, i.e. at the bottom of the prior (about 0.0001 with INLA's default) | 0 |
| r(EBV, a) | correlation of EBVs with the **true** breeding values: the ranking | higher |
| Slope of a on EBV | the scale of the EBVs | 1 (above 1: EBVs shrunk too much, σ̂²g too low) |
| ΔWAIC | WAIC of the model minus WAIC of GHS (+ = worse) | 0 for GHS; > 10 is a clear difference |
| cor(y, ŷ) | Houaga's validation measure: hidden record against its prediction | higher, but mostly measures the environment |

*Columns.* **Measure:** the name used in the result tables. **Meaning:** what it is. **Good value:** what a perfect or
acceptable model gives. (EBV = estimated breeding value from the model; a = true breeding value, known only in a
simulation; WAIC = a score of how well a model predicts the data, lower is better.)

## 5. Expected results, written before the results

These expectations were written on **7 November 2026, before any step-1 result**, and kept unchanged (the code
comments repeat them).

| | Expectation |
|---|---|
| **L1** | GHS, the true model, estimates σ²g without bias (mean within ±15 % of 0.06), but imprecisely. |
| **L2** | Herd and residual trade off (only multi-cow herds tell them apart), more so at 80 % single-cow herds. |
| **L3** | A model without a term misplaces it: G and GS put herd and field into σ²g; GH puts the field into the herd. |
| **L4** | Accuracy of the recorded cows' EBVs, r(EBV, a), of GHS is about 0.35–0.45. |
| **L7** | The EBVs of masked cows are far less accurate than those of training cows (guess 0.15–0.30 against 0.35–0.6); accuracy rises with the number of cows. |
| **L8** | cor(y, ŷ) rises by +0.2 to +0.4 from G to GHS, almost all of it environment. For unseen herds GH is no better than G; only the field helps. |

*Columns.* First column: the label of the expectation, used again in the results (§6.7). **Expectation:** what we
predicted, in words. "Masked cows" (test cows) are defined under *Cross-validation* in §4. "Unseen herds" in L8 were
herds with no cow among the training cows, a split that was later dropped (§4).

After the default-prior results showed collapsed fits, the **PC-prior refit** was planned and its expectations written
on **9 November 2026, before its results**:

| | Expectation |
|---|---|
| **PC1** | The collapsed fits disappear. |
| **PC2** | GHS at 1,000 cows rises from 0.035 to about 0.05–0.06. |
| **PC3** | If the prior causes the −18 % at 2,000 cows, the PC mean lands within 1 SE of 0.060. |
| **PC4** | WAIC is better with the PC prior in most fits; accuracy hardly changes. |

*Columns.* First column: the label of the expectation, used again in §6.7. **Expectation:** what we predicted for the
refit with PC priors (SE = standard error of the mean over the 10 replicates).

## 6. Results

### 6.1 Genetic variance

![GHS genetic variance per data set, default and PC prior](figures/fig1_ghs_genetic_variance.png)

**Variance components with INLA's default prior** (mean over data sets; 40 data sets at 1,000 and 2,000 cows, 5 at
4,000). The first row is the truth.

| Cows | Model | σ̂²g (bias) | 95 % CI covers 0.06 | Collapsed | h² | σ̂²h | σ̂²s | Range km | σ̂²e |
|---|---|---|---|---|---|---|---|---|---|
| **Truth** | | **0.060** | 95 % | 0 | 0.066 | 0.150 | 0.41 | 33 | 0.290 |
| 1,000 | G | 0.019 (−68 %) | 20 % | 31 | 0.022 | – | – | – | 0.832 |
| 1,000 | GH | 0.032 (−47 %) | 46 % | 21 | 0.038 | 0.505 | – | – | 0.323 |
| 1,000 | GS | 0.033 (−45 %) | 56 % | 17 | 0.034 | – | 0.516 | 35.2 | 0.456 |
| 1,000 | **GHS** | **0.035 (−41 %)** | 60 % | 15 | 0.038 | 0.151 | 0.475 | 39.7 | 0.319 |
| 2,000 | G | 0.038 (−37 %) | 57 % | 15 | 0.046 | – | – | – | 0.817 |
| 2,000 | GH | 0.047 (−22 %) | 82 % | 6 | 0.056 | 0.502 | – | – | 0.305 |
| 2,000 | GS | 0.047 (−22 %) | 90 % | 1 | 0.047 | – | 0.510 | 35.9 | 0.453 |
| 2,000 | **GHS** | **0.049 (−18 %)** | 92 % | 0 | 0.051 | 0.156 | 0.473 | 40.4 | 0.300 |
| 4,000 | G | 0.054 (−11 %) | 80 % | 0 | 0.063 | – | – | – | 0.804 |
| 4,000 | GH | 0.055 (−9 %) | 80 % | 0 | 0.065 | 0.511 | – | – | 0.293 |
| 4,000 | GS | 0.058 (−4 %) | 100 % | 0 | 0.062 | – | 0.433 | 37.7 | 0.435 |
| 4,000 | **GHS** | **0.057 (−5 %)** | 100 % | 0 | 0.062 | 0.151 | 0.421 | 41.0 | 0.288 |

*Columns.* **Cows:** cows per data set. **Model:** G genetic only; GH + herd; GS + spatial field; GHS all three (the
true model). **σ̂²g (bias):** estimated genetic variance, mean over the data sets (truth 0.060); in brackets the bias,
estimate ÷ 0.060 − 1. **95 % CI covers 0.06:** share of data sets whose 95 % credible interval for σ²g contains the
truth (target about 95 %). **Collapsed:** number of data sets with σ̂²g below 0.001, the bottom of the prior (target
0). **h²:** estimated heritability, σ̂²g ÷ the sum of all estimated variances (truth 0.066). **σ̂²h:** herd variance
(truth 0.150). **σ̂²s:** spatial-field variance (truth 0.41). **Range km:** distance at which two farms' environments
become almost unrelated (truth 33 km). **σ̂²e:** residual variance (truth 0.290). **–:** the model has no such term.

- **GHS is right only at 4,000 cows** (−5 %; sd over data sets 0.007). At 1,000 cows it is −41 % (sd 0.036) and at
  2,000 cows −18 % (sd 0.019).
- **At 1,000 cows the low mean comes from collapsed fits.** 15 of 40 GHS fits sit at about 0.0001 (G: 31 of 40). The 25
  GHS fits that did not collapse average 0.057, close to the truth. A collapsed fit sends the genetic variance to the
  residual (σ̂²e 0.358 against 0.296) and partly to the field.
- **At 2,000 cows nothing collapses, but the mean is 2.0 SE below the truth** (replicate means 0.034–0.079; SE of the
  mean 0.0053, from the 10 replicate means, because the four cells of a replicate share founders and field).
- Three fits returned no σ̂²g at all (GS n1000_p60 r9, GH n1000_p70 r5, GH n2000_p80 r7); they are counted as collapsed
  and left out of the means.

### 6.2 Where the variance of a missing term goes

![Variance components by model at 2,000 cows](figures/fig2_where_variance_goes.png)

- **G** puts herd and field into the **residual**: σ̂²e 0.82 ≈ 0.29 + 0.15 + 0.35.
- **GS** puts the herd into the **residual**: σ̂²e 0.45 ≈ 0.29 + 0.15.
- **GH** puts the field into the **herd**: σ̂²h 0.50 ≈ 0.15 + 0.35 (the part of the field that differs between herds).
- **No model puts the environment into the genetic variance.** With one record per cow a cow effect cannot hold anything
  that is not shared by relatives, and bulls are used at random in space, so relatives do not share an environment.
- **GHS gets the herd variance right** (0.151 and 0.156 against 0.150) and the field range slightly long (40 km against
  33; the posterior mean of a right-skewed parameter).

### 6.3 Accuracy and model fit

| Cows | Model | r(EBV, a) | Slope of a on EBV (median) | ΔWAIC vs GHS | Best WAIC in |
|---|---|---|---|---|---|
| **Target** | | higher | 1.00 | 0 | GHS in most |
| 1,000 | G | 0.354 | 0.70 | 684.7 | 0 of 40 |
| 1,000 | GH | 0.396 | 0.93 | 99.1 | 3 of 40 |
| 1,000 | GS | 0.430 | 0.94 | 184.8 | 0 of 40 |
| 1,000 | **GHS** | **0.440** | 0.94 | 0.0 | **37 of 40** |
| 2,000 | G | 0.383 | 0.93 | 1,372.7 | 0 of 40 |
| 2,000 | GH | 0.435 | 1.09 | 229.3 | 0 of 40 |
| 2,000 | GS | 0.472 | 1.31 | 360.1 | 0 of 40 |
| 2,000 | **GHS** | **0.484** | 1.24 | 0.0 | **40 of 40** |
| 4,000 | G | 0.410 | 1.01 | 2,794.8 | 0 of 5 |
| 4,000 | GH | 0.454 | 0.96 | 536.1 | 0 of 5 |
| 4,000 | GS | 0.487 | 1.03 | 650.7 | 0 of 5 |
| 4,000 | **GHS** | **0.498** | 1.00 | 0.0 | **5 of 5** |

*Columns.* **r(EBV, a):** accuracy, the correlation between each cow's estimated breeding value (EBV) and her true
breeding value (a), over all cows; it measures ranking, higher is better. **Slope of a on EBV (median):** regression
of the true value on the EBV, median over data sets; it measures the scale of the EBVs (target 1.00; above 1 = EBVs
shrunk too much, σ̂²g too low; below 1 = too spread). **ΔWAIC vs GHS:** the model's WAIC minus GHS's WAIC on the same
data set, mean (WAIC scores prediction, lower is better; positive = worse than GHS; above about 10 is clear).
**Best WAIC in:** in how many data sets the model had the lowest WAIC of the four.

- Accuracy ranks the models: G < GH < GS < GHS. **The field adds +0.04 to +0.05** over GH, and accuracy rises with the
  number of cows.
- **GHS has the best WAIC in 82 of 85 data sets.** The three exceptions are at 1,000 cows (in one of them GHS's σ̂²g had
  collapsed).
- The slope shows the scale error: at 2,000 cows GHS's σ̂²g is 18 % too low, so its EBVs are shrunk too much
  (slope 1.24). At 4,000 cows σ̂²g is right and the slope is 1.00.
- **A collapse ruins σ̂²g and h², not the ranking:** collapsed GHS fits still rank cows with r(EBV, a) 0.424, against
  0.450 in the others.

### 6.4 Cross-validation: cows predicted without their own records

Replicates 1–5. Every cow is a test cow once and is predicted from a fit without her records. Each measure is
computed within each test group and averaged over the 5 groups (§4).

**Accuracy of the test cows' EBVs, r(EBV, a), by share of single-cow herds**

| Cows | Model | 50 % | 60 % | 70 % | 80 % | All shares | Training cows (all shares) |
|---|---|---|---|---|---|---|---|
| 1,000 | G | 0.238 | 0.243 | 0.244 | 0.271 | 0.249 | 0.357 |
| 1,000 | GH | 0.267 | 0.260 | 0.265 | 0.269 | 0.265 | 0.386 |
| 1,000 | GS | 0.298 | 0.295 | 0.285 | 0.292 | 0.293 | 0.427 |
| 1,000 | **GHS** | **0.300** | **0.297** | **0.287** | **0.293** | **0.294** | 0.431 |
| 2,000 | G | 0.279 | 0.262 | 0.264 | 0.269 | 0.268 | 0.365 |
| 2,000 | GH | 0.309 | 0.291 | 0.299 | 0.293 | 0.298 | 0.405 |
| 2,000 | GS | 0.340 | 0.333 | 0.323 | 0.338 | 0.333 | 0.451 |
| 2,000 | **GHS** | **0.346** | **0.335** | **0.334** | **0.340** | **0.339** | 0.459 |

*Columns.* **Cows:** cows per data set. **Model:** as in §4. **50 % … 80 %:** the share of herds that have one cow;
*each cell* is the accuracy of the test cows, the correlation between the EBV each test cow got from the fit without her
records and her true breeding value, computed within each test group, averaged over the 5 groups of a data set and
then over the 5 data sets (replicates) with that share. **All shares:** the mean of the 20 data sets. **Training cows:**
the same accuracy for the cows whose records were used in the same fits, mean of the 20 data sets. Higher is better.

**Houaga's measure, cor(y, ŷ), by share of single-cow herds**

| Cows | Model | 50 % | 60 % | 70 % | 80 % | All shares | Genetic part (all shares) |
|---|---|---|---|---|---|---|---|
| 1,000 | G | 0.281 | 0.309 | 0.290 | 0.262 | 0.286 | 0.286 |
| 1,000 | GH | 0.502 | 0.490 | 0.441 | 0.405 | 0.459 | 0.299 |
| 1,000 | GS | 0.617 | 0.627 | 0.632 | 0.639 | 0.629 | 0.321 |
| 1,000 | **GHS** | **0.643** | **0.645** | **0.645** | **0.636** | **0.642** | 0.315 |
| 2,000 | G | 0.305 | 0.290 | 0.311 | 0.301 | 0.302 | 0.302 |
| 2,000 | GH | 0.511 | 0.502 | 0.487 | 0.434 | 0.483 | 0.315 |
| 2,000 | GS | 0.643 | 0.650 | 0.638 | 0.639 | 0.643 | 0.319 |
| 2,000 | **GHS** | **0.676** | **0.669** | **0.666** | **0.661** | **0.668** | 0.327 |

*Columns.* *Each cell* in the share columns is cor(y, ŷ), the correlation between each hidden record y and its
prediction ŷ (fixed effects + EBV + herd + field), computed within each test group and averaged as above. **All
shares:** the mean of the 20 data sets. **Genetic part:** cor(y, fixed effects + EBV), the part of the prediction
that comes from the fixed effects and the breeding value, mean of the 20 data sets. Higher is better; there is no true
value.

- **A cow without a record** gets an EBV of accuracy 0.29 (1,000 cows) to 0.34 (2,000 cows) from GHS, against 0.43–0.46
  for the training cows.
- **The share of single-cow herds hardly changes the test cows' accuracy** (GHS 0.287–0.300 with 1,000 cows and
  0.334–0.346 with 2,000), because a test cow's EBV comes from relatives spread over many herds.
- **Houaga's cor(y, ŷ) is mostly environment.** GHS reaches 0.64–0.68, but its genetic part, cor(y, fixed + EBV), is only
  0.29–0.33 in every model. A perfect genetic prediction (the true breeding value) would reach about 0.46–0.47.
- **Without a field, the environment of a cow alone in her herd cannot be predicted.** GH's cor(y, ŷ) falls as single-cow
  herds become more common (0.50 → 0.41 with 1,000 cows; 0.51 → 0.43 with 2,000), whereas GS and GHS stay at 0.62–0.68
  at every share: the field predicts that farm's environment from its neighbours.

### 6.5 By share of single-cow herds

Each cell is 10 replicates (cross-validation rows: 5). The 4,000-cow data sets exist only at 70 % and are not repeated
here.

**Genetic variance σ̂²g (bias), INLA's default prior**

| Cows | Model | 50 % | 60 % | 70 % | 80 % |
|---|---|---|---|---|---|
| **Truth** | | **0.060** | **0.060** | **0.060** | **0.060** |
| 1,000 | G | 0.007 (−89 %) | 0.040 (−34 %) | 0.014 (−77 %) | 0.016 (−74 %) |
| 1,000 | GH | 0.021 (−66 %) | 0.052 (−13 %) | 0.030 (−50 %) | 0.025 (−58 %) |
| 1,000 | GS | 0.024 (−60 %) | 0.051 (−14 %) | 0.027 (−56 %) | 0.032 (−46 %) |
| 1,000 | GHS | 0.027 (−55 %) | 0.049 (−19 %) | 0.029 (−52 %) | 0.037 (−38 %) |
| 2,000 | G | 0.048 (−20 %) | 0.029 (−52 %) | 0.032 (−47 %) | 0.043 (−29 %) |
| 2,000 | GH | 0.051 (−15 %) | 0.043 (−29 %) | 0.046 (−23 %) | 0.048 (−21 %) |
| 2,000 | GS | 0.049 (−18 %) | 0.039 (−35 %) | 0.043 (−28 %) | 0.054 (−9 %) |
| 2,000 | GHS | 0.053 (−12 %) | 0.045 (−24 %) | 0.046 (−23 %) | 0.053 (−12 %) |

*Columns.* **50 % … 80 %:** the share of herds that have one cow. *Each cell:* the estimated genetic variance, mean of
the 10 data sets with that share (truth 0.060); in brackets the bias, estimate ÷ 0.060 − 1.

**GHS, every level of the design**

| Cows | Measure | 50 % | 60 % | 70 % | 80 % |
|---|---|---|---|---|---|
| 1,000 | Collapsed GHS fits (of 10), default prior | 5 | 4 | 4 | 2 |
| 1,000 | σ̂²g with the PC prior (truth 0.060) | 0.055 | 0.070 | 0.055 | 0.060 |
| 1,000 | σ̂²h (truth 0.150) | 0.146 | 0.157 | 0.150 | 0.149 |
| 1,000 | r(EBV, a), GHS | 0.446 | 0.430 | 0.447 | 0.437 |
| 1,000 | r(EBV, a), GH | 0.407 | 0.386 | 0.414 | 0.377 |
| 1,000 | Accuracy gain from the field, GHS − GH | +0.039 | +0.044 | +0.033 | **+0.060** |
| 1,000 | Test-cow r(EBV, a), GHS (CV) | 0.300 | 0.297 | 0.287 | 0.293 |
| 2,000 | Collapsed GHS fits (of 10), default prior | 0 | 0 | 0 | 0 |
| 2,000 | σ̂²g with the PC prior (truth 0.060) | 0.060 | 0.054 | 0.055 | 0.060 |
| 2,000 | σ̂²h (truth 0.150) | 0.162 | 0.150 | 0.156 | 0.157 |
| 2,000 | r(EBV, a), GHS | 0.488 | 0.485 | 0.483 | 0.480 |
| 2,000 | r(EBV, a), GH | 0.450 | 0.438 | 0.436 | 0.415 |
| 2,000 | Accuracy gain from the field, GHS − GH | +0.038 | +0.047 | +0.047 | **+0.065** |
| 2,000 | Test-cow r(EBV, a), GHS (CV) | 0.346 | 0.335 | 0.334 | 0.340 |
| 2,000 | cor(y, ŷ), GHS / GH (CV) | 0.676 / 0.511 | 0.669 / 0.502 | 0.666 / 0.487 | 0.661 / 0.434 |

*Columns.* **Measure:** what the row shows (default prior unless stated). **50 % … 80 %:** the share of herds that
have one cow; each cell is the mean of 10 data sets (rows marked CV: 5). The measures:
- **Collapsed GHS fits:** fits with σ̂²g below 0.001.
- **σ̂²g with the PC prior:** the genetic variance when the same data are refitted with the PC prior.
- **σ̂²h:** the herd variance.
- **r(EBV, a):** accuracy, the correlation of estimated with true breeding values.
- **Accuracy gain from the field:** r(EBV, a) of GHS minus that of GH.
- **Test-cow r(EBV, a) (CV):** accuracy of cows whose records were hidden in cross-validation, computed within each
  test group and averaged (§4).
- **cor(y, ŷ), GHS / GH (CV):** Houaga's measure for the test cows with GHS and with GH, computed in the same way.

- **σ̂²g shows no trend with the share of single-cow herds.** Cells move by up to ±20 % in no order: with 10 data sets
  per cell and fits that collapse at random, a cell mean is noisy.
- **GHS gets the herd variance right at every share** (0.146–0.162).
- **The field's gain is larger at 80 % single-cow herds** (+0.060 and +0.065 against +0.033 to +0.047). Paired within
  replicate, 80 % minus the mean of 50–70 % is +0.021 (SE 0.007; larger in 16 of 20 pairs). *This was noticed after the
  results, not expected beforehand.* The fewer herd-mates a cow has, the less the herd term can do, and the more the
  neighbours matter.
- Herd and residual trade off as expected (L2): the correlation of σ̂²h and σ̂²e over data sets at 2,000 cows is −0.11,
  −0.54, −0.62 and −0.70 from 50 % to 80 %; at 1,000 cows −0.12 to −0.52 in no order.

### 6.6 Refit with PC priors (the same 80 data sets)

Each cell shows **default prior → PC prior** on identical data sets (the 4,000-cow data sets, already right, were not
refitted).

| Cows | Model | σ̂²g | Collapsed | CI covers 0.06 | r(EBV, a) | σ̂²e | WAIC better with PC |
|---|---|---|---|---|---|---|---|
| **Truth** | | **0.060** | 0 | about 95 % | – | 0.290 | – |
| 1,000 | G | 0.019 → 0.068 | 31 → 0 | 20 % → 98 % | 0.354 → 0.362 | 0.832 → 0.792 | 39 of 40 |
| 1,000 | GH | 0.032 → 0.063 | 21 → 0 | 46 % → 90 % | 0.396 → 0.403 | 0.323 → 0.291 | 40 of 40 |
| 1,000 | GS | 0.033 → 0.061 | 17 → 0 | 56 % → 95 % | 0.430 → 0.436 | 0.456 → 0.436 | 40 of 40 |
| 1,000 | **GHS** | **0.035 → 0.060** | 15 → 0 | 60 % → 92 % | 0.440 → 0.446 | 0.319 → 0.293 | 40 of 40 |
| 2,000 | G | 0.038 → 0.063 | 15 → 0 | 57 % → 92 % | 0.383 → 0.389 | 0.817 → 0.795 | 38 of 40 |
| 2,000 | GH | 0.047 → 0.060 | 6 → 0 | 82 % → 90 % | 0.435 → 0.437 | 0.305 → 0.293 | 40 of 40 |
| 2,000 | GS | 0.047 → 0.057 | 1 → 0 | 90 % → 98 % | 0.472 → 0.473 | 0.453 → 0.445 | 40 of 40 |
| 2,000 | **GHS** | **0.049 → 0.057** | 0 → 0 | 92 % → 98 % | 0.484 → 0.485 | 0.300 → 0.293 | 40 of 40 |

*Columns.* Each cell: **INLA's default prior → PC prior**, on the same 40 data sets per row. **σ̂²g:** estimated
genetic variance, mean (truth 0.060). **Collapsed:** data sets with σ̂²g below 0.001 (target 0). **CI covers 0.06:**
share of data sets whose 95 % credible interval contains the truth. **r(EBV, a):** accuracy, the correlation of
estimated with true breeding values. **σ̂²e:** residual variance (truth 0.290). **WAIC better with PC:** in how many of
the 40 data sets the PC fit had the lower (better) WAIC.

- **All 106 collapsed fits disappear**, and WAIC is better with the PC prior in **317 of 320 fits**.
- **GHS at 1,000 cows: 0.060**, exactly the truth on average. The 15 data sets that collapsed under the default prior
  give 0.034; the other 25 give 0.076.
- **GHS at 2,000 cows: 0.049 → 0.057.** Every one of the 40 estimates goes up (mean +0.008, SE 0.0005). With the
  replicate-based SE the default prior is 2.0 SE below the truth and the PC prior 0.6 SE: inside the 1-SE rule written
  beforehand. The prior explains about three quarters of the gap; the rest is sampling noise.
- **The prior changes the scale, not the ranking:** accuracy moves by +0.001 to +0.008; GHS's slope moves towards 1
  (1.24 → 1.12 at 2,000 cows). The extra genetic variance comes out of the residual.
- **No prior is neutral when data are weak:** G at 1,000 cows, the weakest model, overshoots to 0.068 (+13 %).

### 6.7 Expected against found

| | Expected (before the results) | Found | Verdict |
|---|---|---|---|
| **L1** | GHS unbiased (±15 %), imprecise | −41 % / −18 % / −5 % at 1,000 / 2,000 / 4,000 cows with the default prior; 0 % and −4 % with the PC prior | Holds only at 4,000 cows, or with the PC prior |
| **L2** | Herd–residual trade-off, worse at 80 % | Trade-off present; stronger at 80 % only at 2,000 cows | Partly |
| **L3** | G and GS put the environment into σ²g; GH puts the field into the herd | GH: yes (0.50). G and GS put it into the **residual** | Half wrong |
| **L4** | GHS accuracy 0.35–0.45 | 0.44 / 0.48 / 0.50 | Holds, at the top end and above |
| **L7** | Masked 0.15–0.30, training 0.35–0.6, rising with N | GHS test cows 0.29 → 0.34, training 0.43 → 0.46 | Holds (test cows at the top end) |
| **L8** | cor(y, ŷ) +0.2 to +0.4 from G to GHS, all environment; unseen herds GH ≈ G | +0.36 to +0.37, genetic part flat (0.29–0.33). The unseen-herd part held in the earlier pooled analysis (GH 0.31 ≈ G 0.30); that split is no longer reported | Holds |
| **PC1** | Collapses disappear | 106 → 0 | Holds |
| **PC2** | GHS at 1,000 cows about 0.05–0.06 | 0.060 | Holds |
| **PC3** | GHS at 2,000 cows within 1 SE of 0.060 | 0.057, 0.6 SE | Holds: the prior caused the −18 % |
| **PC4** | WAIC better, accuracy unchanged | WAIC better in 317 of 320; accuracy +0.001 to +0.008 | Holds |

*Columns.* First column: the label of the expectation (§5). **Expected:** what was written before the results.
**Found:** the result, in the units of §6. **Verdict:** whether the expectation held.

## 7. What step 1 shows

1. **With one record per cow the genetic variance is weakly identified, and INLA's default prior then decides the
   answer.** It collapses weak fits to about 0.0001 (15 of 40 GHS fits at 1,000 cows) and pulls moderately informative
   fits down (−18 % at 2,000 cows). With a PC prior the same data give unbiased estimates. Any analysis with few records
   per cow should use PC priors, or at least refit with a second prior and report how much the estimate moves.
2. **The model structure is not the problem in step 1.** With one record per cow no model puts the environment into
   the genetic variance; a missing herd or field goes to the residual (G, GS) or to the herd (GH). With PC priors all
   four models are within −5 % to +13 % of the true σ²g.
3. **The spatial field earns its place.** GHS has the best WAIC in 82 of 85 data sets, corrects the herd variance
   (0.15 instead of GH's 0.50), adds +0.04 to +0.05 accuracy for recorded cows and +0.03 to +0.04 for test cows, and
   keeps predicting records as single-cow herds become more common (cor(y, ŷ) 0.64–0.68 at every share), where GH falls
   to 0.41–0.43 at 80 %.
4. **The number of cows matters more than the share of single-cow herds** for the genetic variance: noisy at 1,000,
   modestly low at 2,000 (with the default prior), right at 4,000. The share of single-cow herds does not move σ̂²g,
   but the field helps most at 80 %.
5. **Report the accuracy of breeding values, not only cor(y, ŷ).** Houaga's validation measure is about twice the
   genetic accuracy and mostly reflects how well the environment is predicted.
6. **A collapse ruins the variance and heritability, not the ranking.** For choosing the best cows the prior hardly
   matters; for heritability, response to selection and index weights it matters a lot.

## 8. What step 1 cannot tell

- **One record per cow is not the real data**, which have about ten test-day records per cow. Repeated records and a
  permanent environment are steps 2 and 3. Their plan is in [PLAN_STEPS_2_3.md](PLAN_STEPS_2_3.md); their results are
  not in this repository. They change the picture, because a cow effect can then absorb a missing permanent environment.
- **Best-case world:** one breed, bulls used at random in space (no link between genetics and environment), and a field
  simulated from the same Matérn family that is fitted.
- **10 replicates per cell** (5 at 4,000 cows, which exists only at 70 %; no cross-validation at 4,000). Cell means can
  move by about ±20 % by chance; pooled rows are more reliable than single cells.
- **Empirical-Bayes integration** instead of INLA's default (checked on one data set only).
- The real genomic relationship matrix may differ from the simulated one (number of relatives, population structure).

## 9. Repository layout and how to run

<details>
<summary><b>Repository layout</b> (click to open)</summary>

```
README.md                         this file
REPRODUCE.md                      step-by-step replication, from the founder genomes to the tables, with checkpoints
run_step1.sh                      runs every step for any set of data sets on one machine (no Slurm)
docs/POPULATION.md                the simulated population: generations, families, herds, map, relatives close by
founders/                         the 10 founder-genome files of the paper (MaCS, about 2 MB each) and their MD5 sums
synthetic/backbone_synthetic.rds  synthetic stand-in for the real-data backbone (invented locations and records)
scripts/
  make_synthetic_backbone.R       builds synthetic/backbone_synthetic.rds
  questions/sb_simulate.R         simulate one data set (step 1: --records one --spe 0)
  questions/q_fit.R               fit G, GH, GS, GHS (--priors isidore = INLA's default, or --priors pc)
  questions/q_extract.R           compare the fits with the truth (bias, coverage, accuracy, slope, WAIC)
  questions/q_validate.R          5-fold cross-validation (--cv kfold)
  questions/sb_aggregate.R        pool all data sets of a study into tables
  sim/_paths.R, _mesh.R, _fit_setup.R   shared helpers (folders, Houaga's mesh recipe, data and field set-up)
  make_readme_figures.R           the three figures of this README, from results/
  describe_population.R           the numbers of docs/POPULATION.md (also rebuilds the pedigree of replicates 1-3)
cluster/
  q_task.sbatch                   Slurm array job: one task = one data set
  tasks_SB1.txt                   the 80 data sets at 1,000 and 2,000 cows
  tasks_SB1big.txt                the 5 data sets at 4,000 cows
  tasks_SB1pc.txt                 the PC-prior refit of the 80 data sets
results/
  step1_default_prior/            pooled output of sb_aggregate.R: report.md, sb_summary.csv (per cell x model),
                                  sb_by_rep.csv (per data set x model), sb_cv_by_rep.csv (cross-validation),
                                  sb_cv_testpop_by_rep.csv (cross-validation, within each test group: §4, §6.4)
  step1_pc_prior/                 the same for the PC refit, plus waic_default_vs_pc.csv (WAIC per fit, both priors)
  population/                     summary tables behind docs/POPULATION.md (counts and distances, no coordinates)
figures/                          PNG figures of this README
```

</details>

Each script starts with a comment block that explains its role in step 1, why it is needed, and what was expected.
Comments further down sometimes cite project notes (`docs/...`) that are not in this repository; this README holds
their step-1 content. The scripts also serve later steps and other studies of the project, so some options are not
used in step 1; the comment blocks say which.

**One data set, by hand** (from the repository root; first put the founder files in `runs/founders/` and a backbone in
`data/backbone/backbone_tz.rds`, see [REPRODUCE.md](REPRODUCE.md), Steps 1–2; `bash run_step1.sh` then runs all the
steps below):

```bash
RUN=runs_q/SB1_one_record/n1000_p70_r01
Rscript scripts/questions/sb_simulate.R --run $RUN --rep 1 --ncow 1000 --p1 0.70 --records one --spe 0
Rscript scripts/questions/q_fit.R      --run $RUN --qmodels G,GH,GS,GHS --priors isidore
Rscript scripts/questions/q_extract.R  --run $RUN
Rscript scripts/questions/q_validate.R --run $RUN --cv kfold --qmodels G,GH,GS,GHS --priors isidore   # replicates 1-5
Rscript scripts/questions/sb_aggregate.R --question SB1_one_record      # after all data sets are done
```

<details>
<summary><b>All data sets on a Slurm cluster</b> (click to open)</summary>

Edit the `cd` line of `cluster/q_task.sbatch` to your clone; create `logs/`:

```bash
# step 1, INLA's default prior: 80 data sets, cross-validation in replicates 1-5
SIMSCRIPT=scripts/questions/sb_simulate.R FITFLAGS="--priors isidore" VCV=kfold VMODELS=G,GH,GS,GHS VREPS=5 \
  sbatch --export=ALL,TASKS=cluster/tasks_SB1.txt --array=1-80%4 cluster/q_task.sbatch
# 4,000 cows: 5 data sets, no cross-validation (needs more memory and time)
SIMSCRIPT=scripts/questions/sb_simulate.R FITFLAGS="--priors isidore" \
  sbatch --export=ALL,TASKS=cluster/tasks_SB1big.txt --array=1-5%2 --mem=64G --time=7-00:00:00 cluster/q_task.sbatch
# PC-prior refit of the SAME data: copy data/ and truth/, then fit only
for d in runs_q/SB1_one_record/n1000_* runs_q/SB1_one_record/n2000_*; do
  mkdir -p runs_q/SB1_pc/$(basename $d) && cp -r $d/data $d/truth runs_q/SB1_pc/$(basename $d)/
done
SIMSCRIPT=scripts/questions/sb_simulate.R FITFLAGS="--priors pc" \
  sbatch --export=ALL,TASKS=cluster/tasks_SB1pc.txt --array=1-80%4 cluster/q_task.sbatch
Rscript scripts/questions/sb_aggregate.R --question SB1_pc
```

</details>

Set `SIMSCRIPT` as shown: the job script's default simulator belongs to another study and is not in this repository.
Seeds depend only on the replicate number (20261001 + 1000 × replicate), so a replicate gives the same founders,
bulls, cows and field in every cell. Founder genomes (MaCS) are cached in `runs/founders/` on first use.

**Cost** (median per data set, from the `minutes` columns of the results; 5 CPUs, 6 at 4,000 cows): fitting the four
models takes 8 min at 1,000 cows, 46 min at 2,000 and about 4 h at 4,000 (dense G); the 5-fold cross-validation adds
about 45 min at 1,000 cows and 4 h at 2,000.

## 10. Data not included

**Where the real structure comes from.** The herd structure and recording pattern of the simulation come from the
AADGG Tanzania data set of **Houaga et al. (2026)**, Isidore Houaga's file `milk_yield_pheno_cleaned.csv`: 1,894 cows,
1,386 herds, 19,375 test-day records. Our script `scripts/real/00a_build_backbone.R` (run once on a laptop; not in this
repository) turned it into the backbone below.

| Taken from that data set | How |
|---|---|
| Herd locations | GPS projected to km (EPSG:21035, as Houaga et al.), each moved by a random 1 km jitter; original herd IDs replaced |
| Herd sizes | sizes of the real multi-cow herds (single-cow herds are set by the design) |
| Region and ward of each herd | real membership, with new random ward IDs |
| Recording pattern of each cow | number of test-days, calendar month, calving season, lactation group, age, lactation-curve covariates; original cow IDs replaced |
| Fixed-effect sizes | one least-squares fit of the standardised milk yields on Houaga's fixed part; the milk records themselves are not kept |
| True variances of the simulation | chosen to match Houaga et al.'s published estimates (herd 0.15, spatial about 0.32 across herds, range 33 km) |

*Columns.* **Taken from that data set:** which part of the simulation's structure comes from the real data. **How:**
what was done to it, and what was removed or replaced.

The models, fixed part, mesh recipe and field prior follow Houaga's public script `1.Scripts_Spatial_Modelling.R`
(GitHub HighlanderLab/ihouaga_adgg_spatial). His genotypes and genomic relationship matrix were **not** used: the cows,
their genomes and G are simulated.

- **`data/backbone/backbone_tz.rds`** contains real farm locations and is therefore **not shared**. Sharing it would
  need the permission of Isidore Houaga and the owners of the AADGG data. It is an R list with:
  - `herds`: one row per real herd; columns `herd`, `region`, `ward`, `x_km`, `y_km` (projected coordinates in km);
  - `cows`: `animal_id`, `herd`, `dgrp` (breed class);
  - `records`: one row per real test-day record; columns `animal_id`, `herd`, `tyrmn`, `cyrsn`, `dgrp`, `lacgr`,
    `ageZ`, `leg1`, `leg2`, `year`, `Xb_true` (fixed-effect part of the record, from one least-squares fit to the real
    data), `rec_id`;
  - `b`: named fixed-effect coefficients from that fit (only the breed-class terms `dgrp2`–`dgrp4` are used, to remove
    the breed part);
  - `fixed_formula`, `meta`.

  No real milk record enters the simulated data: the simulated phenotypes are built from the simulated effects plus the
  fixed-effect sizes estimated once from the real milk yields.
- **A synthetic stand-in is included:** `synthetic/backbone_synthetic.rds`, made by `scripts/make_synthetic_backbone.R`.
  It has the same format, with invented herd locations and recording patterns. With it the whole pipeline runs end to
  end; the results are similar in kind but not identical to the numbers above (REPRODUCE.md, route B).
- **Per-data-set folders** (`runs_q/…`: simulated records, G and its inverse, fits) are not included: they are large and
  are rebuilt by the scripts. The pooled results in `results/` are what the README reports.

## 11. Software and references

**Software on the cluster:** R 4.5; R-INLA 26.8.7; inlabru 2.15.0; fmesher 0.8.0; AlphaSimR (development version,
commit a3d58f40); dqrng 0.4.1; ggplot2 and ragg for the figures.

**References**
- Houaga I. et al. (2026). *Genetics Selection Evolution* — genomic and spatial models for smallholder dairy data in
  Tanzania; the models, fixed effects, mesh recipe and field prior used here follow this paper.
- Simpson D., Rue H., Riebler A., Martins T. G., Sørbye S. H. (2017). Penalising model component complexity: a
  principled, practical approach to constructing priors. *Statistical Science* 32: 1–28.
- Lindgren F., Rue H., Lindström J. (2011). An explicit link between Gaussian fields and Gaussian Markov random fields:
  the stochastic partial differential equation approach. *Journal of the Royal Statistical Society B* 73: 423–498.
- Bachl F. E., Lindgren F., Borchers D. L., Illian J. B. (2019). inlabru: an R package for Bayesian spatial modelling
  from ecological survey data. *Methods in Ecology and Evolution* 10: 760–766.
- Gaynor R. C., Gorjanc G., Hickey J. M. (2021). AlphaSimR: an R package for breeding program simulations.
  *G3 Genes|Genomes|Genetics* 11: jkaa017.
- Chen G. K., Marjoram P., Wall J. D. (2009). Fast and flexible simulation of DNA sequence data. *Genome Research*
  19: 136–142.
- VanRaden P. M. (2008). Efficient methods to compute genomic predictions. *Journal of Dairy Science* 91: 4414–4423.
