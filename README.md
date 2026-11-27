# ACTG 175 HIV Trial
Reference: [Building and Checking Survival Models by David M. Rocke](https://dmrocke.ucdavis.edu/Class/EPI204-Spring-2021/Lecture12SurvivalAnalysisModelBuilding.pdf)

```
install.packages(c("speff2trial", "survival", "dplyr", "ggplot2", "survminer"))
library(speff2trial)
library(survival)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(survminer)
```
Prepare the data
```
data("ACTG175")
hiv <- ACTG175
View(hiv)
```
Clean variables from ACTG175 (n = 2,139 real HIV patients):
 - days: time to >50% decline in CD4, AIDS-defining event, or death (days)
 - cens: event indicator (1 = event occurred, 0 = censored)
 - arms: treatment arm (0 = ZDV monotherapy, 1 = ZDV+ddI, 2 = ZDV+Zal, 3 = ddI)
 - cd40: baseline CD4 T-cell count (cells/mm^3)
 - age: age in years
 - karnof: Karnofsky performance score (70-100)
 - race: 0 = White, 1 = Non-white

Relabel the arms as 1: ZDV_only, 2: ZDV_ddI, 3: ZDV_Zal, 4: ddI_only
```
hiv$arms <- factor(hiv$arms,
                   labels = c("ZDV_only", "ZDV_ddI", "ZDV_Zal", "ddI_only"))
hiv$karnof_grp <- factor(ifelse(hiv$karnof == 100, "Score_100", "Score_Below_100"))

cat("Patients:", nrow(hiv), "| Events:", sum(hiv$cens), "\n")
```
Output: Patients: 2139 | Events: 521 
Generate the survival curve: Kaplan-Meier
```
fit_km <- survfit(Surv(days, cens) ~ 1, data = hiv)
```
Check the events at some random time points
```
summary(fit_km, times = c(1:7,30,60,90*(1:10)))
```
Output:
<img width="424" height="299" alt="Screenshot 2026-11-18 170056" src="https://github.com/user-attachments/assets/ac9b0de3-ac5b-4127-9d96-3192ca3e4b60" />


Plot the Kaplan-Meier curve
```
plot(fit_km)
```
<img width="873" height="759" alt="image" src="https://github.com/user-attachments/assets/65fdb4e5-51b9-4dbe-8e79-26d2d0f052a5" />

Alternative plot with ggplot2 
```
autoplot(fit_km) + theme_bw() # theme_bw() is a predesigned "theme" which makes the plot prettier
```
<img width="873" height="759" alt="image" src="https://github.com/user-attachments/assets/75416e02-d0ba-454d-8d8a-90b0c026c17b" />


## 1. Naive Cox Model
```
fit_naive <- coxph(Surv(days, cens) ~ arms + cd40 + age + karnof_grp + race,
                   data = hiv)
summary(fit_naive)
```
<img width="499" height="332" alt="Screenshot 2026-11-18 170726" src="https://github.com/user-attachments/assets/f611c045-f664-4fdb-867f-145dac8011e5" />

<img width="440" height="173" alt="Screenshot 2026-11-18 170930" src="https://github.com/user-attachments/assets/76023335-a104-47e7-a781-808aad7b82fd" />


We estimate beta^ by maximizing the partial log-likelihood. The covariance matrix of beta^ is the inverse of the observed Fisher Information matrix: var(beta^) = I^(-1)(beta^)

### Observations
age, race have no effect (the confidence intervals of age and race contains 1).

All other variables have significant effects on the model.

Concordance is very high (0.67), which is a positive sign, likelihood is also very high, looks great.

ZDV_ddI: HR = 0.46 (95% CI: 0.36–0.59) — reduces risk by 54%.

ZDV_Zal: HR = 0.51 (95% CI: 0.40–0.65) — reduces risk by 49%.

ddI_only: HR = 0.56 (95% CI: 0.45–0.70) — reduces risk by 44%.

Baseline Functional Status (karnof_grp): Patients with a Karnofsky score below 100 have a 39% higher hazard (HR = 1.39, 95% CI: 1.17–1.65) compared to patients entering the trial at full functional capacity (score = 100).

cd40 has an HR of 0.9960, very close to zero, but still the confidence interval (0.9951, 0.9968) does not contain 1, which shows that cd40 has significant effect.
For example, we can scale the coefficient -0.004042 by some big number, say, 100, the hazard ratio becomes exp(-0.004042*100) = 0.668.
This shows that cd40 reduces the risk of progression or death by ~33%.

Even though age (p = 0.087) and race (p = 0.124) have confidence intervals crossing 1.0, they may appear as confounders.

It also can happen that age/race have a non-linear relationship (which Martingale residuals will reveal), which is causing these variables seem non-significant.

Note: Every p-value and confidence interval above assumes proportional hazards and linearity. Next, running cox.zph(fit_naive) will test whether those hazard ratios (especially arms and cd40) actually stay constant over the ~3-year follow-up or drift over time.


## 2. Schoenfeld Residuals: Test for Proportional Hazard (cox.zph)
```
zph_naive <- cox.zph(fit_naive, terms = TRUE)
print(zph_naive)
```
<img width="296" height="131" alt="Screenshot 2026-11-18 171250" src="https://github.com/user-attachments/assets/8d1fa69d-0f09-4f56-98b8-f075dc75426c" />
 
 - 'arms' and 'cd40' show strong non-proportionality (p < 0.01).
 - Global p < 0.001: Naive constant-HR Cox model violates PH assumption.

### Plot 1 (beta for arms), p = 0.03161: 


<img width="873" height="759" alt="image" src="https://github.com/user-attachments/assets/e641394c-70a2-4432-b4bd-8b1d7e03b67c" />

Downward Trend Over Time: 
 - starts high (around +2.3 at Day 150–200)
 - drops down to a plateau around +0.8 to +1.0 between Days 460 and 800
 - tends downward toward 0 by Day 1,000.
Because arms is a 4-level categorical factor (df = 3), so the plot combines the multi-degree-of-freedom term into a single linear predictor contrast (which is why the individual residual points split into distinct horizontal bands at the top and bottom).

### Plot 2 (beta for cd40), p = 0.00024:

<img width="873" height="759" alt="image" src="https://github.com/user-attachments/assets/d35119f3-38ed-4587-84ac-c0c1f50e2a37" />

 - From Day 200 to Day 460, the smoothed line sits steadily around −0.0065

 - Between Day 460 and Day 800, the curve climbs steadily upward to roughly −0.001, where the upper 95% confidence band touches 0.00 (meaning baseline CD4 count temporarily loses much of its predictive power among patients who survived past ~2 years), before dipping slightly again toward −0.004 near Day 990.

cd40 Had the Highest χ² (13.49) because cd40 is a continuous variable with tight standard errors (se = 0.00043), even a modest drift from −0.0065 to −0.001 is detected with extremely high statistical power (p = 0.00024).

Unlike arms (which formed horizontal bands), the dots here form a continuous cloud because cd40 takes on hundreds of distinct numerical values.

### Plot 3 (beta for age), p = 0.28211:

<img width="873" height="759" alt="image" src="https://github.com/user-attachments/assets/1acb8e47-3e71-4b9e-a600-d9e0575177f1" />

Aside from a tiny dip in the first ~300 days, the smoothed line runs almost flat right along 0.0 across the entire 1,000-day follow-up, with the 95% confidence band enclosing 0 throughout.

This visually confirms the non-significant cox.zph test (p = 0.282)—whatever small effect age has does not drift over time.

### Plot 4 (beta for karnof_grp), p = 0.88497

<img width="873" height="759" alt="image" src="https://github.com/user-attachments/assets/64c3d0ee-64ee-4ea0-87c7-e7b4da175d65" />

The true beta (the solid curve) stays around the estimated beta (= 0.3300611) across the entire 1,000-day trial, with zero slope.

Because karnof_grp is a binary variable (Score_Below_100 vs. Score_100), the individual residuals neatly split into two parallel horizontal clouds—one near +2.3 (events in the Score_Below_100 group) and one near −1.6 (events in the Score_100 group)—while the smoothed average sits cleanly in between.


### Plot 5 (beta for race), p = 0.97102

<img width="873" height="759" alt="image" src="https://github.com/user-attachments/assets/d079d5db-8812-4df9-b007-0ed492ce5785" />

True beta sits right around −0.16 (matching coef = -0.158) with nearly zero slope across the entire follow-up, and the 95% confidence band includes 0 everywhere.

Just like karnof_grp, the residuals split into two parallel horizontal clouds—with a denser bottom band: 0 = White, 1 = Non-white, White: 71.2%, black: 28.8%, who represent the majority of events) and a sparser top band (Non-white patients).

We want to find whether our assumption that there is a time-varying term g(t) hiding inside the true beta holds. If the covariate X_j is time-varying, and we only consider its subject-varying trend, we ignore its time component; therefore, its coefficient beta_j in the regression model must contain that time component. So, we hope that even if we ignore that time-varying part of X_j, it should be present in the Schoenfeld residuals of X_j (the raw X's used in the regression model minus their risk-weighted average values). So this can be reflected in the (Pearson) correlation coefficient between the hypothesized time-varying component and the raw observed time-varying component- the residuals: that is, t_j. But we do not use the raw Schoenfeld residuals; instead, we scale the residuals by the covariance matrix of beta^ at time t.

<img width="155" height="53" alt="Screenshot 2026-11-27 170552" src="https://github.com/user-attachments/assets/d23b9c96-9ced-42c8-88bb-1449cfae82fc" />

<img width="518" height="341" alt="Screenshot 2026-11-27 164202" src="https://github.com/user-attachments/assets/10c2923a-ad3e-4777-bb04-386524b0e93b" />


Further, let us examine the proportional hazards criterion for the categorical variables arms and karnof_grp by plotting the log(-log(survival)) vs log(time) that is log(cumulative hazard) vs log(time).


ggsurvplot() expects a Kaplan–Meier fit created by survfit() as its first argument.
```
km_arms <- survfit(Surv(days, cens) ~ arms, data = hiv)
```
Plot for 'arms':
```
ggsurvplot(km_arms, data = hiv, fun = "cloglog",
           xlab = "Log(Time in Days)",
           ylab = "Log(-Log(Survival))",
           title = "Log-Cumulative Hazard: Treatment Arms")

```
<img width="1462" height="759" alt="image" src="https://github.com/user-attachments/assets/89da3f4d-560c-46d1-b720-5db5644afcfc" />

### Observations: 

Log-Cumulative Hazard (cloglog) plot for Treatment Arms visually confirms why arms failed the proportional hazards test (cox.zph p = 0.0316): the four curves cross early in follow-up (before Day 200) and gradually narrow toward the end of the trial.

From Day 200 onward, the three active regimens (ZDV_ddI, ZDV_Zal, ddI_only) merge into a single overlapping bundle, whereas ZDV_only (red) stays distinctly higher throughout the rest of the trial. So all these three alternative regimens achieved a similar level of protection; their cumulative hazard curves overlap tightly from Day 200 onward.

This happens because ZDV monotherapy is clinically inferior to the other three treatments—patients on ZDV_only experience disease progression or death at roughly twice the rate (HR ≈ 0.46 to 0.56 for the three active arms vs. ZDV_only). On a cloglog plot, a higher curve means a higher cumulative hazard (worse survival).


 - Days 200-300: the vertical gap is wide (~1.0 log-hazard units apart).
 - Days 1,000+: however, that vertical gap shrinks to about ~0.5 log-hazard units—meaning ZDV_only slowly converges partway toward the active bundle late in follow-up.



```
km_karnof <- survfit(Surv(days, cens) ~ karnof_grp, data = hiv)

```
Plot for 'karnof_grp':
```
ggsurvplot(km_karnof, data = hiv, fun = "cloglog",
           xlab = "Log(Time in Days)",
           ylab = "Log(-Log(Survival))",
           title = "Log-Cumulative Hazard: Karnofsky Group")
```

<img width="1462" height="759" alt="image" src="https://github.com/user-attachments/assets/8e37b8dd-8fe1-4ca1-9a8f-cd1c904dfcdb" />

### Observations:

The two curves for the two Karnofsky groups never cross and stay virtually parallel across the entire follow-up period - proportional hazards assumption is satisfied.


## 3. Martingale Residuals: Check Non-linearity of continuous predictors

For each time interval, the martingale residual measures the difference between what actually happened (whether the event occurred in that slice, 0 or 1) and what the model expected to happen based on the patient's covariates and strata.
```
fit_null  <- coxph(Surv(days, cens) ~ 1, data = hiv)
mart_null <- residuals(fit_null, type = "martingale")
```

### Martingale Residuals vs Baseline cd40

```
plot(hiv$cd40, mart_null, xlab = "Baseline cd40 Count", ylab = "Martingale Residuals",
     main = "Martingale Linearity Check: Baseline cd40", pch = 16, col = rgb(0, 0, 0, 0.25))
lines(lowess(hiv$cd40, mart_null), col = "red", lwd = 2)
```

<img width="873" height="759" alt="image" src="https://github.com/user-attachments/assets/f3160d62-4b69-4317-9f55-75b9a07ed667" />


A clear "elbow" pattern is visible: risk drops sharply as cd40 increases from 0 to ~250–300, then completely flattens out above 300.

So the variable cd40 needs a spline (pspline(cd40)) or log/sqrt transformation.


### Martingale Residuals vs Age

```
plot(hiv$age, mart_null, xlab = "Age (years)", ylab = "Martingale Residuals",
     main = "Martingale Linearity Check: Age", pch = 16, col = rgb(0, 0, 0, 0.25))
lines(lowess(hiv$age, mart_null), col = "blue", lwd = 2)
```
<img width="873" height="759" alt="image" src="https://github.com/user-attachments/assets/04821428-109d-4677-9ae3-2d9c364fb206" />


The plot shows a flat, linear trend across the entire age range (10 to 70 years), confirming that age has no hidden non-linear curve and requires no transformation.
Earlier, we wondered whether age looked statistically non-significant (p = 0.087) because of a hidden non-linear curve (like risk spiking only in older patients). This plot rules that out—age genuinely has a weak, steady linear relationship with the outcome in this HIV cohort, so keeping it as a simple linear term (age) is completely appropriate.



### 4. Deviance resudual & dfbeta: Outliers & influencial observations

Deviance residuals measure how poorly the model predicted each patient's outcome, reshaped so that a well-fitting model scatters evenly between roughly −2.5 and +2.5.

Just like regular residuals in linear regression, values close to 0 mean the model's prediction matched what actually happened to that patient.
Positive values (> 0): The patient experienced the event earlier than the model expected (the model was too optimistic).
Negative values (< 0): The patient survived longer than the model expected (the model was too pessimistic).
Raw Martingale residuals are squished on top (capped at +1) and stretched on the bottom (down to -∞), making it hard to spot outliers. Deviance residuals fix that lopsidedness so both positive and negative outliers stand out on the same balanced scale.

Deviance Residuals (Outliers)
```
dev_res <- residuals(fit_naive, type = "deviance")
plot(1:nrow(hiv_real), dev_res, xlab = "Patient Index", ylab = "Deviance Residuals",
     main = "Deviance Residuals (Outlier Detection)", pch = 16, col = "steelblue")
abline(h = c(-2.5, 2.5), col = "red", lty = 2)
```
<img width="873" height="759" alt="image" src="https://github.com/user-attachments/assets/afccd0de-f44c-4463-bb94-6f1c5aa2dabb" />

### Observations

Deviance residuals are balanced mostly within ±2.5, with only ~20–25 mild positive outliers (~1% of n=2,139, expected by chance), confirming no extreme observation anomalies.
Overall strong fit. The vast majority of the 2,139 patients sit comfortably inside the red dashed lines (−2.5 to +2.5), with no bizarre spikes (like +5 or −5) that would point to corrupted data or severe model failure.
Normal count of mild outliers: Around 20–25 patients poke slightly above +2.5 (up to ~3.3). With 2,139 people in the study, having ~1% of patients fall just past 2.5 is completely expected by pure chance—these are simply a handful of patients who looked healthy on paper (high CD4, good treatment arm) but suffered an unexpectedly early event.
Zero negative outliers: Nobody falls below −2.5, meaning there are no high-risk patients who survived absurdly long in defiance of the model.


dfbeta (Influence on Baseline CD4 coefficient):

In the dfbeta plot, each dot shows exactly how much the regression coefficient beta^ would change if we removed that single patient and refit the mode beta^ - beta^_(-1)
```
dfbeta_res <- residuals(fit_naive, type = "dfbeta")
colnames(dfbeta_res) <- names(coef(fit_naive))
plot(1:nrow(hiv), dfbeta_res[, "cd40"], xlab = "Patient Index",
     ylab = "dfbeta (Baseline CD4)", main = "Observation Influence (dfbeta) on CD4",
     pch = 16, col = "darkorange")
abline(h = 0, lty = 2)
```
<img width="873" height="759" alt="image" src="https://github.com/user-attachments/assets/5c647dd0-97f5-4be0-a002-6558fc42fb9f" />


Max dfbeta for cd40 is ~0.00008 (<2% of coef -0.00404 and <0.2 SE), confirming no single patient unduly influences the model.



## Fixing Violations: Stratification + Spline + Time Interaction (cd40 * log(t))
### Step A: Stratify by 'karnof_grp' and 'race' to absorb confounding baseline hazards and model non-linear 'cd40' using a penalized spline pspline(cd40).

```
fit_strat <- coxph(Surv(days, cens) ~ arms + pspline(cd40) + age + strata(karnof_grp, race),
                   data = hiv)
print(cox.zph(fit_strat))
```
<img width="199" height="74" alt="Screenshot 2026-11-18 173444" src="https://github.com/user-attachments/assets/cf90f6dc-075b-4d3d-bca6-ac8749e8e8e1" />



Model Summary (Stratified + Penalized Spline for CD4):

Stratifying by Karnofsky score and race (strata(karnof_grp, race)) absorbs their
baseline hazard differences, while pspline(cd40) accounts for the non-linear L-shaped
risk curve of baseline cd40 identified in the Martingale residuals. This improves the
global proportional hazards test (GLOBAL chisq drops from 24.72 to 20.71, p = 0.008)
and substantially reduces the CD4 violation (p improves from 0.00024 to 0.039).
However, both arms (p = 0.036) and pspline(cd40) (p = 0.039) still show mild
non-proportionality (p < 0.05), confirming that accounting for non-linearity and
stratification alone is not enough—we must also add time-varying interactions (tt()).


Split follow-up at every observed event time so we can model cd40 * log(time) explicitly

```

hiv_split <- survSplit(Surv(days, cens) ~ ., data = hiv,
                       cut = event_times, episode = "tgroup")

```


Create the explicit time-varying interaction term: cd40 * log(t)
```
hiv_split$cd40_logt <- hiv_split$cd40 * log(hiv_split$days)

fit <- coxph(Surv(tstart, days, cens) ~ arms + cd40 + cd40_logt + age + strata(karnof_grp, race),
                   data = hiv_split,
                   id = patient_id)
summary(fit)

```
<img width="426" height="367" alt="Screenshot 2026-11-18 174536" src="https://github.com/user-attachments/assets/2be5a67c-6105-4742-9015-a5835fd99594" />


Final Time-Varying Stratified Cox Model Summary (fit):

 - Time-Varying CD4 Effect: Both the baseline cd40 main effect (coef = -0.0216, p = 7.59e-06) and its log-time interaction tt(cd40) (coef = +0.00281, p = 0.000245) are highly significant. Because the interaction is positive, the strong early protective effect of higher baseline CD4 counts diminishes over time (e.g., for +100 CD4 cells, HR = 0.36 at Day 60 vs. HR = 0.72 at
Day 700), resolving the non-proportionality observed in cox.zph().
 - Treatment Efficacy: After stratifying by Karnofsky score and race and accounting for time-varying CD4 risk, all alternative regimens significantly outperform ZDV monotherapy (p < 0.0001): ZDV+ddI reduces hazard by 53.8% (HR = 0.46, 95% CI: 0.36-0.59), ZDV+Zal by 48.8% (HR = 0.51, 95% CI: 0.40-0.65), and ddI alone by 43.9% (HR = 0.56, 95% CI: 0.45-0.70).
 - Age & Model Fit: Age remains non-significant (HR = 1.008 per year, p = 0.099).
Overall model fit is highly significant (Likelihood Ratio Test = 159.2 on 6 df, p < 2e-16; Concordance = 0.654).

## 6. Cox-Snell Residuals: Overall Goodness of Fit (Fixed for tt() interval expansion) 



<img width="419" height="307" alt="Screenshot 2026-11-18 192959" src="https://github.com/user-attachments/assets/f87c005c-7427-4df3-b977-a05b16c6cf12" />

<img width="416" height="307" alt="Screenshot 2026-11-18 193551" src="https://github.com/user-attachments/assets/eacf6a54-9c03-4c6a-9ec1-60741b2641cf" />

<img width="414" height="301" alt="Screenshot 2026-11-18 193627" src="https://github.com/user-attachments/assets/83917650-db8f-4fd2-9275-80a50fdcca7c" />


In standard linear regression, a residual e_i = y_i - y^_i is useful because if the model is right, the residuals follow a known, universal reference distribution—specifically N(0, sigma^2) —regardless of what the predictors X_i were. We check model fit by verifying whether the residuals look like a sample from that reference distribution.
In survival analysis, event times T_i are positive, right-skewed, and often censored, so simple subtraction t_i - t^_i does not produce a clean reference distribution. However, transforming t_i through the estimated cumulative hazard function r_Ci = H^(t_i|X_i) strips away the effect of the covariates X_i and the baseline hazard shape. If the model is correct, every subject's r_Ci is mapped onto the same universal reference distribution: Exponential(1).
Because r_Ci plays the exact diagnostic role of a regression residual (standardizing observations so they can be checked against a known theoretical distribution), Cox and Snell named it a generalized residual.

Note that subtracting this from the event indicator delta_i gives the Martingale residual M_i = delta_i – r_Ci, which has mean 0 and looks even more like a traditional "observed minus expected" residual.

When learning that r_Ci = H^(t_i|X_i) ~ Exponential(1), the most instinct is to plot a histogram of r_Ci and overlay the Exponential(1) density curve f(r) = e^(-r) to see if they match. 

But a histogram treats every number in r_Ci as a completed observation. If a subject was censored early (say at r_Ci = 0.2), their true value might have been 1.8, but the histogram puts a bar at 0.2. Because censored subjects are artificially cut short, the histogram piles up near 0 and drops off too fast, making a completely valid model look broken. To fix this, we need an estimator that accounts for right-censoring—which is why we use the Nelson-Aalen estimator instead of a histogram.

Instead, we use the fact that: the cumulative hazard of the cumulative hazard is the identity function.

 - The x-axis r_Ci: This is the model-based cumulative hazard for subject i at their observed time t_i, calculated using our fitted parametric or Cox model: r_Ci = H^(t_i|X_i).

 - The y-axis H^(r_Ci) = H^( H^(t_i|X_i)): Here, we temporarily pretend that r_Ci is a "survival time" and apply the non-parametric Nelson-Aalen estimator to the pairs (r_Ci, delta_i). This estimates the cumulative hazard of the residuals themselves directly from the data, properly adjusting for which r_Ci values were censored.

<img width="518" height="158" alt="Screenshot 2026-11-27 180846" src="https://github.com/user-attachments/assets/6c34c0e1-96b7-47d7-98af-34e8e1ad36bd" />


---

Calculate the martingale residual M_ij for every single time interval.

Because coxph() automatically drops any rows with missing values (NA), extracting residuals first ensures we only work with the exact rows the model used.

```
mart_interval <- residuals(fit, type = "martingale")
```
Extract the exact patient_id and event status from the fitted model's data
model.frame(fit_fixed) grabs the internal data table stored inside the fitted model object fit, and pulls out the patient IDs and the cens (status) column.

```
model_frame <- model.frame(fit)
head(model_frame)
```

<img width="536" height="119" alt="Screenshot 2026-11-18 201933" src="https://github.com/user-attachments/assets/b8a68087-0527-43f6-9252-a96829f6cb78" />

A new column called id has been added; this is the subject ID

```
used_subjid <- mf[["(id)"]]
used_event <- mf[["Surv(tstart, days, cens)"]][, "status"]
```

Pulling used_subjid and used_event directly from the fitted model guarantees that every residual lines up 1-to-1 with its true patient and event status.

Collapse both to one value per patient_id

Sum the mart_intervals to get the patient's total martingale residual across their entire follow-up: M_i = sum_j(M_ij).

```
mart_id <- tapply(mart_interval, used_subjid, sum)
```

Summing used_event gives the patient's final overall status: 1 if they experienced the event at the end of follow-up, or 0 if they were censored throughout (since earlier intervals all have cens = 0).

```
event_id <- tapply(event, used_id, sum)

```
Compute Cox-Snell residuals: r_Ci = delta_i - M_i
```
cs_residual <- event_id - mart_id
```

Estimate cumulative hazard and plot against the x = y line.

For a true exponential distribution with mean: Exponential(1), the hazard rate is constant at lambda(t) = 1, which means its cumulative hazard function is simply H(t) = t (a straight line with intercept 0 and slope 1. Comparing H^(r_Ci) visually against the 45-degree line lets us verify whether the residuals actually behave like an Exponential(1) sample and therefore whether the Cox model adequately describes the data.

```
fit_cs <- coxph(Surv(cs_residual, event_id) ~ 1, ties = "efron")
```


basehaz() extracts the Nelson-Aalen cumulative hazard estimate H^(CS) from the null model.
na_haz <- basehaz(fit_cs, centered = FALSE)

Finally, let us plot a step function of the estimated cumulative hazard against the residual values.

```
plot(na_haz$time, na_haz$hazard,
     type = "s",
     xlab = "Cox-Snell Residuals",
     ylab = "Estimated Cumulative Hazard (Nelson-Aalen)",
     main = "Cox-Snell Residual Goodness-of-Fit Plot")
abline(0, 1, col = "red", lty = 2, lwd = 2)
```

<img width="877" height="759" alt="image" src="https://github.com/user-attachments/assets/28aac544-c76a-4824-b23d-504ccd40dbe8" />


### Observations:
 - Lower-to-Mid Range (0.0 to 0.5): The empirical cumulative hazard aligns almost seamlessly with the theoretical y = x line, meaning the model predicts risk accurately where the bulk of the observations and events occur.
 - Upper-Middle Range (0.5 to 0.8): There is a very slight downward dip below the red line, which is minor and well within normal sampling variability rather than a systematic model failure.
 - Upper Tail (0.8 to 1.1): The larger vertical jumps near 0.85 and 1.0 are expected: at large residual values, fewer patients remain at risk due to prior events and censoring, so single events create larger step increments in the Nelson-Aalen estimator.

