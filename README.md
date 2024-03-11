# riskModel

The *riskModel* pkg is developed based on several publicly available packages and codes for developing risk prediction models using survival data.

# Collections of risk prediction models

## BMJ series papers

* **2024** 
  * [Part 1. Evaluation of clinical prediction models: from development to external validation](https://www.bmj.com/content/384/bmj-2023-074819), [R code]( https://github.com/gscollins1973/validationCRASH)
  * [Part 2. Evaluation of clinical prediction models: how to undertake an external validation study](https://www.bmj.com/content/384/bmj-2023-074820)
  * [Part 3. Evaluation of clinical prediction models: calculating the sample size required for an external validation study](https://www.bmj.com/content/384/bmj-2023-074821), [R code](https://www.prognosisresearch.com/software)

* **2023**
  * [Transparent reporting of multivariable prediction models developed or validated using clustered data (TRIPOD-Cluster): explanation and elaboration](https://www.bmj.com/content/380/bmj-2022-071058)
  * [Transparent reporting of multivariable prediction models developed or validated using clustered data: TRIPOD-Cluster checklist](https://www.bmj.com/content/380/bmj-2022-071018)

* **Before 2023**
  * 2022 [Validation of prediction models in the presence of competing risks: a guide through modern methods](https://www.bmj.com/content/377/bmj-2021-069249), [R code](https://github.com/survival-lumc/ValidationCompRisks)
  * 2020 [Calculating the sample size required for developing a clinical prediction model](https://www.bmj.com/content/368/bmj.m441), [R code](https://www.bmj.com/content/368/bmj.m441)
  * 2019 [Guide to presenting clinical prediction models for use in clinical settings](https://www.bmj.com/content/365/bmj.l737), [R code](https://cran.r-project.org/web/packages/rms/index.html)
  * 2016 [Adding tests to risk based guidelines: evaluating improvements in prediction for an intermediate risk group](https://www.bmj.com/content/354/bmj.i4450)
  * 2015 [How to develop a more accurate risk prediction model when there are few events](https://www.bmj.com/content/351/bmj.h3868)
  * 2012 [Comparing risk prediction models](https://www.bmj.com/content/344/bmj.e3186)
  * 2012 *Heart* [Part 1. Risk prediction models: I. Development, internal validation, and assessing the incremental value of a new biomarker](https://heart.bmj.com/content/98/9/683)
  * 2012 *Heart* [Part 2. Risk prediction models: II. External validation, model updating, and impact assessment](https://heart.bmj.com/content/98/9/691)


## Methodological papers

* 2022 *Ann Intern Med* [Assessing Performance and Clinical Usefulness in Prediction Models With Survival Outcomes: Practical Guidance for Cox Proportional Hazards Models](https://www.acpjournals.org/doi/10.7326/M22-0844), [R code](https://github.com/danielegiardiello/Prediction_performance_survival)
* 2016 *J Clin Oncol* [Assessing the clinical impact of risk prediction models with decision curves: guidance for correct interpretation and appropriate Use](https://ascopubs.org/doi/10.1200/JCO.2015.65.5654), [R code](https://cran.r-project.org/web/packages/dcurves/)
* 2014 *Ann Intern Med* [Net reclassification improvement: compuataion, interpretation, and controversies](https://www.acpjournals.org/doi/10.7326/M13-1522)
* 2010 *Epidemiology* [Assessing the performance of prediction models: A framework for traditional and novel measures](https://journals.lww.com/epidem/fulltext/2010/01000/assessing_the_performance_of_prediction_models__a.22.aspx), [R code 1](https://links.lww.com/EDE/A355), [R code 2](https://github.com/danielegiardiello/ValLogRegMod)
* 2003 *J Clin Epidemiol* [Internal and external validation of predictive models: a simulation study of bias and precision in small samples](https://linkinghub.elsevier.com/retrieve/pii/S0895435603000477)
* 2001 *J Clin Epidemiol* [Internal validation of predictive models: efficiency of some procedures for logistic regression analysis](https://linkinghub.elsevier.com/retrieve/pii/S0895435601003419)


## Advanced topics

### Predictors selections

* 2023 *J Clin Epidemiol* [Practical guide to the typical analysis of prognostic factors and biomarkers without the use of P-values](https://linkinghub.elsevier.com/retrieve/pii/S0895435623000768)

### Super learner

* 2018 *Epidemiology* [Using super learner prediction modeling to improve high-dimensional propensity score estimation](https://journals.lww.com/epidem/fulltext/2018/01000/using_super_learner_prediction_modeling_to_improve.13.aspx), [R code 1](https://github.com/lendle/hdps); [R code 2](https://github.com/lendle/TargetedLearning.jl)

### Dynamic prediction

* 2020 *Int J Epidemiol* [Reflections on modern methods: Dynamic prediction using joint models of longitudinal and time-to-event data](https://academic.oup.com/ije/article/50/5/1731/6174516?login=false), [R code](https://erandrinopoulou.github.io/EducationalCorner_JMpred/index.html)

* 2015 *Int J Epidemiol* [Joint modelling of repeated measurement and time-to-event data: an introductory tutorial](https://academic.oup.com/ije/article/44/1/334/657852?login=false), [R code](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/ije/44/1/10.1093_ije_dyu262/4/dyu262_Supplementary_Data.zip?Expires=1713186801&Signature=R3cLsdcvgZ0cJAxnpDi8ZYoKkAKxUGEF4JSAvtG3GmdKQ4nn0LHv30BIkmQHYJUCnyrgNUyx2nUthGP8W4hcsa9Xi6Hr9T5I1Lq6HbaIUkFQlIcjKU0Gl1asnr7bUwF7M23sdUjvymzsc1HAgTqUgMSgk1g2CrEWp8j1uVrexd3StsaRR6hBPSONYO72mZO3IVfpQ7en0fVo8nUVwatMku28rZaHFm7dmhV3gi6CrblZ0jQf~i6kLrMmVGYZBKSGNyQ7E4~Bw~zyzxHCV-rOpOoiAdLV8q7Bhv6p20-CLcLHP9IZUlHNmDlF-6SN0nmAOOjWNTl8HjboweaovJY7lQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

* 2013 *Biometrics* [Dynamic pseudo-observations: a robust approach to dynamic prediction in competing risks](https://academic.oup.com/biometrics/article/69/4/1043/7492355?login=false), [R code](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/biometrics/69/4/10.1111_biom.12061/5/biom12061-sm-0001-supinfo-s1.pdf?Expires=1713186954&Signature=1edWKfgynE-cL-7qyr7fHvS8izB29WPN89pDi4qSk5548tyLZbFfxUnoT1BwDuB~RqS8ukHv2-6lw2FqZguU60uMrY8D-lZ02xalQWjPd0LmYlw1SUblFbetY58v~6tPE2~pcFspHXUnGaN8rZVGuiSpv4Acb4Vkpp4bPtqNh8Xicu376oiXCV9YRWmWVggda3jx5gQ9Sp2rG3OTCkjXsVrL1jxVltdqeyBOBoUpuStSpFK-r-~4c5ChijDim49y4LD2OxNYavA0C0hDeqdSt9hnLJDiir0AzZKNF4yU9YrMC5TMGyL0PvQaE3oOm7rTIPqu14FkSjWZBNlZUZHTcA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

### Prediction of counterfactual risk

* 2022 *Eur J Epidemiol* [Predicting counterfactual risks under hypothetical treatment strategies: an application to HIV](https://link.springer.com/article/10.1007/s10654-022-00855-8)

* 2022 *Eur J Epidemiol* [A structural characterization of shortcut features for prediction](https://link.springer.com/article/10.1007/s10654-022-00892-3)

* 2020 *Eur J Epidemiol* [Counterfactual prediction is not only for causal inference](https://link.springer.com/article/10.1007/s10654-020-00659-8)

* 2020 *Eur J Epidemiol* [Prediction meets causal inference: the role of treatment in clinical prediction models](https://link.springer.com/article/10.1007/s10654-020-00636-1)

### Lifetime risk

* 2021 *Eur J Epidemiol* [A comparison of statistical methods to predict the residual lifetime risk](https://link.springer.com/article/10.1007/s10654-021-00815-8), [R code](https://github.com/s-conner/lifetimerisk)

* 2000 *Stat Med* [Computing estimates of incidence, including lifetime risk: Alzheimer’s disease in the framingham study. the practical incidence estimators (pie) macro](https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1097-0258(20000615/30)19:11/12%3C1495::AID-SIM441%3E3.0.CO;2-E)



