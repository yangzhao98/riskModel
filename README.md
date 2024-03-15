# riskModel

The *riskModel* pkg is developed based on several publicly available packages and codes for developing risk prediction models using survival data.

# Collections of risk prediction models

## 1. BMJ series papers

* **2024** 
  * [Part 1. Evaluation of clinical prediction models: from development to external validation](https://www.bmj.com/content/384/bmj-2023-074819), [R code]( https://github.com/gscollins1973/validationCRASH)
  
  * [Part 2. Evaluation of clinical prediction models: how to undertake an external validation study](https://www.bmj.com/content/384/bmj-2023-074820)
  
  * [Part 3. Evaluation of clinical prediction models: calculating the sample size required for an external validation study](https://www.bmj.com/content/384/bmj-2023-074821), [R code](https://www.prognosisresearch.com/software)

* **2023**
  * [Transparent reporting of multivariable prediction models developed or validated using clustered data (TRIPOD-Cluster): explanation and elaboration](https://www.bmj.com/content/380/bmj-2022-071058)
  
  * [Transparent reporting of multivariable prediction models developed or validated using clustered data: TRIPOD-Cluster checklist](https://www.bmj.com/content/380/bmj-2022-071018)

  * [Development and internal-external validation of statistical and machine learning models for breast cancer prognostication: cohort study](https://www.bmj.com/content/381/bmj-2022-073800)

* **Before 2023**
  * 2022 [Validation of prediction models in the presence of competing risks: a guide through modern methods](https://www.bmj.com/content/377/bmj-2021-069249), [R code](https://github.com/survival-lumc/ValidationCompRisks)
  
  * 2020 [Calculating the sample size required for developing a clinical prediction model](https://www.bmj.com/content/368/bmj.m441), [R code](https://www.bmj.com/content/368/bmj.m441)
  
  * 2019 [Guide to presenting clinical prediction models for use in clinical settings](https://www.bmj.com/content/365/bmj.l737), [R code](https://cran.r-project.org/web/packages/rms/index.html)
  
  * 2019 [When and how to use data from randomised trials to develop or validate prognostic models](https://www.bmj.com/content/365/bmj.l2154)
  
  * 2016 [Adding tests to risk based guidelines: evaluating improvements in prediction for an intermediate risk group](https://www.bmj.com/content/354/bmj.i4450)
  
  * 2015 [How to develop a more accurate risk prediction model when there are few events](https://www.bmj.com/content/351/bmj.h3868)
  
  * 2012 [Comparing risk prediction models](https://www.bmj.com/content/344/bmj.e3186)
  
  * 2012 *Heart* [Part 1. Risk prediction models: I. Development, internal validation, and assessing the incremental value of a new biomarker](https://heart.bmj.com/content/98/9/683)
  
  * 2012 *Heart* [Part 2. Risk prediction models: II. External validation, model updating, and impact assessment](https://heart.bmj.com/content/98/9/691)


## 2. Methodological papers

* 2022 *Int J Epidemiol* [Lessons learnt when accounting for competing events in the external validation of time-to-event prognostic models](https://academic.oup.com/ije/article/51/2/615/6468864?login=false), [R code](https://github.com/survival-lumc/ValidationCompRisks)

* 2022 *Ann Intern Med* [Assessing Performance and Clinical Usefulness in Prediction Models With Survival Outcomes: Practical Guidance for Cox Proportional Hazards Models](https://www.acpjournals.org/doi/10.7326/M22-0844), [R code](https://github.com/danielegiardiello/Prediction_performance_survival)

* 2016 *J Clin Oncol* [Assessing the clinical impact of risk prediction models with decision curves: guidance for correct interpretation and appropriate Use](https://ascopubs.org/doi/10.1200/JCO.2015.65.5654), [R code](https://cran.r-project.org/web/packages/dcurves/)

* 2014 *Ann Intern Med* [Net reclassification improvement: compuataion, interpretation, and controversies](https://www.acpjournals.org/doi/10.7326/M13-1522)

* 2010 *Epidemiology* [Assessing the performance of prediction models: A framework for traditional and novel measures](https://journals.lww.com/epidem/fulltext/2010/01000/assessing_the_performance_of_prediction_models__a.22.aspx), [R code 1](https://links.lww.com/EDE/A355), [R code 2](https://github.com/danielegiardiello/ValLogRegMod)

* 2003 *J Clin Epidemiol* [Internal and external validation of predictive models: a simulation study of bias and precision in small samples](https://linkinghub.elsevier.com/retrieve/pii/S0895435603000477)

* 2001 *J Clin Epidemiol* [Internal validation of predictive models: efficiency of some procedures for logistic regression analysis](https://linkinghub.elsevier.com/retrieve/pii/S0895435601003419)


## 3. Advanced topics

### 3.1 Predictors selections

* 2023 *J Clin Epidemiol* [Practical guide to the typical analysis of prognostic factors and biomarkers without the use of P-values](https://linkinghub.elsevier.com/retrieve/pii/S0895435623000768)

* 2021 *J Stat Softw* [glmulti: An R package for easy automated model selection with (Generalized) linear models](https://www.jstatsoft.org/article/view/v034i12), [R code](https://cran.r-project.org/web/packages/glmulti/index.html)

### 3.2 Super learner

* 2018 *Epidemiology* [Using super learner prediction modeling to improve high-dimensional propensity score estimation](https://journals.lww.com/epidem/fulltext/2018/01000/using_super_learner_prediction_modeling_to_improve.13.aspx), [R code 1](https://github.com/lendle/hdps); [R code 2](https://github.com/lendle/TargetedLearning.jl)

### 3.3 Dynamic prediction

* 2021 *Am J Epidemiol* [Dynamical modeling as a tool for inferring causation](https://academic.oup.com/aje/article/191/1/1/6358337?login=false)

* 2021 *R Journal* [JMcmprsk: An R Package for Joint Modelling of Longitudinal and Survival Data with Competing Risks](https://journal.r-project.org/archive/2021/RJ-2021-028/RJ-2021-028.pdf), [R code](https://github.com/cran/JMcmprsk)

* 2020 *Int J Epidemiol* [Reflections on modern methods: Dynamic prediction using joint models of longitudinal and time-to-event data](https://academic.oup.com/ije/article/50/5/1731/6174516?login=false), [R code](https://erandrinopoulou.github.io/EducationalCorner_JMpred/index.html)

* 2015 *Int J Epidemiol* [Joint modelling of repeated measurement and time-to-event data: an introductory tutorial](https://academic.oup.com/ije/article/44/1/334/657852?login=false), [R code](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/ije/44/1/10.1093_ije_dyu262/4/dyu262_Supplementary_Data.zip?Expires=1713186801&Signature=R3cLsdcvgZ0cJAxnpDi8ZYoKkAKxUGEF4JSAvtG3GmdKQ4nn0LHv30BIkmQHYJUCnyrgNUyx2nUthGP8W4hcsa9Xi6Hr9T5I1Lq6HbaIUkFQlIcjKU0Gl1asnr7bUwF7M23sdUjvymzsc1HAgTqUgMSgk1g2CrEWp8j1uVrexd3StsaRR6hBPSONYO72mZO3IVfpQ7en0fVo8nUVwatMku28rZaHFm7dmhV3gi6CrblZ0jQf~i6kLrMmVGYZBKSGNyQ7E4~Bw~zyzxHCV-rOpOoiAdLV8q7Bhv6p20-CLcLHP9IZUlHNmDlF-6SN0nmAOOjWNTl8HjboweaovJY7lQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

* 2013 *Biometrics* [Dynamic pseudo-observations: a robust approach to dynamic prediction in competing risks](https://academic.oup.com/biometrics/article/69/4/1043/7492355?login=false), [R code](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/biometrics/69/4/10.1111_biom.12061/5/biom12061-sm-0001-supinfo-s1.pdf?Expires=1713186954&Signature=1edWKfgynE-cL-7qyr7fHvS8izB29WPN89pDi4qSk5548tyLZbFfxUnoT1BwDuB~RqS8ukHv2-6lw2FqZguU60uMrY8D-lZ02xalQWjPd0LmYlw1SUblFbetY58v~6tPE2~pcFspHXUnGaN8rZVGuiSpv4Acb4Vkpp4bPtqNh8Xicu376oiXCV9YRWmWVggda3jx5gQ9Sp2rG3OTCkjXsVrL1jxVltdqeyBOBoUpuStSpFK-r-~4c5ChijDim49y4LD2OxNYavA0C0hDeqdSt9hnLJDiir0AzZKNF4yU9YrMC5TMGyL0PvQaE3oOm7rTIPqu14FkSjWZBNlZUZHTcA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

### 3.4 Prediction of counterfactual risk

* 2022 *Eur J Epidemiol* [Predicting counterfactual risks under hypothetical treatment strategies: an application to HIV](https://link.springer.com/article/10.1007/s10654-022-00855-8)

* 2022 *Eur J Epidemiol* [A structural characterization of shortcut features for prediction](https://link.springer.com/article/10.1007/s10654-022-00892-3)

* 2020 *Eur J Epidemiol* [Counterfactual prediction is not only for causal inference](https://link.springer.com/article/10.1007/s10654-020-00659-8)

* 2020 *Eur J Epidemiol* [Prediction meets causal inference: the role of treatment in clinical prediction models](https://link.springer.com/article/10.1007/s10654-020-00636-1)

### 3.5 Lifetime risk

* 2021 *Eur J Epidemiol* [A comparison of statistical methods to predict the residual lifetime risk](https://link.springer.com/article/10.1007/s10654-021-00815-8), [R code](https://github.com/s-conner/lifetimerisk)

* 2000 *Stat Med* [Computing estimates of incidence, including lifetime risk: Alzheimer’s disease in the framingham study. the practical incidence estimators (pie) macro](https://onlinelibrary.wiley.com/doi/10.1002/(SICI)1097-0258(20000615/30)19:11/12%3C1495::AID-SIM441%3E3.0.CO;2-E)


### 3.6 Multistate model

* 2016 *Circ Heart Fail* [Multistate model to predict heart failure hospitalizations and all-cause mortality in outpatients with heart failure with reduced ejection fraction: Model derivation and external validation](https://www.ahajournals.org/doi/10.1161/CIRCHEARTFAILURE.116.003146)

### 3.7 Deep learning 

* 2024 *arXiv* [Deep learning for survival analysis: A review](https://arxiv.org/abs/2305.14961), [R code](https://survival-org.github.io/DL4Survival/)

* 2017 *arXiv* [Machine Learning for Survival Analysis: A Survey](https://arxiv.org/abs/1708.04649), [R code](https://github.com/MLSurvival)

## 4 Well-established risk prediction models

### 4.1 Framingham Heart Study risk scores

* 1988 *Circulation* Original FHS CHD score [Prediction of coronary heart disease using risk factor categories](https://www.ahajournals.org/doi/10.1161/01.CIR.97.18.1837)

* 1999 *Lancet* [Lifetime risk of developing coronary heart disease](https://linkinghub.elsevier.com/retrieve/pii/S0140673698102799)  

* 2006 *Circulation* USA-PRC [Estimation of 10-year risk of fatal and nonfatal ischemic cardiovascular diseases in Chinese adults](https://pubmed.ncbi.nlm.nih.gov/17088464/)

* 2009 *Circulation* FHS 30-year CVD risk prediction model [Predicting the 30-year risk of cardiovascular disease: The Framingham Heart Study](https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.108.816694)


### 4.2 QRISK

* 2007 *BMJ* QRISK [Derivation and validation of QRISK, a new cardiovascular disease risk score for the United Kingdom: prospective open cohort study](https://www.bmj.com/content/335/7611/136)

* 2008 *BMJ* QRISK2 [Predicting cardiovascular risk in England and Wales: prospective derivation and validation of QRISK2](https://www.bmj.com/content/336/7659/1475)

* 2010 *BMJ* QRISK-lifetime risk [Derivation, validation, and evaluation of a new QRISK model to estimate lifetime risk of cardiovascular disease: cohort study using QResearch database](https://www.bmj.com/content/341/bmj.c6624)

* 2017 *BMJ* QRISK3 [Development and validation of QRISK3 risk prediction algorithms to estimate future risk of cardiovascular disease: prospective cohort study](https://www.bmj.com/content/357/bmj.j2099)

* 2024 *Eur J Prev Cardiol* QRISK2-PRS [A polygenic risk score added to a QRISK®2 cardiovascular disease risk calculator demonstrated robust clinical acceptance and clinical utility in the primary care setting](https://academic.oup.com/eurjpc/advance-article/doi/10.1093/eurjpc/zwae004/7577877)


### 4.3 CMCS Risk prediction models

* 2003 *中华心血管病杂志* [中国35~64岁人群心血管病危险因素与发病危险因素预测模型的前瞻性研究](https://rs.yiigle.com/CN115399202004/978228.htm)

* 2004 *JAMA* FHS CHD score in Chinese population [Predictive Value for the Chinese Population of the Framingham CHD Risk Assessment Tool Compared With the Chinese Multi-provincial Cohort Study](https://jamanetwork.com/journals/jama/fullarticle/198843)

* 2013 *Eur J Prev Cardiol* CMCS CVD lifetime risk [Lifetime risk for cardiovascular disease in a Chinese population: the Chinese Multi-Provincial Cohort Study](https://academic.oup.com/eurjpc/article/22/3/380/5984131?login=false)

* 2016 *J Hypertens* [Lifetime risk of stroke in young-aged and middle-aged Chinese population: the Chinese Multi-Provincial Cohort Study](https://journals.lww.com/jhypertension/fulltext/2016/12000/lifetime_risk_of_stroke_in_young_aged_and.19.aspx)

* 2018 *中华心血管病杂志* [中国动脉粥样硬化性心血管病发病危险评估的新方案](https://rs.yiigle.com/CN2021/1027278.htm)

* 2023 *CVIA* [Addition of Risk-enhancing Factors Improves Risk Assessment of Atherosclerotic Cardiovascular Disease in Middle-aged and Older Chinese Adults: Findings from the Chinese Multi-provincial Cohort Study](https://www.scienceopen.com/hosted-document?doi=10.15212/CVIA.2023.0036)


### 4.4 Pooled Cohort ASCVD Risk Equations [ASCVD Risk Calculator](https://clincalc.com/cardiology/ascvd/pooledcohort.aspx)

* 2013 *Circulation* PCE [2013 ACC/AHA guideline on the assessment of cardiovascular risk: a report of the American College of Cardiology/American Heart Association Task Force on Practice Guidelines](https://www.ahajournals.org/doi/full/10.1161/01.cir.0000437741.48606.98)

* 2013 *J Am Coll Cardiol* PCE [2013 ACC/AHA guideline on the treatment of blood cholesterol to reduce atherosclerotic cardiovascular risk in adults: a report of the American College of Cardiology/American Heart Association Task Force on Practice Guidelines](https://www.sciencedirect.com/science/article/pii/S0735109713060282?via%3Dihub)

* 2017 *Circulation* PCE lifetime risk [Estimating Longitudinal Risks and Benefits From Cardiovascular Preventive Therapies Among Medicare Patients: The Million Hearts Longitudinal ASCVD Risk Assessment Tool: A Special Report From the American Heart Association and American College of Cardiology](https://www.ahajournals.org/doi/10.1161/CIR.0000000000000467)

* 2020 *JAMA Netw Open* PCE by BMI [Performance of the Pooled Cohort Equations to Estimate Atherosclerotic Cardiovascular Disease Risk by Body Mass Index](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2772343)

* 2021 *JAMA Cardiol* PCE by self-reported physical activity [Performance of the American Heart Association/American College of Cardiology Pooled Cohort Equations to Estimate Atherosclerotic Cardiovascular Disease Risk by Self-reported Physical Activity Levels](https://jamanetwork.com/journals/jamacardiology/fullarticle/2779383) 


### 4.5 China-PAR, [心脑血管风险评估网站](https://cvdrisk.com.cn/ASCVD/Eval)

* 2016 *Circulation* China-PAR ASCVD risk [Predicting the 10-Year Risks of Atherosclerotic Cardiovascular Disease in Chinese Population: The China-PAR Project (Prediction for ASCVD Risk in China)](https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.116.022367)

* 2016 *Chronic Dis Transl Med* China-PAR ASCVD risk stratification [Risk stratification of atherosclerotic cardiovascular disease in Chinese adults](https://onlinelibrary.wiley.com/doi/10.1016/j.cdtm.2016.10.001)

* 2018 *Sci Bull (Beijing)* China-PAR ASCVD lifetime risk [Predicting lifetime risk for developing atherosclerotic cardiovascular disease in Chinese population: the China-PAR project](https://linkinghub.elsevier.com/retrieve/pii/S2095927318302421)

* 2022 *Eur Heart J* China-PAR PRS [A polygenic risk score improves risk stratification of coronary artery disease: a large-scale prospective Chinese cohort study](https://academic.oup.com/eurheartj/article/43/18/1702/6534984?login=false)

* 2019 *中国循环杂志* [中国心血管病风险评估和管理指南](https://chinacirculation.org/UploadFile/SiteContent/FJList/1ncfl2hj.pdf)


### 4.6 SCORE project

* 2003 *Eur Heart J* SCORE [Estimation of ten-year risk of fatal cardiovascular disease in Europe: the SCORE project](https://academic.oup.com/eurheartj/article/24/11/987/427645?login=false)

* 2021 *Eur Heart J* [SCORE2 risk prediction algorithms: new models to estimate 10-year risk of cardiovascular disease in Europe](https://academic.oup.com/eurheartj/article/42/25/2439/6297709)

* 2021 *Eur Heart J* [SCORE2-OP risk prediction algorithms: estimating incident cardiovascular event risk in older persons in four geographical risk regions](https://academic.oup.com/eurheartj/article/42/25/2455/6297711?login=false)


### 4.7 Others 

* 2016 *Eur Heart J* CHD GRS [Genomic prediction of coronary heart disease](https://doi.org/10.1093/eurheartj/ehw450)

* 2019 *Lancet Glob Health* WHO CVD Risk Charts [World Health Organization cardiovascular disease risk charts: revised models to estimate risk in 21 global regions](https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(19)30318-3/fulltext)

* 2020 *Eur Heart J* LIFE-CVD [Prediction of individualized lifetime benefit from cholesterol lowering, blood pressure lowering, antithrombotic therapy, and smoking cessation in apparently healthy people](https://doi.org/10.1093/eurheartj/ehz239)

* 2023 *Chin Med J (Engl)* CKB PRS [Minimal improvement in coronary artery disease risk prediction in Chinese population using polygenic risk scores: evidence from the China Kadoorie Biobank](https://journals.lww.com/cmj/fulltext/2023/10200/minimal_improvement_in_coronary_artery_disease.10.aspx)

* 2024 *Eur Heart J* CVD risk communication [Cardiovascular disease risk communication and prevention: a meta-analysis](https://academic.oup.com/eurheartj/advance-article/doi/10.1093/eurheartj/ehae002/7578334?login=false)
