# JMirror_Experiment

This project incorporates the codes for the simulations and real data analysis of "Joint Mirror Procedure: Controlling False Discovery Rate for Identifying Simultaneous Signals".

## Simulation

* To reproduce the simulation results in Section 4, input the following codes in the terminal:

```bash
bash sh/Simu1_SepY_Dep_pi_POW_VaryAlterProp.sh &
bash sh/Simu1_SepY_Dep_pi_POW.sh &
bash sh/Simu2_RepK_001.sh &
bash sh/Simu2_RepK_003.sh &
```

* To create Figures 3 to 4 in Section 4 of the main paper, run `PlotFiles/Plot_Res_POW_VaryAlterProp.R` and `PlotFiles/Plot_Res_POW_RepK.R` (Please also modify the root). 
* To create Tables S2-S4 in Section S.IV.2 of the supplementary materials, run `PlotFiles/Plot_Res_POW.R`. 

## Real Data

* To reproduce the results in Section 5.1, run `RealData/NASAnalysis.R` (See the instruction in it, one need produce the result in [[1]](#1) firstly). 
* To reproduce the results in Section 5.2, run `RealData/CDAnalysis.R`. 


## Other Figures

* To create Figures 1-2 in the main paper, run `Illustration/IlluTwo.R` and `Illustration/IlluThree.R`.  
* To create Figures S1 in the supplementary materials, run  `Illustration/IlluFour.R`.  


## References
<a id="1">[1]</a> 
Liu, Z., Shen, J., Barfield, R., Schwartz, J., Baccarelli, A. A., & Lin, X. (2021). Large-scale hypothesis testing for causal mediation effects with applications in genome-wide epigenetic studies. Journal of the American Statistical Association, 1-15.

<a id="2">[2]</a> 
Dai, J. Y., Stanford, J. L., & LeBlanc, M. (2020). A multiple-testing procedure for high-dimensional mediation hypotheses. Journal of the American Statistical Association, 1-16. 

<a id="3">[3]</a> 
Zhao, S. D., & Nguyen, Y. T. (2020). Nonparametric false discovery rate control for identifying simultaneous signals. Electronic Journal of Statistics, 14(1), 110-142, 133. 

<a id="4">[4]</a> 
Wang, J., Gui, L., Su, W. J., Sabatti, C., & Owen, A. B. (2022). Detecting multiple replicating signals using adaptive filtering procedures. The Annals of Statistics, 50(4), 1890--1909. 

<a id="5">[5]</a> 
Dickhaus, T., Heller, R., & Hoang, A.-T. (2021). Multiple testing of partial conjunction null hypotheses with conditional p-values based on combination test statistics. arXiv preprint arXiv:2110.06692. 

<a id="6">[6]</a> 
Benjamini, Y., & Heller, R. (2008). Screening for Partial Conjunction Hypotheses. Biometrics, 64(4), 1215-1222. 
