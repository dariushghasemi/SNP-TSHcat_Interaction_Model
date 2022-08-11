# SNP-TSHcat_Interaction_Model


#### 1. We first compute the eGFR using CKD-Epi formula in Nephro package in R using the Standardized Serum Creatinine.
```Rscript
### Computing eGFR.SCr
library(nephro)
table(chris$Sex) #check if SEX needs to be recoded

chris$ethnicity <- 0 # 0 for Europeans, and 1 for Africans
chris$Sex      <- as.numeric(as.character(chris$Sex))
chris$eGFR     <- CKDEpi.creat(chris$SerumCreatinine.Std, chris$Sex, chris$Age, chris$ethnicity)
chris$eGFR.log <- log(chris$eGFR)
```Rscript

#### Winsorizing lower tail of eGFR distribution
```Rscript
chris$eGFRw      <- chris$eGFR
chris[chris$eGFR < 15 & is.na(chris$eGFR) != TRUE, "eGFRw"] <- 15
chris$eGFRw.log  <- log(chris$eGFRw)
```
#### 2.Then we adjust eGFR for the Age and Sex of the CHRIS participants
```Rscript
# eGFR.log residuals by new scheme
chris$Sex           <- as.factor(chris$Sex)
meGFRw.log          <- lm(eGFRw.log ~ Age + Sex, data = chris)
myresidmeGFRw.log   <- data.frame(resm = meGFRw.log$residuals)
chris$eGFRw.log.Res <- NA
chris$eGFRw.log.Res[as.numeric(rownames(myresidmeGFRw.log))] <- myresidmeGFRw.log$resm
```


#### 3. Using GWAS, we replicated 162 SNPs in 10 Loci out of the 147 associated Loci with kidney function discovered by Mattias Wuttke Meta-GWAS paper in 2019. 
 


#### 4. For CHRIS participants, we derived dosage level for those 162 replicated SNPs from TOPMed imputed VCF files and merged it with our phenotype using unique AID number.



#### 5. We then standardized TSH using quantile transformation function in Caret R package. Then we categorized TSH based on the imperical cutpoints as follows:

We've used these intervals for categorizing TSH:
| Status | Lower bound, Upper bound |
| ------ | ------------------------ |
| Hyperthyrodism |  (0, 0.400)      |
| Normal         |  (0.401, 3.799)  |
| Hypothyrodism  |  (3.800, Inf)    |

```Rscript
library(caret)
---
nq <- normalize2Reference(chris[chris$TSH.Ins == 0, "TSH"], 
                          refData = quantile(chris[chris$TSH.Ins == 1, "TSH"], 
                                             probs = seq(0,1, length.out=length(chris[chris$TSH.Ins == 0, "TSH"])), 
                                             na.rm = TRUE, names = TRUE, type = 7, digits = 7), ties = TRUE)
---
#TSH quantile normalized
chris$TSH.q <- chris$TSH
chris[chris$TSH.Ins == 0, "TSH.q"] <- nq

summary(cbind(nq, chris[chris$TSH.Ins == 0, "TSH"]))
---
#categorized TSH
chris <- chris %>% mutate(TSH_cat = cut(TSH.q, breaks = c(-Inf, 0.401, 3.799, Inf), labels = c("0", "1", "2")))
chris$TSH_cat <- as.character(chris$TSH_cat)
---
#changing the reference level to normal TSH
chris$TSH_cat <- relevel(chris$TSH_cat, ref = 2)
```
#### 6. Re-assignment of the individuals who have taken any mediacation/treatment for their thyroid dysfunction: 
We assigned these people to their actual TSH level which we expect they belong before taking the treatment and regulaing their TSH.

#### 7. And Finally here is the summary of Regression Model:
```Rscript
summary(lm(eGFRw.log.Res ~ `chr1:10599281`:TSH_cat + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = vcfReg_TSHmod))
---
Call:
lm(formula = eGFRw.log.Res ~ `chr1:10599281`:TSH_cat + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = vcfReg_TSHmod)

Residuals:
     Min       1Q   Median       3Q      Max 
-1.75039 -0.07521  0.01555  0.10180  0.48003 

Coefficients:
                           Estimate Std. Error t value    Pr(>|t|)   
(Intercept)              -0.0003528  0.0014046  -0.251    0.801662    
PC1                      -0.1037080  0.1572179  -0.660    0.509497    
PC2                      -0.5212966  0.1817797  -2.868    0.004143 ** 
PC3                      -0.7526514  0.1528275  -4.925    8.58e-07 ***
PC4                      -0.0193142  0.1637463  -0.118    0.906108    
PC5                      -0.2164566  0.1447958  -1.495    0.134971    
PC6                       0.3265085  0.1405957   2.322    0.020236 *  
PC7                      -0.3482386  0.1407388  -2.474    0.013364 *  
PC8                      -0.2431564  0.1557854  -1.561    0.118594    
PC9                       0.2113223  0.1436013   1.472    0.141164    
PC10                     -0.8620174  0.1392735  -6.189    6.28e-10 ***
`chr1:10599281`:TSH_cat1 -0.0322001  0.0090456  -3.560    0.000373 ***
`chr1:10599281`:TSH_cat0  0.0879747  0.0631732   1.393    0.163774    
`chr1:10599281`:TSH_cat2 -0.0241069  0.0263290  -0.916    0.359898    
                         Pr(>|t|)    
(Intercept)              0.801662    
PC1                      0.509497    
PC2                      0.004143 ** 
PC3                      8.58e-07 ***
PC4                      0.906108    
PC5                      0.134971    
PC6                      0.020236 *  
PC7                      0.013364 *  
PC8                      0.118594    
PC9                      0.141164    
PC10                     6.28e-10 ***
`chr1:10599281`:TSH_cat1 0.000373 ***
`chr1:10599281`:TSH_cat0 0.163774    
`chr1:10599281`:TSH_cat2 0.359898    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
---
Residual standard error: 0.1345 on 9705 degrees of freedom
  (6 observations deleted due to missingness)
Multiple R-squared:  0.01141,	Adjusted R-squared:  0.01008 
F-statistic: 8.613 on 13 and 9705 DF,  p-value: < 2.2e-16
```
