# Overview
---
## Analysis pipeline to test and replicate the interaction of the kidney-associated SNPs with thyroid function traits in SHIP study
As we have seen a considerable shift in the distribution of the P-values of the interaction coefficients in the below model, we would like to ask your team at University of Greifswald to kindly replicate the same anslyis in your cohort study. 

##### Table of Contents  
1. [Conceptual model](#Conceptual Model)  
2. [Steps](#emphasis)  
...snip...    
<a name="conceptual model"/>
# Conceptual model

This is the conceptual model for the replication analysis:
Model: log(eGFR) ~ SNP + TSH_cat + SNP*TSH_cat + PC1 + ... + PC10 

Please note that the natural logarithm of the estimated Glomerular Filteration Rate, the depend variable in the model, log(eGFR), is the residuals of the model which indirictly adjusted eGFR for the both covariates, AGE and SEX of the individuals. 

<a name="steps"/>
# Steps
Therefore, we need several steps prior performing the interaction tests to make the outcome variable for the above model:
1. Compute eGFR using CKD-Epi formula in Nephro package in R.
2. Take natural log of the eGFR.
3. Winsorize the eGFR to calibrate into this interval: (15, 200)
4. Adjust eGFR for the Age and Sex using a multiple linear regression model.
5. Taking the residuals of the above regression model.

## SNP-TSHcat_Interaction_Model


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

#### Winsorizing lower tail and upper tail of eGFR distribution
```Rscript
chris$eGFRw      <- chris$eGFR
chris[chris$eGFR < 15  & is.na(chris$eGFR) != TRUE, "eGFRw"] <- 15
chris[chris$eGFR > 200 & is.na(chris$eGFR) != TRUE, "eGFRw"] <- 200
chris$eGFRw.log  <- log(chris$eGFRw)
```
#### 2.Then we adjust eGFR for the Age and Sex of the CHRIS participants
```Rscript
#eGFR.log residuals by the indirect scheme
chris$Sex           <- as.factor(chris$Sex)
meGFRw.log          <- lm(eGFRw.log ~ Age + Sex, data = chris)
#taking the residuals
myresidmeGFRw.log   <- data.frame(resm = meGFRw.log$residuals)
chris$eGFRw.log.Res <- NA
chris$eGFRw.log.Res[as.numeric(rownames(myresidmeGFRw.log))] <- myresidmeGFRw.log$resm
```


#### 3. Using GWAS, we replicated 162 SNPs in 10 Loci out of the 147 associated Loci with kidney function discovered by Mattias Wuttke Meta-GWAS paper in 2019. 
 


#### 4. For CHRIS participants, we derived dosage level for those 162 replicated SNPs from TOPMed imputed VCF files and merged it with our phenotype using unique AID number.



#### 5. We then categorized TSH based on the imperical cutpoints as follows:

We've used these intervals for categorizing TSH:
| Status | Lower bound, Upper bound |
| ------ | :----------------------: |
| Hyperthyrodism |  (0, 0.400]      |
| Normal         |  (0.401, 3.799]  |
| Hypothyrodism  |  (3.800, Inf)    |

```Rscript

#categorized TSH
chris <- chris %>% mutate(TSH_cat = cut(TSH.q, breaks = c(-Inf, 0.401, 3.799, Inf), labels = c("0", "1", "2")))
chris$TSH_cat <- as.character(chris$TSH_cat)
---
#changing the reference level to normal TSH
chris$TSH_cat <- relevel(chris$TSH_cat, ref = 2)
```
#### 6. Re-assignment of the individuals who have taken any mediacation/treatment for their thyroid dysfunction: 
We assigned these people to their actual TSH level which we expected they belong before having taken the treatment and regulaing their TSH like the following table:

<img width="539" alt="Screenshot 2022-08-15 010539" src="https://user-images.githubusercontent.com/47204821/184558182-9d0df21f-1f9d-4660-85f1-e2900f30c247.png">


#### 7. And Finally here is the summary of Regression Model:
```Rscript
summary(lm(eGFRw.log.Res ~ `chr1:10599281` * TSH_cat + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = vcfReg_TSHmod))
---
Call:
lm(formula = eGFRw.log.Res ~ `chr1:10599281` * TSH_cat + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = vcfReg_TSHmod)
---
Residuals:
     Min       1Q   Median       3Q      Max 
-1.75128 -0.07534  0.01664  0.10113  0.47925 
---
Coefficients:
                           Estimate Std. Error t value Pr(>|t|)    
(Intercept)               0.0004206  0.0014824   0.284 0.776611    
`chr1:10599281`          -0.0331346  0.0090480  -3.662 0.000251 ***
TSH_cat0                  0.0486289  0.0109097   4.457 8.39e-06 ***
TSH_cat2                 -0.0172469  0.0048188  -3.579 0.000346 ***
PC1                      -0.1036469  0.1569589  -0.660 0.509048    
PC2                      -0.5159122  0.1814875  -2.843 0.004483 ** 
PC3                      -0.7451690  0.1525813  -4.884 1.06e-06 ***
PC4                      -0.0192509  0.1634859  -0.118 0.906266    
PC5                      -0.2101364  0.1445884  -1.453 0.146161    
PC6                       0.3277149  0.1403827   2.334 0.019593 *  
PC7                      -0.3540776  0.1405188  -2.520 0.011758 *  
PC8                      -0.2429022  0.1555288  -1.562 0.118372    
PC9                       0.2156111  0.1433751   1.504 0.132659    
PC10                     -0.8510523  0.1390583  -6.120 9.71e-10 ***
`chr1:10599281`:TSH_cat0  0.0619086  0.0649935   0.953 0.340850    
`chr1:10599281`:TSH_cat2  0.0283585  0.0282496   1.004 0.315474    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
---
Residual standard error: 0.1343 on 9703 degrees of freedom
  (6 observations deleted due to missingness)
Multiple R-squared:  0.01486,	Adjusted R-squared:  0.01334 
F-statistic:  9.76 on 15 and 9703 DF,  p-value: < 2.2e-16
```
