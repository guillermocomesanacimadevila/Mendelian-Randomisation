##Install Packages
#TwoSample MR package - documentation available at https://mrcieu.github.io/TwoSampleMR/
install.packages("remotes") #Installing a predefined package <- "remotes"
remotes::install_github("MRCIEU/TwoSampleMR@0.4.26") #Installing GitHub derived MR function package

force = TRUE #?

##load packages
library(TwoSampleMR)


##set working directory
setwd("/Users/guillermocomesanacimadevila/Desktop/(MSc) HUMAN NUTRITION/GENETIC EPIDEMIOLOGY PROJECT/R DATA")

#########################################################
#set the following parameters:
#pval<- 5e-08 
#r2<-0.001

#In this case, we don't have to do this because we are using a pre-selected list of genetic instruments for iron status. If we didn't have this, we would have
#to select genetic instruments ourselves on the basis of the P-value from the GWAS (<5E-08) and an r2 threshold of 0.001 to ensure that we are selecting independent genetic instruments
#This is usually done through "clumping" - this will be explained in recommended reading. Note we skip the clumping step in this code.

#Here we are pre-allocating the "results" data frame
results<-data.frame(
  exposure=character(),outcome=character(),snps=numeric(),
  ivw.OR=numeric(),ivw.l=numeric(),ivw.u=numeric(),ivw.p=numeric(),
  wm.OR=numeric(),wm.l=numeric(),wm.u=numeric(),wm.p=numeric(),
  egger.OR=numeric(),egger.l=numeric(),egger.u=numeric(),egger.p=numeric(),
  Egger.intercept.p=numeric(),Q.p=numeric(),I2=numeric(),F.statistic=numeric(),
  stringsAsFactors = F
)


# Load exposure dataset: serum iron levels - you can download the data here: https://www.decode.com/summarydata/ #https://www.nature.com/articles/s42003-020-01575-zhttps://www.nature.com/articles/s42003-020-01575-z
#In this case, we have already selected the genetic instruments based on published work: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10773400/
exposure <- read.csv("/Users/guillermocomesanacimadevila/Desktop/(MSc) HUMAN NUTRITION/GENETIC EPIDEMIOLOGY PROJECT/R DATA/iron-instruments.csv", header = TRUE)

  # We need to convert the loaded exposure dataframe into a format that the "TwoSampleMR" package can understand and feed in to the MR computation
  exp<-format_data(
    exposure,
    type= "exposure",
    snps = NULL,
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "p")
  
  # name exposure correctly
  exp$exposure <- "serum iron levels"
  
  
  # Load outcome dataset: Pernicious Anaemia - you can download the data here: http://www.geenivaramu.ee/tools/pernicious_anemia_Laisketal2021_sumstats.gz Paper:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8213695/ 
  outcome <- read.table("/Users/guillermocomesanacimadevila/Desktop/(MSc) HUMAN NUTRITION/GENETIC EPIDEMIOLOGY PROJECT/R DATA/b12_meta_ukbb_estbb_finngen_061020.out", header=TRUE, sep="\t")

  ##Extract Exposure SNPs from Outcome data
  outc <- outcome[outcome$rs_number %in% exp$SNP, ]
  
    ## Convert outcome to "TwoSample" Format
    outc <- format_data(
      dat = outc,
      type = "outcome",
      snps = outc$rs_number,
      snp_col = "rs_number",
      beta_col = "beta",
      se_col = "se",
      eaf_col = "eaf",
      effect_allele_col = "reference_allele",
      other_allele_col = "other_allele",
      pval_col = "p.value")
    
    
    ## name outcome
    outc$outcome <- "Pernicious Anaemia"
    
    # Harmonise data
    dat <- harmonise_data(exp, outc, action = 1)
    
    
    # Perform MR (IVW, WM, Egger)
    # F statistic using chi-squared approximation - F statistic will be explained in recommended readings
    dat$Fstat <- dat$beta.exposure^2 / dat$se.exposure^2
    
    # Perform MR
    res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
    
    # MR-Egger intercept
    pleio <- mr_pleiotropy_test(dat)
    
    # heterogeneity test
    het <- mr_heterogeneity(dat)
    het$I2 <- ifelse((het$Q - het$Q_df ) / het$Q < 0, 0, 100*(het$Q - het$Q_df ) / het$Q)
    
    # assign values in table
    results[1,1]<- "serum iron" #exposure name
    results[1,2]<-"PA Laisk et al" #outcome name
    results[1,3]<-res[1,6] #N SNPs
    
    results[1,4]<-exp(res[1,7]) #IVW OR
    results[1,5]<-exp(res[1,7]+qnorm(.025)*res[1,8])#IVW lower 95%
    results[1,6]<-exp(res[1,7]+qnorm(.975)*res[1,8])#IVW upper 95%
    results[1,7]<-res[1,9]#IVW p-value
    
    results[1,8]<-exp(res[3,7])#WM beta
    results[1,9]<-exp(res[3,7]+qnorm(.025)*res[3,8])#WM lower 95%
    results[1,10]<-exp(res[3,7]+qnorm(.975)*res[3,8])#WM upper 95%
    results[1,11]<-res[3,9]#WM p-value
    
    results[1,12]<-exp(res[2,7])#Egger beta
    results[1,13]<-exp(res[2,7]+qnorm(.025)*res[2,8])#Egger lower 95%
    results[1,14]<-exp(res[2,7]+qnorm(.975)*res[2,8])#Egger upper 95%
    results[1,15]<-res[2,9]#Egger p-value
    
    results[1,16]<-pleio[1,7]#Egger intercept p
    results[1,17]<-het[2,8]#IVW Q statistic p
    results[1,18]<-het[2,9]# I2 statistic
    results[1,19]<-mean(dat$Fstat)#mean F statistic
    
    # write results: adjust your directory accordingly
    write.table(results, file = "/Users/guillermocomesanacimadevila/Desktop/(MSc) HUMAN NUTRITION/GENETIC EPIDEMIOLOGY PROJECT/R DATA/MR-results-iron-PA-Laisketal.tsv",sep = "\t", col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE)
    
