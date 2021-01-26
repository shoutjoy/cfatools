# CFA result output function
library(qgraph)
library(dplyr)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

library(processR)
library(predict3d)
library(mediation)
library(bda)

library(knitr)
library(stargazer)
library(lavaan)
library(semTools)
library(tibble)
library(semPlot)

library(forcats) #빈도계산 함수
library(FactoMineR)
library(factoextra)
library(HH)
library(emmeans)
library(tidyr)






cfa1 <- function(x){

  library(dplyr)
  library(knitr)
  library(stargazer)
  library(lavaan)
  library(semTools)
  library(tibble)
  library(semPlot)
  library(semTools)


options(scipen = 10)

fit.indices=c("chisq","pvalue", "df","rmsea","gfi","agfi","srmr","cfi","tli","nfi","aic","bic")
fitMeasures <- round(fitMeasures(x,fit.indices),2)
fitMeasures_s <- round(fitMeasures(x,fit.indices),3)
fitMeasures_s1<- round(fitMeasures(x,c("chisq","df","pvalue","rmsea","rmsea.ci.lower","rmsea.ci.upper","rmsea.pvalue", "cfi","tli","srmr")),3) %>% kable (format="pandoc",caption = "Model fit information")

fitMeasures <- as.data.frame(fitMeasures(x,fit.indices))
fitMeasures$chiq_df <- c("","",round(fitMeasures[1,1]/fitMeasures[3,1],2),"","","","","","","","","")
fitMeasures$critera <- c("","p.value >= 0.05(fair)","chisq/df <= 3(<5(ok)",
                         "RMSEA<0.08(fair),<0.1(ok)","GFI >= 0.9","AGFI >= 0.9","SRMR< 0.08(good),<0.1(ok)",
                         "CFI >= 0.9","TLI >= 0.9","NFI >= 0.9",
                         " lower","lower")
fitMeasures$Ref <-c("-","-","Wheaton et al.(1977)",
                    "Browne & Cudeek(1993)","Joreskog & Sorbom(1984)","Tanaka & Huba(1985)","HU & Bentler(1999)",
                    "Bentler(1990)","Bentler & Bonett(1980)","Bollen(1989)",
                    "Akaike(1973)","-")
fitMeasures$model.fit.index  <- c("absolute fit","","absolute fit",
                                  "absolute fit ","absolute fit ","absolute fit ","absolute fit ",
                                  "incremental fit","incremental fit","incremental fit",
                                  "parsimonious fit","parsimonious fit")

  #fitMeasures %>%t() %>%kable(digits=3, format="pandoc", caption="fitmeasure")

  #factor loading

  options(knitr.kable.NA="")

  factorloading <- parameterEstimates(x, standardized=TRUE) %>%
    filter(op=="=~") %>%
    mutate(stars=ifelse(pvalue < 0.001, "***",
                        ifelse(pvalue < 0.01, "**",
                               ifelse(pvalue < 0.05, "*", "")))) %>%
    mutate(label=ifelse(std.all>0.7,"Yes(Good)",ifelse(std.all>0.5,"Yes(fair)","No"))) %>%
    select("Latent"=lhs, Indicator=rhs, Estimate=est, S.E.=se,
           c.r=z, Sig.=stars, "p-value"=pvalue, std.lamda=std.all, Accept=label) %>%
    kable(digits=3, format="pandoc", caption="Factor Loadings:: (1) c.r(=Estimate/S.E) p<0.05, (2) std.damda >= 0.5")


  #reliablity Construct Reliability & Average Variance Extracted
  alpha<- reliability(x) %>%
    t() %>%
    as.data.frame() %>%rownames_to_column("Construct") %>%
        select(Construct, "Composite Reliability"=omega,
           "AVE"=avevar, "Cronbach's alpha"=alpha) %>%
    kable(digits=3, format="pandoc",
          caption="Convergent Validity: Composite validity(C.R>0.7):Bollen(1980)& Raykov(2001), Average Variance Extracted (AVE >0.5):Fornell & Larcker (1981)")

  #criterion-related validity

  cri.val <- inspect(x,"std")$psi #%>%
  #  kable(digits=3, format="pandoc", caption="criterion-related validity")

  #Discriminant Validity"
  validity <- lavInspect(x, what="cor.lv") %>%
    as.data.frame() %>%
    rownames_to_column("Construct") %>%
    cbind(sqrt.AVE=
            sqrt(reliability(x)["avevar", ])) %>%
    kable(digits=3, format="pandoc", caption="Discriminant Validity: Square Root of(AVE) > rho :: Fornell & Lacker(1981)")

  all.reuslt <-list(fitMeasures,fitMeasures_s,fitMeasures_s1, factorloading, alpha, cri.val, validity)
  all.reuslt
}


cfa2 <- function(x){

  library(dplyr)
  library(knitr)
  library(stargazer)
  library(lavaan)
  library(semTools)
  library(tibble)
  library(semPlot)


#01 fit table
  options(scipen = 100)

  fit.indices=c("chisq","pvalue", "df","rmsea","gfi","agfi","srmr","cfi","tli","nfi","aic","bic")
  fitMeasures <- round(fitMeasures(x,fit.indices),2)
  fitMeasures_s <- round(fitMeasures(x,fit.indices),2)


  fitMeasures <- as.data.frame(fitMeasures(x,fit.indices))
  fitMeasures$chiq_df <- c("","",round(fitMeasures[1,1]/fitMeasures[3,1],2),"","","","","","","","","")
  fitMeasures$critera <- c("","p.value >= 0.05(fair)","chisq/df <= 3(<5(ok)",
                           "RMSEA<0.08(fair),<0.1(ok)","GFI >= 0.9","AGFI >= 0.9","SRMR< 0.08(good),<0.1(ok)",
                           "CFI >= 0.9","TLI >= 0.9","NFI >= 0.9",
                           " lower","lower")
  fitMeasures$Ref <-c("-","-","Wheaton et al.(1977)",
                      "Browne & Cudeek(1993)","Joreskog & Sorbom(1984)","Tanaka & Huba(1985)","HU & Bentler(1999)",
                      "Bentler(1990)","Bentler & Bonett(1980)","Bollen(1989)",
                      "Akaike(1973)","-")
  fitMeasures$model.fit.index  <- c("absolute fit","","absolute fit",
                                    "absolute fit ","absolute fit ","absolute fit ","absolute fit ",
                                    "incremental fit","incremental fit","incremental fit",
                                    "parsimonious fit","parsimonious fit")
  fit <- fitMeasures %>% kable(digits=3, format="pandoc", caption="FitMeasure and criterian")

  #fitMeasures %>%t() %>%kable(digits=3, format="pandoc", caption="fitmeasure")

  #modelfit
  fitdata <- fitMeasures(x,c("chisq","df","pvalue","rmsea","rmsea.ci.lower","rmsea.ci.upper","rmsea.pvalue", "cfi","tli","srmr"))
  criteria_data = c("chisq","df","pvalue<0.05","rmsea<0.1","rmsea.ci.lower","rmsea.ci.upper","rmsea.pvalue<0.05", "cfi>0.9","tli>0.9","srmr<0.1")
  modelfitdata <-cbind("criterian"=criteria_data, "data"=round(fitdata,3))
  fitMeasures_s1 <- modelfitdata %>% kable (format="pandoc",caption = "Model fit information")




  #04 factor loading
  options(knitr.kable.NA="")
  factorloading <- parameterEstimates(x, standardized=TRUE) %>%
    filter(op=="=~") %>%
    mutate(stars=ifelse(pvalue < 0.001, "***",
                        ifelse(pvalue < 0.01, "**",
                               ifelse(pvalue < 0.05, "*", "")))) %>%
    mutate(label=ifelse(std.all>0.7,"Yes(Good)",ifelse(std.all>0.5,"Yes(fair)","No"))) %>%
    select("Latent"=lhs, Indicator=rhs, Estimate=est, S.E.=se,
           c.r=z, Sig.=stars, "p-value"=pvalue, std.lamda=std.all, Accept=label) %>%
    kable(digits=3, format="pandoc", caption="Convergent Validity(1)-Factor Loadings:: (1) c.r(=Estimate/S.E) p<0.05, (2) std.damda >= 0.5(Bagozzi & Yi(1988)")




  #
  alpha.1 <- reliability(x,return.total = T) %>%
    t() %>%
    as.data.frame() %>%
    select( "Cronbach's alpha"=alpha)#," average variance extracted(AVE)"=avevar)

  #05 Reprort cronbach, AVE, C.R
  FL.1 <- cbind(alpha.1)

  FL<-FL.1%>%kable(digits=3, format="pandoc",
                   caption="The coefficient alpha (Cronbach, 1951)")

  #criterion-related validity



  #reliablity Construct Reliability & Average Variance Extracted#
  #AVE cal

  l.matrix<- inspect(x,"std")$lambda
  l.matrix[l.matrix==0]<-NA
  AVE<-apply(l.matrix^2,2,mean,na.rm=T)

  #sqrt.AVE
  sqrt.AVE <- sqrt(AVE)
  #latent correlation
  rho <- inspect(x,"std")$psi
  #Composit reliability
  t.matrix<- inspect(x,"std")$theta
  t.matrix[t.matrix==0]<-NA
  t.matrix
  cr1<-apply(l.matrix,2,sum,na.rm=T)
  d.sum <-apply(1-l.matrix^2,2,sum,na.rm=T)
  C.R<- cr1^2/(cr1^2+d.sum)


  AVE_CR <- cbind( C.R,AVE,sqrt.AVE ) %>%kable(digits = 3,format = "pandoc", caption = "AVE(>0.5) & CR(>0.7): Fornell & Lacker(1981)")

  #06 Reprort
    FornellNacker <- cbind(rho, sqrt.AVE)
    validity <- FornellNacker%>%
    kable(digits=3, format="pandoc", caption="Discriminant Validity: rho < Square Root of(AVE)   By Fornell & Lacker(1981)")

    lv.cor <-inspect(x,"std")$psi#%>% kable(digits=3, format="pandoc", caption="latent correlation")


  all.reuslt <-list(fit,fitMeasures_s, fitMeasures_s1, factorloading, FL, AVE_CR, validity,lv.cor )
  all.reuslt
}


#cfa.plot###############################3333
cfa.plot <- function(x){
  semPaths(x, whatLabels = "std", nCharNodes = 8, #text number
           rotation = 2,
           style = "lisrel",residScale = 10, curvePivot=T,curve = 3,
           layout = "tree2", #shape
           layoutSplit = F,subScale = 0.9,subScale2 = 1, # manifest 1 row, 2column
           subRes = 4, #Default=4
           sizeLat =8,sizeLat2 = 6, sizeMan = 5, sizeMan2=3 ,
           color = list(lat="skyblue", man="Gold", int="gray80"),label.cex=1,
           edge.label.cex = 1,edge.color = "gray10", edge.label.position=0.65,
           edge.width= 0.8, asize=2,
           mar=c(4,15,4,15), bg="gray90",
           structural = F)
}

#cfa.plot2##########################333
cfa.plot2 <- function(x){
  semPaths(x, what= "std", nCharNodes = 8, #text number
           rotation = 2,
           style = "lisrel",residScale = 10, curvePivot=F,curve = 1.3,
           layout = "tree2", #shape
           layoutSplit = F,subScale = 0.9,subScale2 = 1, # manifest 1 row, 2column
           subRes = 4, #Default=4
           sizeLat =8,sizeLat2 = 6, sizeMan = 5, sizeMan2=3 ,
           color = list(lat="skyblue", man="Gold", int="gray80"),label.cex=2,
           edge.label.cex = 1.3,edge.color = "gray10", edge.label.position=0.55,
           edge.width= 0.8, asize=2,
           mar=c(4,15,4,15), bg="gray90",
           structural = F)
}

#model comparing ############################
fit_index <-function(x){

  a <-round(fitMeasures(x ,c("logl","unrestricted.logl","aic","bic","bic2","chisq","df","pvalue","rmsea","rmsea.ci.lower","rmsea.ci.upper","rmsea.pvalue", "cfi","tli","srmr")),3) %>% kable (format="pandoc",caption = "Model fit information")
  return(a)
}



CompareFit <- function(...) {
  library(magrittr)
  library(stargazer)
  library(knitr)
  m <- list(...)

  result <-sapply(m, fitMeasures) %>%
    set_colnames(paste0("Model-", 1:length(m))) %>%
    as.data.frame() %>%
    rownames_to_column("Fit_Measures") %>%
    slice(match(c("chisq", "df", "pvalue",
                   "rmsea","gfi","agfi", "srmr","cfi","tli","nfi","ecvi"), Fit_Measures)) %>%
    mutate(Fit_Measures=c("Chi-square", "df", "p-value",
                          "RMSEA(<0.08)","GFI(>0.9)","AGFI(>0.9)","SRMR(<0.08)", "CFI(>0.9)","TLI(NNFI)(>0.9)","NFI(>0.9)","ECVI(lower)"))%>%  kable(digits=3, format="pandoc", caption="model comparison")

  return(result)

  }





##validity#########################3
Validity.FL <-function(x){

#AVE cal
l.matrix<- inspect(x,"std")$lambda
l.matrix[l.matrix==0]<-NA
AVE<-apply(l.matrix^2,2,mean,na.rm=T)
#sqrt.AVE
sqrt.AVE <- sqrt(AVE)
#latent correlation
rho <- inspect(x,"std")$psi
#Composit reliability
t.matrix<- inspect(x,"std")$theta
t.matrix[t.matrix==0]<-NA
t.matrix
cr1<-apply(l.matrix,2,sum,na.rm=T)
d.sum <-apply(1-l.matrix^2,2,sum,na.rm=T)
CR<- cr1^2/(cr1^2+d.sum)
#Reprort
FornellNacker1 <- cbind(CR, AVE)
FornellNacker2 <- cbind(rho,sqrt.AVE)
FornellNacker3 <- cbind(rho,sqrt.AVE, AVE,CR)
result.1<- list(FornellNacker1, FornellNacker2,FornellNacker3)
result<- result.1 %>% kable(digits=3, format="pandoc", caption="Discriminant Validity: Square Root of(AVE) > rho :: Fornell & Lacker(1981)")
return(result)
}

#stargazer(type = "text", title = "Discriminant Validity: Square Root of(AVE) > rho :: Fornell & Lacker(1981)")



#AVE#####################33333
AVE<-function(x){
  l.matrix<- inspect(x,"std")$lambda
  l.matrix[l.matrix==0]<-NA
  AVE.1<-apply(l.matrix^2,2,mean,na.rm=T)
  AVE<- AVE.1 %>% kable(digits=3, format="pandoc", caption="AVE(average variance extracted)>0.5")
  return(AVE)
}


#CR#################33333
CR <- function(x){
inspect(x,"std")
  #lambda
inspect(x,"std")$lambda
l.matrix<- inspect(x,"std")$lambda
l.matrix[l.matrix==0]<-NA
#theta: error coefficient
t.matrix<- inspect(x,"std")$theta
t.matrix[t.matrix==0]<-NA
t.matrix
cr1<-apply(l.matrix,2,sum,na.rm=T)
d.sum <-apply(1-l.matrix^2,2,sum,na.rm=T)
CR.c<- cr1^2/(cr1^2+d.sum)
CR <-CR.c  %>% kable(digits=3, format="pandoc", caption="Construct Reliability > 0.7(Bagozzi and Yi 1988)")
return(CR)}


CR_AVE <- function(x){
  l.matrix<- inspect(x,"std")$lambda
  l.matrix[l.matrix==0]<-NA
  AVE<-apply(l.matrix^2,2,mean,na.rm=T)

  inspect(x,"std")
  #lambda
  inspect(x,"std")$lambda
  l.matrix<- inspect(x,"std")$lambda
  l.matrix[l.matrix==0]<-NA
  #theta: error coefficient
  t.matrix<- inspect(x,"std")$theta
  t.matrix[t.matrix==0]<-NA
  t.matrix
  cr1<-apply(l.matrix,2,sum,na.rm=T)
  d.sum <-apply(1-l.matrix^2,2,sum,na.rm=T)
  Contruct_Reliablity<- cr1^2/(cr1^2+d.sum)

  rel.1<- cbind(Contruct_Reliablity, AVE)
  rel<- rel.1%>% kable(digits=3, format="pandoc", caption="CR(>0.7) & AVE(>0.5)")
  return(rel)

}






loading <- function(x){

#factor loading

options(knitr.kable.NA="")

factorloading <- parameterEstimates(x, standardized=TRUE) %>%
  filter(op=="=~") %>%
  mutate(stars=ifelse(pvalue < 0.001, "***",
                      ifelse(pvalue < 0.01, "**",
                             ifelse(pvalue < 0.05, "*", "")))) %>%
  mutate(label=ifelse(std.all>0.7,"Yes(Good)",ifelse(std.all>0.5,"Yes(fair)","No"))) %>%
  select("Latent"=lhs, Indicator=rhs, Estimate=est, S.E.=se,std.lamda=std.all,
         c.r=z, Sig.=stars, "p-value"=pvalue) %>%
  kable(digits=3, format="pandoc", caption="Factor Loadings:: (1) c.r(=Estimate/S.E) p<0.05, (2) std.damda >= 0.5")

return(factorloading)

}


covariance <- function(x){

  #corvariance

  options(knitr.kable.NA="")

  covari<- parameterEstimates(x, standardized=TRUE) %>%
    filter(op=="~~") %>%
    mutate(stars=ifelse(pvalue < 0.001, "***",
                        ifelse(pvalue < 0.01, "**",
                               ifelse(pvalue < 0.05, "*", "")))) %>%
    mutate(label=ifelse(std.all>0.7,"Yes(Good)",ifelse(std.all>0.5,"Yes(fair)","No"))) %>%
    select("Error"=lhs, variable =rhs, Covariance=est,correlation=std.all,
           c.r=z, Sig.=stars, "p-value"=pvalue,CI.lower=ci.lower ,CI.upper=ci.upper) %>%
    kable(digits=3, format="pandoc", caption=" Covariance")

  return(  covari)

}

R2 <- function(x){
  # R-square

  srs <-  summary(x, rsquare=T, standardized=T)
  sr <- srs$PE %>% filter(op=="r2") %>% select(lhs,rhs, R2 =est) %>%
    arrange(lhs)%>% kable(digits=3, format="pandoc",caption = "R-square rhs(manifest)->lhs(Latent)")
  return(sr)

}



validity <- function(x){


  #
  alpha.1 <- reliability(x,return.total = T) %>%
    t() %>%
    as.data.frame() %>%
    select( "Cronbach's alpha"=alpha)#," average variance extracted(AVE)"=avevar)

  #05 Reprort cronbach, AVE, C.R
  FL.1 <- cbind(alpha.1)

  FL<-FL.1%>%kable(digits=3, format="pandoc",
                   caption="The coefficient alpha (Cronbach, 1951)")

  #criterion-related validity
  #reliablity Construct Reliability & Average Variance Extracted#
  #AVE cal

  l.matrix<- inspect(x,"std")$lambda
  l.matrix[l.matrix==0]<-NA
  AVE<-apply(l.matrix^2,2,mean,na.rm=T)

  #sqrt.AVE
  sqrt.AVE <- sqrt(AVE)
  #latent correlation
  rho <- inspect(x,"std")$psi
  #Composit reliability
  t.matrix<- inspect(x,"std")$theta
  t.matrix[t.matrix==0]<-NA
  t.matrix
  cr1<-apply(l.matrix,2,sum,na.rm=T)
  d.sum <-apply(1-l.matrix^2,2,sum,na.rm=T)
  C.R<- cr1^2/(cr1^2+d.sum)


  AVE_CR <- cbind(C.R,AVE,sqrt.AVE) %>%kable(digits = 3,format = "pandoc", caption = "AVE(>0.5) & CR(>0.7): Fornell & Lacker(1981)")

  #06 Reprort
  FornellNacker <- cbind(rho, sqrt.AVE)
  validity <- FornellNacker%>%
    kable(digits=3, format="pandoc", caption="Discriminant Validity: Square Root of(AVE) > rho :: Fornell & Lacker(1981)")
  aaa<- list(validity,AVE_CR)

aaa
}


#alpha####
Alpha <- function(x){

#reliablity Construct Reliability & Average Variance Extracted
alpha<- reliability(x) %>%
  t() %>%
  as.data.frame() %>%rownames_to_column("Construct") %>%
  select(Construct, "Composite Reliability"=omega,
         "AVE"=avevar, "Cronbach's alpha"=alpha) %>%
  kable(digits=3, format="pandoc",
        caption="Convergent Validity: Composite validity(C.R>0.7):Bollen(1980)& Raykov(2001), Average Variance Extracted (AVE >0.5):Fornell & Larcker (1981)")
return(alpha)
}


#model fit view####
model.fit<- function(x){

  options(scipen = 100)

  fit.indices=c("chisq","pvalue", "df","rmsea","gfi","agfi","srmr","cfi","tli","nfi","aic","bic")
  fitMeasures <- round(fitMeasures(x,fit.indices),2)

  fitMeasures_s <- round(fitMeasures(x,fit.indices),2)
  fitMeasures <- as.data.frame(fitMeasures(x,fit.indices))

  fitMeasures$chiq_df <- c("","",round(fitMeasures[1,1]/fitMeasures[3,1],2),"","","","","","","","","")
  fitMeasures$critera <- c("","p.value >= 0.05(fair)","chisq/df <= 3(<5(ok)",
                           "RMSEA<0.08(fair),<0.1(ok)","GFI >= 0.9","AGFI >= 0.9","SRMR< 0.08(good),<0.1(ok)",
                           "CFI >= 0.9","TLI >= 0.9","NFI >= 0.9",
                           " lower","lower")
  fitMeasures$Ref <-c("-","-","Wheaton et al.(1977)",
                      "Browne & Cudeek(1993)","Joreskog & Sorbom(1984)","Tanaka & Huba(1985)","HU & Bentler(1999)",
                      "Bentler(1990)","Bentler & Bonett(1980)","Bollen(1989)",
                      "Akaike(1973)","-")
  fitMeasures$model.fit.index  <- c("absolute fit","","absolute fit",
                                    "absolute fit ","absolute fit ","absolute fit ","absolute fit ",
                                    "incremental fit","incremental fit","incremental fit",
                                    "parsimonious fit","parsimonious fit")


  #modelfit ####
  fitdata <- fitMeasures(x,c("chisq","df","pvalue","rmsea","rmsea.ci.lower","rmsea.ci.upper","rmsea.pvalue", "cfi","tli","srmr"))
  criteria_data = c("chisq","df","pvalue<0.05","rmsea<0.1","rmsea.ci.lower","rmsea.ci.upper","rmsea.pvalue<0.05", "cfi>0.9","tli>0.9","srmr<0.1")
  modelfitdata <-cbind("criterian"=criteria_data, "data"=round(fitdata,3))
  fitMeasures_s1 <- modelfitdata %>% kable (format="pandoc",caption = "Model fit information")



  all.reuslt <-list(fitMeasures,fitMeasures_s,fitMeasures_s1)
  all.reuslt
}



#CFA
#source("C:/Users/shout/OneDrive/00_R/2020_04_park/00source/cfatools.R")

#useage
#   cfa(fit1)
#effect####
effect<- function(x){

library(dplyr)
library(stargazer)
standardizedsolution(x) %>%
  filter(op=="~"|op==":=") %>%
  mutate(stars=ifelse(pvalue<0.001,"***",
                      ifelse(pvalue<0.01,"**",
                             ifelse(pvalue<0.05,"*","")))) %>%
  mutate(op=ifelse(op=="~","<--",
                   ifelse(op==":=","effect",""))) %>%
  select(Dependent=lhs,Path=op, Independent=rhs, "Coefficient"= est.std, c.r=z,
         Sig.=stars, "p-value"=pvalue) %>%
  stargazer(type="text", title="Regression & Total Effect .", summary = FALSE,
            digits = 3, digits.extra = 0, rownames = FALSE)
}

#effect(fitM2)

regression<- function(x){

  library(dplyr)
  library(stargazer)
  library(knitr)
  standardizedsolution(x) %>%
    filter(op=="~") %>%
    mutate(stars=ifelse(pvalue<0.001,"***",
                        ifelse(pvalue<0.01,"**",
                               ifelse(pvalue<0.05,"*","")))) %>%
    mutate(op=ifelse(op=="~","<--",
                     ifelse(op==":=","effect",""))) %>%
    select(Dependent=lhs,Path=op, Independent=rhs, "Coefficient"= est.std, Z=z,
           Sig.=stars, "p-value"=pvalue) %>%
    stargazer(type="text", title="Latent Regression coefficients.", summary = FALSE,
             digits = 3, digits.extra = 0, rownames = FALSE) #%>%
  #kable(digits=3, format="pandoc", caption="Regression coefficients & Interaction term check")
}



regression.gruop <- function(x){

  library(dplyr)
  library(stargazer)
  library(knitr)

  standardizedsolution(x) %>%
    filter(op=="~") %>%
    mutate(stars=ifelse(pvalue<0.001,"***",
                        ifelse(pvalue<0.01,"**",
                               ifelse(pvalue<0.05,"*","")))) %>%
    mutate(op=ifelse(op=="~","<--",
                     ifelse(op==":=","effect",""))) %>%
    select(Group=group, Dependent=lhs,Path=op, Independent=rhs, "Coefficient"= est.std, Z=z,
           Sig.=stars, "p-value"=pvalue) %>%
    stargazer(type="text", title="Group comparing ", summary = FALSE,
              digits = 3, digits.extra = 0, rownames = FALSE) #%>%
  #kable(digits=3, format="pandoc", caption="Regression coefficients & Interaction term check")
}



#latent variable correlation
lvcor <- function(x){
  inspect(x, "std")$psi %>%
    kable(digits=3, format="pandoc", caption="latent variable correlation")
   # stargazer(type="text", title="latent variable correlation  ", summary = FALSE,
    #          digits = 3, digits.extra = 0, rownames = FALSE)
}

#
#var1= c("sux1","sux2","sux3","sux4","sux5","sux6")
#var2 = c("cux1","cux2","cux3","cux4","cux5")
mod.maker <- function(data,var1,var2){
  library(semTools)

tam1.mod <- indProd(data, var1 =var1 , var2 = var2, match=FALSE, meanC=TRUE, residualC=FALSE, doubleMC=TRUE)

names(tam1.mod)
}


#standard Diagram
Diagram2 <- function(x){

  semPaths(x, whatLabels = "std", nCharNodes = 4, #text number
           rotation = 2,
           style = "lisrel",residScale = 4, curvePivot=TRUE,
           layout = "tree", #shape
           layoutSplit = F,subScale = 0.9, subScale2 = 1, # manifest 1 row, 2column
           subRes = 2, #Default=4
           sizeLat =6 , sizeMan = 4, color = list(lat="skyblue", man="yellow", int="gray80"),label.cex=2,
           edge.label.cex = 1.3, edge.color = "gray10", edge.label.position=0.5, edge.width= 3, asize=1,
           mar=c(4,4,4,4), bg="gray86",
           structural = F
  )
}


#standard Diagram
Diagram1 <- function(x){

  semPaths(x, whatLabels = "est", nCharNodes = 4, #text number
           rotation = 1,
           style = "lisrel",residScale = 4, curvePivot=TRUE,
           layout = "tree", #shape
           layoutSplit = F,subScale = 0.9,subScale2 = 1, # manifest 1 row, 2column
           subRes = 2, #Default=4
           sizeLat =6 , sizeMan = 4, color = list(lat="skyblue",
                                                  man="yellow", int="gray80"),
           label.cex=2,
           edge.label.cex = 1.3, edge.color = "gray10",
           edge.label.position=0.5, edge.width= 3, asize=1,
           mar=c(4,4,4,4), bg="gray86",nDigits = 3,
           structural = F
  )
}



Diagram.t3 <- function(x,lat,man,labelcex, est_std,rotate,layout){


  semPaths(x, whatLabels = est_std, nCharNodes = 8, #text number
           rotation = rotate,
           style = "lisrel",residScale = 4, curvePivot=TRUE,
           layout = layout, #shape
           layoutSplit = F ,subScale = 0.9,subScale2 = 1, # manifest 1 row, 2column
           subRes = 4, #Default=4
           sizeLat = lat , sizeMan = man, color = list(lat="skyblue", man="orange", int="gray80"),label.cex=2,
           edge.label.cex = labelcex, edge.color = "gray10", edge.label.position=0.5, edge.width= 3, asize=1,
           mar=c(4,2,4,2), bg="gray80",
           structural = F
  )

}


Diagram.park <- function(x,lat,man,labelcex, est_std,rotate,layout,label_position,curve ,rscale, struct_TF,title_name){

    semPaths(x, whatLabels = est_std, nCharNodes = 8, #text number
           rotation = rotate, #default=1
           style = "lisrel",residScale = rscale, #Default=6
           curvePivot=TRUE,curve = curve,intercepts = F, residuals = T,
           layout = layout, #shape
           layoutSplit = F ,subScale = 0.6 ,subScale2 = 1, # manifest 1 row, 2column
           subRes = 4, #Default=4
           sizeLat = lat , sizeLat2 = 4, sizeMan = man, sizeMan2 = 3,
           color = list(lat="skyblue", man="Gold", int="gray80"),
           label.cex=2,
           edge.label.cex = labelcex, edge.color = "gray10", edge.label.position=label_position,
           edge.width= 2, asize=1.2,
           optimizeLatRes = T,
           mar=c(6,4,6,4), bg="gray80", negCol="red",cut=0,
           structural = struct_TF ,nDigits = 2
           )
  title(title_name,line=3)
}



Diagram.park2 <- function(x,lat,man,labelcex, est_std,rotate,layout,label_position,cur ,rscale, struct_TF,title_name){

  semPaths(x, what = est_std, nCharNodes = 8,fade=T, #text number
           rotation = rotate, #default=1
           style = "lisrel",residScale = rscale, #Default=6
           curvePivot=TRUE,curve = cur,intercepts = F,
           layout = layout, #shape
           layoutSplit = F ,subScale = 0.6 ,subScale2 = 1, # manifest 1 row, 2column
           subRes = 4, #Default=4
           sizeLat = lat, sizeMan = man, sizeLat2 = 4,sizeMan2 = 3,
           color = list(lat="skyblue", man="Gold", int="gray80"),
           label.cex=2,
           edge.label.cex = labelcex, edge.color = "gray20", edge.label.position=label_position,
           edge.width= 1, asize=1.2,
           optimizeLatRes = T,
           mar=c(6,4,6,4), bg="gray80",negCol="red",cut=0,
           structural = struct_TF ,nDigits = 2 )
  title(title_name,line=3)

}




Diagram.park3 <- function(x,lat,lat2, man, man2,labelcex, est_std,rotate,layout,label_position,curve, residual, rscale, struct_TF,title_name,mar1,mar2,mar3,mar4){

  semPaths(x, whatLabels = est_std, nCharNodes = 8, #text number
           rotation = rotate, #default=1
           style = "lisrel",residScale = rscale, #Default=6
           curvePivot=TRUE,curve = curve,intercepts = F, residuals = residual,
           layout = layout, #shape
           layoutSplit = F ,subScale = 0.6 ,subScale2 = 1, # manifest 1 row, 2column
           subRes = 4, #Default=4
           sizeLat = lat , sizeMan = man, sizeLat2 = lat2,sizeMan2 = man2,
           color = list(lat="skyblue", man="Gold", int="gray80"),
           label.cex=2,
           edge.label.cex = labelcex, edge.color = "gray10", edge.label.position=label_position,
           edge.width= 2, asize=1.2,
           optimizeLatRes = T,
           mar=c(mar1,mar2,mar3,mar4), bg="gray80", negCol="red",cut=0,
           structural = struct_TF ,nDigits = 2  )
  title(title_name,line=3)
}




Diagram.what <- function(x, edge_width){

  semPaths(x, what = "std", nCharNodes = 4, #text number
           rotation = 2,
           style = "lisrel",residScale = 4, curvePivot=TRUE,
           layout = "tree2", #shape
           layoutSplit = F,subScale = 0.9,subScale2 = 1, # manifest 1 row, 2column
           subRes = 2, #Default=4
           sizeLat =6 , sizeMan = 4, color = list(lat="skyblue", man="orange", int="gray80"),label.cex=2,
           edge.label.cex = 1.3, edge.color = "darkblue", edge.label.position=0.5, edge.width= edge_width, asize=2,weighted=T,
           mar=c(4,4,4,4), bg="gray81",
           structural = F
  )
}



Diagram.w2 <- function(x,l_cex,titlename){


  semPaths(x, what="std",edge.label.cex=l_cex,style="lisrel",intercepts = TRUE,fade=FALSE,
           layout = "tree",
           rotation = 1,curvePivot = TRUE, curve=0.6,optimizeLatRes = TRUE,
           residScale= 6,
           sizeLat=5, sizeMan=4,nCharNodes=8,label.cex=2,edge.label.position=0.6,
           residuals=TRUE,asize=1.5,
           color = list(lat="skyblue", man="orange", int="gray80"),
           bg="gray90" )
  title(titlename,line=3)

}




Diagramcfa <- function(x,rotate){

  semPaths(x, what= "std", nCharNodes = 8, #text number
           rotation = rotate,
           style = "lisrel",residScale = 6, curvePivot=TRUE,
           layout = "tree", #shape
           layoutSplit = F,subScale = 0.9,subScale2 = 1, # manifest 1 row, 2column
           subRes = 2, #Default=4
           sizeLat =12,sizeLat2 = 5, sizeMan = 8, sizeMan2=3 ,
           color = list(lat="skyblue", man="Gold", int="gray80"),label.cex=1.5,
           edge.label.cex = 1.1,edge.color = "gray10", edge.label.position=0.7,
           edge.width= 0.8, asize=2,
           mar=c(4,10,4,10), bg="gray80", cutoff= 0.6,
           structural = F
  )
}


Diagram.full <- function(x,lat,man,labelcex, est_std,rotate,layout, edge_label_pos,strut_T_F){

  semPaths(x, whatLabels = est_std, nCharNodes = 8, #text number
           rotation = rotate,
           style = "lisrel",residScale = 4, curvePivot=TRUE,
           layout = layout, #shape
           layoutSplit = F ,subScale = 1,subScale2 = 1, # manifest 1 row, 2column
           subRes = 2, #Default=4
           sizeLat = lat , sizeMan = man, color = list(lat="skyblue", man="yellow", int="gray80")
           ,label.cex=2,
           edge.label.cex = labelcex, edge.color = "gray10", edge.label.position= edge_label_pos,
           edge.width= 3, asize=1,
           mar=c(4,4,4,4), bg="gray84",
           structural = strut_T_F  )
}


Diagram <- function(x,rotate,lay){
  semPaths(x, what= "std",
           style = "lisrel",residScale = 7,  subScale = 0.3, curvePivot=TRUE, curve = 0.3,
           layout = lay , layoutSplit = F,
           optimizeLatRes = F,
           intercepts = F,
           nCharNodes = 8,
           rotation=rotate,
           subRes = 4,
           edge.label.cex = 1.2, edge.color = "gray10", edge.label.position=0.6,edge.width= 1, asize=1,
           sizeLat = 8 ,sizeLat2 =5 , sizeMan = 4,sizeMan2 = 2,
           color = list(lat="skyblue", man="Gold"),
           bg="gray84",negCol="red", label.cex=2, mar=c(6,6,6,6), title=T)
}




#des ---- descrive data####
des<- function(input_data){
  library(psych)
  library(knitr)

  input_data<-as.data.frame(input_data)
  output_data <- describeBy(input_data,group = NULL)
  output_data$skew.z <- output_data$skew/sqrt(6/output_data$n)
  output_data$kurt.z <- output_data$kurtosis/sqrt(24/output_data$n)
  output_data
  output_data<-as.data.frame(output_data)
  a1 <-output_data%>% kable(digits=3, format="pandoc", caption="Describe Statistics data")
  output<-round(output_data[,c("vars","n","mean","sd","skew", "skew.z", "kurtosis","kurt.z")],3)
  #data refor.
  output[,"skew_TF"]<- "Not"
  output[output$skew.z < 3,"skew_TF"]<- "fair.norm"
  output[output$skew.z < 2,"skew_TF"]<- "Good.norm"
  output[output$skew.z < -2,"skew_TF"]<- "fair.norm"
  output[output$skew.z < -3,"skew_TF"]<- "Not"

  output[,"kurt_TF"]<- "Not"
  output[output$kurt.z < 10,"kurt_TF"]<- "fair.norm"
  output[output$kurt.z < 7,"kurt_TF"]<- "Good.norm"
  output[output$kurt.z < -7,"kurt_TF"]<- "fair.norm"
  output[output$kurt.z < -10,"kurt_TF"]<- "Not"

    a2 <-output%>% kable(digits=3, format="pandoc", caption="Describe Skew.z=[Skew/sqrt(6/n)]) and Kurt.z=[kurt/sqrt(24/n)]")
    a3 <-output
    #Kline(2011)
    all<- list(a1,a2,a3)
    all
}

#des2----descrive data#####
des.1 <- function(input_data){
  library(psych)
    library(knitr)
  input_data<-as.data.frame(input_data)
  output_data <- describeBy(input_data,group = NULL)
  output_data$skew.z <- output_data$skew/sqrt(6/output_data$n)
  output_data$kurt.z <- output_data$kurtosis/sqrt(24/output_data$n)
  output_data
  output_data<-as.data.frame(output_data)
  a1 <-output_data
  #%>% kable(digits=3, format="pandoc", caption="Describe Statistics data")
  output<-round(output_data[,c("n","mean","sd","skew", "skew.z", "kurtosis","kurt.z")],3)
  #data refor.
  output[,"skew_TF"]<- "Not"
  output[output$skew.z < 3,"skew_TF"]<- "fair.norm"
  output[output$skew.z < 2,"skew_TF"]<- "Good.norm"
  output[output$skew.z < -2,"skew_TF"]<- "fair.norm"
  output[output$skew.z < -3,"skew_TF"]<- "Not"

  output[,"kurt_TF"]<- "Not"
  output[output$kurt.z < 10,"kurt_TF"]<- "fair.norm"
  output[output$kurt.z < 7,"kurt_TF"]<- "Good.norm"
  output[output$kurt.z < -7,"kurt_TF"]<- "fair.norm"
  output[output$kurt.z < -10,"kurt_TF"]<- "Not"

  a2 <-output
  #%>% kable(digits=3, format="pandoc", caption="Describe Skew.z=[Skew/sqrt(6/n)]) and Kurt.z=[kurt/sqrt(24/n)]")
  all<- list(a1,a2)
  all
}



#descrive data########
#왜도 첨도 표준화는 극한값을 이용함
des.2<- function(input_data, digit){
  library(psych)
  library(knitr)
  input_data<-as.data.frame(input_data)
  output_data <- describeBy(input_data, group = NULL)
  output_data$skew.z <- output_data$skew/sqrt(6/output_data$n)
  output_data$kurt.z <- output_data$kurtosis/sqrt(24/output_data$n)
  output_data
  output_data<-as.data.frame(output_data)
  a1 <-output_data%>% kable(digits=digit, format="pandoc", caption="Describe Statistics data")
  output<-round(output_data[,c("n","mean","sd","skew", "skew.z", "kurtosis","kurt.z")],digit)
  #data refor.
  output[,"skew_TF"]<- "Not"
  output[output$skew.z < 3,"skew_TF"]<- "fair.norm"
  output[output$skew.z < 2,"skew_TF"]<- "Good.norm"
  output[output$skew.z < -2,"skew_TF"]<- "fair.norm"
  output[output$skew.z < -3,"skew_TF"]<- "Not"

  output[,"kurt_TF"]<- "Not"
  output[output$kurt.z < 10,"kurt_TF"]<- "fair.norm"
  output[output$kurt.z < 7,"kurt_TF"]<- "Good.norm"
  output[output$kurt.z < -7,"kurt_TF"]<- "fair.norm"
  output[output$kurt.z < -10,"kurt_TF"]<- "Not"

  a2 <-output%>% kable(digits=digit, format="pandoc", caption="Describe Normality Kline(2011):skew<3, kurt<10, Crran,West & Finch(1997):skew<2, kurt<7")

 return(a2)
}


#왜도와 첨도에 의한 정규성 검정  데이터 테이블을 확인, 정밀한 왜도와 첨도기준
#계산 방법은 SPSS와 동일함
#descrive data 정규성 판단####
des.3<- function(input_data, digit){
  library(psych)
  library(knitr)
  input_data<-as.data.frame(input_data)
  output_data <- describe(input_data)
  N<-output_data$n
  output_data$skew.z <- output_data$skew/sqrt((6*N*((N-1))/((N-2)*(N+1)*(N+3))))
  output_data$kurt.z <- output_data$kurtosis/sqrt((24*N*(N-1)*(N-1))/((N-3)*(N-2)*(N+3)*(N-5)))
  output_data
  output_data<-as.data.frame(output_data)
  a1 <-output_data%>% kable(digits=digit, format="pandoc", caption="Describe Statistics data")

  output<-round(output_data[,c("n","mean","sd","skew", "kurtosis", "skew.z","kurt.z")],digit)
  #data nomality check making .
  output[,"skew_TF"]<- "Not"
  output[output$skew.z < 3,"skew_TF"]<- "fair.norm"
  output[output$skew.z < 1.96,"skew_TF"]<- "Good.norm"
  output[output$skew.z < -1.96,"skew_TF"]<- "fair.norm"
  output[output$skew.z < -3,"skew_TF"]<- "Not"

  output[,"kurt_TF"]<- "Not"
  output[output$kurt.z < 3,"kurt_TF"]<- "fair.norm"
  output[output$kurt.z < 1.96,"kurt_TF"]<- "Good.norm"
  output[output$kurt.z < -1.96,"kurt_TF"]<- "fair.norm"
  output[output$kurt.z < -3,"kurt_TF"]<- "Not"

  a2 <-output%>% kable(digits=digit, format="pandoc", caption="Describe Normality, |skew.Z|<1.96,|Krut.z|<1.96 ")
  cat("ref: Kline(2011):skew<3, kurt<10, Crran,West & Finch(1997):skew<2, kurt<7 ","\n",
      "|skew.Z|<1.96,|Krut.z|<1.96인 경우 정규성을 충족한다","\n",
      "H0: 정규성을 충족한다, H1:정규성을 충족하지 않는다")
  a2

  #ss<- list(a2,pairs.panels(input_data))
  #ss
}




#상관분석 데이터 ####
Cortest<- function(data){
  library(psych)
  library(tibble)
 cor.p <-corr.test(data,method="pearson")
 cor.s <-corr.test(data,method="spearman")
 #cor.k <-corr.test(data,method="kendall")

 allR <- list(cor.p, cor.s)
allR
}


cortest.panel<- function(data){
  library(psych)
  library(tibble)
  cor <-corr.test(data)
  allR <- list(cor, pairs.panels(data,pch = 21,hist.col = "Gold", main="Correaltions Plot"))
  allR
}



#PLS-SEM####
#Internal Reliability

PLS.inter.reliaility <- function(x_pls){
  #internal reliability
  library(dplyr)
  library(knitr)
  library(ggplot2)

  pls.reliability <- x_pls$unidim
  colnames(pls.reliability)=c("Mode","MVs","Cronbach.alpha","C.R(=DG.rho)","Eig.1","Eig.2")
  a1 <-  pls.reliability%>%  kable(format = "pandoc", digits = 3,caption="Cronbach's alpha(>0.7) & Composite reliability(CR)=DG's rho(>0.7)")

  # check outer model
  a2 <- x_pls$outer_model%>% mutate("> 0.7"= ifelse(loading > 0.7,"Y","Problem"))%>% kable(format = "pandoc", digits = 3, caption = "지표 신뢰도(loading > 0.7, communality > 0.5)")



  # barchart of loadings
  g <-ggplot(data = x_pls$outer_model,
             aes(x = name, y = loading, fill = block)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    # threshold line (to peek acceptable loadings above 0.7)
    geom_hline(yintercept = 0.7, color = 'gray50') +
    # add title
    ggtitle("Barchart of Loadings") +
    # rotate x-axis names
    theme(axis.text.x = element_text(angle = 90))

  all<-list(a1,a2,g)

  all
}

#convegent validity
PLS.Convergent.Validity <- function(x_boot){
  library(knitr)
  library(dplyr)
  #t value t>1.96 --> ok
  #loding signification
  x_boot$boot$loadings$t_value <- with(x_boot$boot$loadings, Original/Std.Error)
  bootloading<- x_boot$boot$loadings
  manifest_Vaables <-rownames(bootloading)
  bootloading<- bootloading %>%
    filter(t_value > 0) %>%
    mutate(stars=ifelse(t_value > 3.29,"***",
                        ifelse(t_value > 2.56,"**",
                               ifelse(t_value> 1.96,"*",""))))
  bootloading<- cbind(manifest_Vaables,bootloading)
  a1 <-bootloading %>% kable(format = "pandoc", digits = 3, caption = "Convergent Validity-01:Indicator Loading with Bootsraping(>0.7), t > 1.96")


  #Structure.V
  x_boot$inner_summary
  #AVE
  Structure.V <-rownames(x_boot$inner_summary)
  plsAVE <-x_boot$inner_summary%>% mutate(Accept=ifelse(AVE>0.5,"Accept","No"))

  rel1 <- cbind(Structure.V ,plsAVE)
  a2 <-rel1 %>%kable(format = "pandoc", digits = 3, caption = "Convergent Validity-02: (AVE>0.5) ")


  all <- list(a1,a2)
  all
}


#discriminant Validity
PLS.discriminant.Validity<- function(x_pls){
  #판별 타당성
  # check cross loadings####
  a1 <-  x_pls$crossloadings %>% kable(format = "pandoc", digits = 3, caption = "Cross loadings")

  # load ggplot2 and reshape
  library(ggplot2)
  library(reshape)

  # reshape crossloadings data.frame for ggplot
  xloads = melt(x_pls$crossloadings, id.vars = c("name", "block"),
                variable_name = "LV")

  # bar-charts of crossloadings by block
  g <-ggplot(data = xloads,
             aes(x = name, y = value, fill = block)) +
    # add horizontal reference lines
    geom_hline(yintercept = 0, color = "gray75") +
    geom_hline(yintercept = 0.5, color = "gray70", linetype = 2) +
    # indicate the use of car-charts
    geom_bar(stat = 'identity', position = 'dodge') +
    # panel display (i.e. faceting)
    facet_wrap(block ~ LV) +
    # tweaking some grahical elements
    theme(axis.text.x = element_text(angle = 90),
          line = element_blank(),
          plot.title = element_text(size = 12)) +
    # add title
    ggtitle("Crossloadings")




  #판별타당도 2nd####
  #x_pls$scores

  # Gefen-Straub(2005): AVE의 제곱근이 그 잠재변수와 다른  잠재변수 사이의 상관계수보다 높아야 함####
  la.score = x_pls$scores
  la.cor = cor(la.score, use = "complete.obs", method = "pearson");la.cor
  sqr.AVE = with(x_pls$inner_summary, sqrt(AVE));sqr.AVE

  ##lower triangle matrix______________
  la.cor_lower <- la.cor  # 변환을 위해 기존것을 남기고 하나 더 복사
  upper.tri(la.cor_lower, diag=FALSE) # 0으로 만들 것들을 판단
  la.cor_lower[upper.tri(la.cor_lower, diag=FALSE)] <- c(0) # True를 모두 0으로 넣기
  round(la.cor_lower,3)
  #_____________________________________
  cbind(la.cor_lower, sqr.AVE)
  AVE2<- cbind(la.cor_lower, sqr.AVE)
  AVE3 <- round(AVE2,3)    # sqr.AVE  >  pearson coeff compara   ===> OK
  #round(AVE3,3)
  a2<-AVE3 %>% kable(format = "pandoc", digits = 3, caption = "Discriminant Validity: Gefen-Strab(2005)- correlation coeff(rho) < sqrt(AVE)")

  rel <- list(a1,a2,g)
  rel

}


#PLS bootstrap total effect paper form and visualization###
PLS_boot.effect<- function(edu_val){
  #bootsrap effect signification
  edu_val$boot$total.efs <- edu_val$boot$total.efs %>% filter(Original>0)


  path.tvalue = with(edu_val$boot$total.efs, Original/Std.Error)
  #path.tvalue
  cbind(edu_val$boot$total.efs, path.tvalue)
  boottotal.efs<- round(cbind(edu_val$boot$total.efs, path.tvalue),3)
  boottotal.efs
  #boottotal.efs[is.na(boottotal.efs)]<-"0"

  #t.value에 star * 표시
  boottotal.efs[is.na(boottotal.efs)]<- "0"
  boottotal.efs$star <- ""  #변수 생성
  boottotal.efs[boottotal.efs$path.tvalue > 1.96,"star"]= "*"
  boottotal.efs[boottotal.efs$path.tvalue > 2.58,"star"]= "**"
  boottotal.efs[boottotal.efs$path.tvalue > 3.29 ,"star"]= "***"
  boottotal.efs
  #t-value에 선택에 ㄷ애한 것
  #boottotal.efs
  boottotal.efs[,"Accept"]="No"
  boottotal.efs[boottotal.efs$path.tvalue >= 1.96,"Accept"]="Yes"

  a2<- boottotal.efs %>% kable(format = "pipe", caption = "Bootstrap-직접, 간접, 총효과의 유의성 검정 ")
  a1<- edu_val$effects %>% filter(total>0) %>%
    kable(format="pandoc", digits = 3,caption = "직접효과, 간접효과, 총효과")

  res <- list(a1, a2)
  res
}


#effect barplot and paper form
PLS.effectBar<- function(x_pls){
  #논문에 넣기 위해 엑셀로 편집하기 쉽도록 구성한 표
  library(dplyr)

  a1 <-x_pls$effects %>% filter(total>0) %>% kable(format ="pandoc",digits =3,
                                                   caption = "Total Effect(Direct & Indirect)")


  #path_effs<- x_pls$effects %>% filter(total>0) %>% kable(format ="pipe",
   #                                                       digits=3,
    #                                                      caption = "Total Effect(Direct & Indirect)")


  path_effs <- x_pls$effects %>% filter(total>0)
  path_effs1 <- as.matrix(path_effs[,2:3])
  # add rownames to path_effs
  rownames(path_effs1) = path_effs[,1]

  # setting margin size
  op = par(mar = c(8, 3, 1, 0.5))
  # barplots of total effects (direct + indirect)
 barplot(t(path_effs1), border = NA, col = c("#9E9AC8", "#DADAEB"),
              las = 2, cex.names = 0.8, cex.axis = 0.8,
              legend = c("Direct", "Indirect"),
              args.legend = list(x = "top", ncol = 2, border = NA,
                                 bty = "n", title = "Effects"))
  # resetting default margins
  par(op)

  return(a1)

}

#path coefficeint siginification

#path coefficeint siginification
PLS_boot.paths <- function(x_pls){
  #경로계수의 유의성
  path.tvalue = with(x_pls$boot$paths, Original/Std.Error)
  #path.tvalue
  #cbind(x_pls$boot$paths, path.tvalue)
  boot.paths<- round(cbind(x_pls$boot$paths, path.tvalue),3)
  boot.paths
  #boot.paths[is.na(boot.paths)]<-"0"

  #t.value에 star * 표시
  boot.paths[is.na(boot.paths)]<- "0"
  boot.paths$star <- ""  #변수 생성
  boot.paths[boot.paths$path.tvalue > 1.96,"star"]= "*"
  boot.paths[boot.paths$path.tvalue > 2.58,"star"]= "**"
  boot.paths[boot.paths$path.tvalue > 3.29 ,"star"]= "***"
  boot.paths
  #t-value에 선택에 ㄷ애한 것
  #boot.paths
  boot.paths[,"Accept"]="No"
  boot.paths[boot.paths$path.tvalue >= 1.96,"Accept"]="Yes"

  res<-boot.paths %>% kable(format = "pipe", caption = "Bootstrap 경로계수의 유의성 검정 ")

  plot(x_pls,arr.pos = 0.3, box.col = "skyblue", cex.txt = 1,
       txt.col = "black", main="Bootstrap Paths coefficient by Park Joonghee")
  res

  return(res)
}

#Bootstrap를 안한경우 경로계수와 효과(effects)
PLS.path.coef <-function(x_pls){
  #bootstrap를 안한 경우 경로계수와 effect
  #경로계수 나타내기
  library(dplyr)
  library(knitr)
  a1 <- x_pls$effects %>% filter(total>0) %>% filter(direct>0)  %>%select(relationships ,path.coeff=direct) %>%
    kable(format="pandoc", digits = 3,caption="direct = 경로계수")

  #effects
  a2<- x_pls$effects %>% filter(total>0) %>%
    mutate(구분=ifelse(indirect == 0,".",ifelse(direct ==0, "간접효과","."))) %>%
    kable(format="pandoc", digits = 3,caption="직접효과,간접효과, 총효과 ")

  #regression 유의성 검정
  options(scipen = 100)
  a3 <- x_pls$inner_model

  #x_pls$path_coefs
  #library(Diagram)
  #par(mar=c(1,1,1,1))
  #x_pls$path_coefs %>% round(4) %>% plotmat(box.size = 0.08, box.prop = 0.5,box.col = "Gold",box.lcol = 0,
  #                         cex.txt = 1, curve = 0.01, arr.pos = 0.45, arr.lcol = "gray70", arr.col = "gray60",
  #                                             shadow.size = 0)
  res <- list(a1,a2)

  plot(x_pls,arr.pos = 0.35, box.col = "Gold", cex.txt = 1,
       txt.col = "black", main="Regression coefficient by Park Joonghee")
  res
}


#부트스트랩 없는 효과, 경로 검정
PLS.path.coef_sig <-function(x_pls){
  #bootstrap를 안한 경우 경로계수와 effect
  #경로계수 나타내기
  library(dplyr)
  library(knitr)

  a1 <- x_pls$effects %>% filter(total>0) %>% filter(direct>0)  %>%select(relationships ,path.coeff=direct) %>%
    kable(format="pandoc", digits = 3, caption="direct = 경로계수")

  #effects
  a2<- x_pls$effects %>% filter(total>0) %>%
    mutate(구분=ifelse(indirect == 0,".",ifelse(direct ==0, "간접효과","direct"))) %>%
    kable(format="pandoc", digits = 3,caption="직접효과,간접효과, 총효과 ")

  #regression 유의성 검정
  options(scipen = 1)
  a3 <- x_pls$inner_model

  #x_pls$path_coefs
  #library(Diagram)
  #par(mar=c(1,1,1,1))
  #x_pls$path_coefs %>% round(4) %>% plotmat(box.size = 0.08, box.prop = 0.5,box.col = "Gold",box.lcol = 0,
  #                         cex.txt = 1, curve = 0.01, arr.pos = 0.45, arr.lcol = "gray70", arr.col = "gray60",
  #                                             shadow.size = 0)
  res <- list(a1,a2,a3)

  plot(x_pls,arr.pos = 0.35, box.col = "Gold", cex.txt = 1,
       txt.col = "black", main="Regression coefficient by Park Joonghee")
  res
}


#멀티그룹비교 boostratp t-test 시각화
PLS.groupcompare_ttest <- function(x_ttest){
  a1 <- x_ttest
  a2<- x_ttest$test %>% kable(format="pandoc", digits = 3)
  res <- list(a1,a2)

  barplot(t(as.matrix(x_ttest$test[,2:3])), beside = TRUE,
          col = c("#FEB24C","#74A9CF"), las = 2, ylim = c(-0.1, 1),
          cex.names = 0.8, col.axis = "gray30", cex.axis = 0.8)
  abline(h=0, col="gray50")
  title(" Path coefficient of Female and Male ")
  legend("top", legend=c("femal","male"),col = c("orange","skyblue"),ncol=2,bty="n",pch=22,
         pt.bg = c("orange","skyblue"))
  res
}



#경로도표


PLS.str_plot <- function(x){
  op1=par(mar=c(1,1,1,1), oma=c(0,0,0,0)+0.1)
  # plotting results (inner model)
  # plot(x)
  # matrix of path coefficients
  #x$path_coefs

  # plotting results (inner model):위치 조절
  #plot(x, arr.pos = 0.35)

  # matrix of path coefficients
  Paths = x$path_coefs
  # matrix with values based on path coeffs
  arrow_lwd = 10 * round(Paths, 2)
  # how does it look like?
  arrow_lwd
  #par(mfrow=c(1,1))
  # arrows of different sizes reflecting the values of the path coeffs
  plot(x, arr.pos = 0.35, arr.lwd = arrow_lwd, cex.txt = 1,
       box.size = 0.1, box.col = "gray85", box.prop = 0.3, box.cex = 1.5, txt.col = "black")
}



PLS.outer_plot <- function(x){

  # plotting results (inner model)
  #  plot(x)
  # matrix of path coefficients
  # x$path_coefs

  # plotting results (inner model):위치 조절
  #plot(x, arr.pos = 0.35)

  # matrix of path coefficients
  Paths = x$path_coefs
  # matrix with values based on path coeffs
  arrow_lwd = 10 * round(Paths, 2)
  # how does it look like?
  arrow_lwd
  #par(mfrow=c(1,1))
  # arrows of different sizes reflecting the values of the path coeffs

  plot(x, "loadings",arr.lwd = arrow_lwd, cex.txt = 1,
       box.size = 0.1, box.col = "Gold", box.prop = 0.5, box.cex = 1.5, txt.col = "black")
}

PLS.weights_plot <- function(x){

  # plotting results (inner model)
  #  plot(x)
  # matrix of path coefficients
  # x$path_coefs

  # plotting results (inner model):위치 조절
  #plot(x, arr.pos = 0.35)

  # matrix of path coefficients
  Paths = x$path_coefs
  # matrix with values based on path coeffs
  arrow_lwd = 10 * round(Paths, 2)
  # how does it look like?
  arrow_lwd
  #par(mfrow=c(1,1))
  # arrows of different sizes reflecting the values of the path coeffs

  plot(x, "weights",arr.lwd = arrow_lwd, cex.txt = 1,
       box.size = 0.1, box.col = "Gold", box.prop = 0.5, box.cex = 1.5, txt.col = "black")
}












# latent growth modeling - basicL variable define -  Slope, Inter

lgm <- function(x){
  library(knitr)
  #estimate
  parameterEstimates(x,standardized = T)


  #filter latent growth data
  parameterEstimates(x,standardized = T) %>% filter(op=="~~"|op=="~1") %>% filter(est>0)
  #define data
  para.data <- parameterEstimates(x,standardized = T) %>% filter(op=="~~"|op=="~1") %>% filter(est>0)



  #Extraction Slope
  SLOPTE<- para.data %>% filter(lhs=="Slope") %>% filter(op=="~1")%>% select(est)

  #Extraction intercept
  INTERCEPT<- para.data %>% filter(lhs=="Inter")%>% filter(op=="~1")%>% select(est)


  #Extraction Slope
  SLOPTE.1<- para.data %>% filter(lhs=="Slope") %>% filter(op=="~1")%>% select(est, pvalue)

  #Extraction intercept
  INTERCEPT.1<- para.data %>% filter(lhs=="Inter")%>% filter(op=="~1")%>% select(est, pvalue)



  #slope  intercept
  cov_T_S <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~~") %>% filter(rhs=="Slope") %>% select(est,pvalue)

  #resute
  data_1<- cbind(INTERCEPT.1,SLOPTE.1,cov_T_S)
  colnames(data_1)=c("Intercepts","p.value(Inter)", "Slope", "p.value(Slope)","Cov(I<->S)","p.value(Cov)")
  data_2 <- data_1 %>% kable(digits = 3, format = "pandoc", caption="Slope & intercept & Covariance sig.")

  #extraction 0,1,2,3 time ####
  TIME <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Slope") %>% filter(op=="=~") %>% select(est)

  #1~4 ; error
  ERROR <- para.data %>% slice(1:4) %>% filter(op=="~~")%>% select(est)
  varname <- c("1st","2nd","3rd","4st")

  #cbind data
  DATA.TOTAL <- cbind(varname , INTERCEPT,TIME ,SLOPTE,  ERROR)

  DATA.TOTAL<- as.data.frame(DATA.TOTAL)
  DATA.TOTAL<-DATA.TOTAL[,1:5]
  #names
  colnames(DATA.TOTAL)=c("Measure","Intercepts", "time","Slope","error")

  #DATA.TOTAL
  DATA.TOTAL_1 <- DATA.TOTAL %>% mutate(Mean_with_Error=Intercepts+time*Slope+error) %>% mutate(Mean_No_Error=Intercepts+time*Slope)%>%  kable(digits=3, format="pandoc", caption="latent growth data and mean result")

  #des.data <- des(x)[[3]]
  #des.data1 <-des.data[,c("vars","mean")]

  all.r <- list(data_2, DATA.TOTAL_1 )
  all.r

}

#LGM mean plot#####

lgm.plot <- function(x){

  #estimate
  parameterEstimates(x,standardized = T)


  #filter latent growth data
  parameterEstimates(x,standardized = T) %>% filter(op=="~~"|op=="~1") %>% filter(est>0)
  #define data
  para.data <- parameterEstimates(x,standardized = T) %>% filter(op=="~~"|op=="~1") %>% filter(est>0)

  #Extraction Slope
  SLOPTE<- para.data %>% filter(lhs=="Slope") %>% filter(op=="~1")%>% select(est)

  #Extraction intercept
  INTERCEPT<- para.data %>% filter(lhs=="Inter")%>% filter(op=="~1")%>% select(est)

  #slope  intercept
  cov_T_S <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~~") %>% filter(rhs=="Slope") %>% select(est,std.all,pvalue)

  #resute
  data_1<- cbind(INTERCEPT,SLOPTE,cov_T_S)
  colnames(data_1)=c("Intercepts", "Slope","covariance","correlation","p.value")
  data_2 <- data_1 %>% kable(digits = 3, format = "pandoc", cpation="Slope & intercept & Covariance")

  #extraction 0,1,2,3 time ####
  TIME <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Slope") %>% filter(op=="=~") %>% select(est)

  #1~4 ; error
  ERROR <- para.data %>% slice(1:4) %>% filter(op=="~~")%>% select(est)
  varname <- para.data %>% slice(1:4) %>% filter(op=="~~")%>% select(lhs)

  #cbind data
  DATA.TOTAL <- cbind(varname , INTERCEPT,TIME ,SLOPTE,  ERROR)

  DATA.TOTAL<- as.data.frame(DATA.TOTAL)
  DATA.TOTAL<-DATA.TOTAL[,1:5]
  #names
  colnames(DATA.TOTAL)=c("Measure","Intercepts", "time","Slope","error")

  #DATA.TOTAL
  DATA.TOTAL_1 <- DATA.TOTAL %>% mutate(result_with_Error=Intercepts+time*Slope+error) %>% mutate(result_No_Error=Intercepts+time*Slope)

  #des.data <- des(x)[[3]]
  #des.data1 <-des.data[,c("vars","mean")]

  plot(result_No_Error ~ time, data=DATA.TOTAL_1,type="b", col="red", pch=16,cex=1.5,main="잠재성장모델 분석결과 예측된 평균변화 plot" )
  abline(v=0,lty=3,col=3,lwd=1)
  abline(v=1,lty=3,col=4,lwd=1)
  abline(v=2,lty=3,col=5,lwd=1)
  abline(v=3,lty=3,col=6,lwd=1)
  text(0.1,DATA.TOTAL_1[1,7],"1st")
  text(1.1,DATA.TOTAL_1[2,7],"2nd")
  text(2.1,DATA.TOTAL_1[3,7],"3rd")
  text(2.9,DATA.TOTAL_1[4,7],"4st")

}


#one variable lgm #####

lgm_1 <-function(x,Gen_name){
  library(knitr)
# #################################

#estimate
#parameterEstimates(x,standardized = T)

#filter latent growth data
#parameterEstimates(x,standardized = T) %>% filter(op=="~~"|op=="~1") %>% filter(est>0)
#define data
para.data1 <- parameterEstimates(x,standardized = T) %>% filter(op=="~~"|op=="~1") %>% filter(est>0)
#para.data1
#Extraction Slope
SLOPTE<- para.data1 %>% filter(lhs=="Slope") %>% filter(op=="~1")%>% select(est)
#SLOPTE
#Extraction intercept
INTERCEPT <- para.data1 %>% filter(lhs=="Inter")%>% filter(op=="~1")%>% select(est)
#INTERCEPT

SLOPTE.1<- para.data1 %>% filter(lhs=="Slope") %>% filter(op=="~1")%>% select(est,pvalue)
#SLOPTE
#Extraction intercept
INTERCEPT.1 <- para.data1 %>% filter(lhs=="Inter")%>% filter(op=="~1")%>% select(est,pvalue)
#INTERCEPT



#slope  intercept
cov_T_S <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~~") %>% filter(rhs=="Slope") %>% select(est)
#cov_T_S

#slope  intercept
cov_T_S.1 <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~~") %>% filter(rhs=="Slope") %>% select(est,std.all ,pvalue)
#cov_T_S.1


#cov_G_N<- parameterEstimates(x,standardized = T) %>% filter(lhs=="Gender") %>% filter(op=="~~") %>% filter(rhs=="NR") %>% select(est)
#cov_G_N

#cov_G_N.1<- parameterEstimates(x,standardized = T) %>% filter(lhs=="Gender") %>% filter(op=="~~") %>% filter(rhs=="NR") %>% select(est, pvalue)
#cov_G_N.1

#resute
data_1 <- cbind(INTERCEPT.1,SLOPTE.1,cov_T_S.1)
colnames(data_1)=c("Intercepts","p (Inter)", "Slope","p (Slope)", " Inter <-> Slope","cor(I-S)","p-value(I_S)")
data_2 <- data_1 %>% kable(digits = 3, format = "pandoc", caption="Slope & intercept & Covariance")

#data_1
#data_2  #결과




#extraction 0,1,2,3 time ####
TIME <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Slope") %>% filter(op=="=~") %>% select(est)
##TIME


#1~4 ; error
ERROR <- para.data %>% slice(1:4) %>% filter(op=="~~")%>% select(est)
varname <- c("1st","2nd","3rd","4st")
#ERROR
#varname


#parameterEstimates(x,standardized = T)

#regression
Inter_Gen <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~") %>% filter(rhs==Gen_name) %>% select(est)
Slope_Gen <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Slope") %>% filter(op=="~") %>% filter(rhs==Gen_name)%>% select(est)

#Inter_NR <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~") %>% filter(rhs=="NR")%>% select(est)
#Slope_NR <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Slope") %>% filter(op=="~") %>% filter(rhs=="NR")%>% select(est)



#Inter_Gen
#Slope_Gen

Gen_intercept <- parameterEstimates(x,standardized = T) %>% filter(lhs==Gen_name) %>% filter(op=="~1")%>% select(est)
#NR_intercept <- parameterEstimates(x,standardized = T) %>% filter(lhs=="NR") %>% filter(op=="~1")%>% select(est)

#Gen_intercept
#NR_intercept

#cbind data
DATA.TOTAL.m <- cbind(varname, INTERCEPT ,TIME ,SLOPTE,  ERROR, Inter_Gen,Slope_Gen,Gen_intercept)

#DATA.TOTAL.m

#DATA.TOTAL.m<- as.data.frame(DATA.TOTAL.m)
#DATA.TOTAL.m<-DATA.TOTAL.m[,1:5]

#names
colnames(DATA.TOTAL.m)=c("name", "Intercepts", "time","Slope","error","Grp_inter","Grp_slope","Grp_intercept")
#DATA.TOTAL.m

#plot ####
#dt_1 <- DATA.TOTAL.m %>% mutate(Resuilt_Mean = (Intercepts + Gen_inter*Gen_intercept) + (time)*(Slope+ Gen_slope*Gen_intercept ) +  error)
#plot(Resuilt_Mean ~ time, dt_1, type="b", pch=16,cex=1.5,  col="red", main="평균변화")
#dt_1


#DATA.TOTAL_1 ####
DATA.TOTAL_1 <- DATA.TOTAL.m %>% mutate(Mean_with_Error = (Intercepts + Grp_inter*Grp_intercept ) + time*(Slope + Grp_slope*Grp_intercept ) +  error) %>% mutate(Mean_No_error = (Intercepts + Grp_inter*Grp_intercept ) + time*(Slope + Grp_slope*Grp_intercept )) %>% kable(digits=3, format="pandoc", caption="latent growth data of Group variable")

#DATA.TOTAL_1

result=list(data_2 , DATA.TOTAL_1)
result
}


#one variable lgm #####
#https://blog.naver.com/shoutjoy/222031547222
lgm_1.plot <-function(x,Gen_name){
  library(knitr)


  #estimate
  #parameterEstimates(x,standardized = T)

  #filter latent growth data
  #parameterEstimates(x,standardized = T) %>% filter(op=="~~"|op=="~1") %>% filter(est>0)
  #define data
  para.data1 <- parameterEstimates(x,standardized = T) %>% filter(op=="~~"|op=="~1") %>% filter(est>0)
  #para.data1
  #Extraction Slope
  SLOPTE<- para.data1 %>% filter(lhs=="Slope") %>% filter(op=="~1")%>% select(est)
  #SLOPTE
  #Extraction intercept
  INTERCEPT <- para.data1 %>% filter(lhs=="Inter")%>% filter(op=="~1")%>% select(est)
  #INTERCEPT

  SLOPTE.1<- para.data1 %>% filter(lhs=="Slope") %>% filter(op=="~1")%>% select(est,pvalue)
  #SLOPTE
  #Extraction intercept
  INTERCEPT.1 <- para.data1 %>% filter(lhs=="Inter")%>% filter(op=="~1")%>% select(est,pvalue)
  #INTERCEPT



  #slope  intercept
  cov_T_S <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~~") %>% filter(rhs=="Slope") %>% select(est)
  #cov_T_S

  #slope  intercept
  cov_T_S.1 <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~~") %>% filter(rhs=="Slope") %>% select(est, pvalue)
  #cov_T_S.1


  #cov_G_N<- parameterEstimates(x,standardized = T) %>% filter(lhs=="Gender") %>% filter(op=="~~") %>% filter(rhs=="NR") %>% select(est)
  #cov_G_N

  #cov_G_N.1<- parameterEstimates(x,standardized = T) %>% filter(lhs=="Gender") %>% filter(op=="~~") %>% filter(rhs=="NR") %>% select(est, pvalue)
  #cov_G_N.1

  #resute
  data_1 <- cbind(INTERCEPT.1,SLOPTE.1,cov_T_S.1)
  colnames(data_1)=c("Intercepts","p (Inter)", "Slope","p (Slope)", " Inter <-> Slope","p-value(I_S)")
  data_2 <- data_1 %>% kable(digits = 3, format = "pandoc", caption="Slope & intercept & Covariance")

  #data_1
  #data_2  #결과




  #extraction 0,1,2,3 time ####
  TIME <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Slope") %>% filter(op=="=~") %>% select(est)
  ##TIME


  #1~4 ; error
  ERROR <- para.data %>% slice(1:4) %>% filter(op=="~~")%>% select(est)
  varname <- c("1st","2nd","3rd","4st")
  #ERROR
  #varname


  #parameterEstimates(x,standardized = T)

  #regression
  Inter_Gen <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~") %>% filter(rhs==Gen_name) %>% select(est)
  Slope_Gen <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Slope") %>% filter(op=="~") %>% filter(rhs==Gen_name)%>% select(est)

  #Inter_NR <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~") %>% filter(rhs=="NR")%>% select(est)
  #Slope_NR <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Slope") %>% filter(op=="~") %>% filter(rhs=="NR")%>% select(est)



  #Inter_Gen
  #Slope_Gen

  Gen_intercept <- parameterEstimates(x,standardized = T) %>% filter(lhs==Gen_name) %>% filter(op=="~1")%>% select(est)
  #NR_intercept <- parameterEstimates(x,standardized = T) %>% filter(lhs=="NR") %>% filter(op=="~1")%>% select(est)

  #Gen_intercept
  #NR_intercept

  #cbind data
  DATA.TOTAL.m <- cbind(varname, INTERCEPT ,TIME ,SLOPTE,  ERROR, Inter_Gen,Slope_Gen,Gen_intercept)

  #DATA.TOTAL.m

  #DATA.TOTAL.m<- as.data.frame(DATA.TOTAL.m)
  #DATA.TOTAL.m<-DATA.TOTAL.m[,1:5]

  #names
  colnames(DATA.TOTAL.m)=c("name", "Intercepts", "time","Slope","error","Grp_inter","Grp_slope","Grp_intercept")
  #DATA.TOTAL.m

  #plot ####
  #dt_1 <- DATA.TOTAL.m %>% mutate(Resuilt_Mean = (Intercepts + Gen_inter*Gen_intercept) + (time)*(Slope+ Gen_slope*Gen_intercept ) +  error)
  #plot(Resuilt_Mean ~ time, dt_1, type="b", pch=16,cex=1.5,  col="red", main="평균변화")
  #dt_1


  #DATA.TOTAL_1 ####
  DATA.TOTAL_1 <- DATA.TOTAL.m %>% mutate(Mean_with_Error = (Intercepts + Grp_inter*Grp_intercept ) + time*(Slope + Grp_slope*Grp_intercept ) +  error) %>% mutate(Mean_No_error = (Intercepts + Grp_inter*Grp_intercept ) + time*(Slope + Grp_slope*Grp_intercept ))

  #DATA.TOTAL_1


  plot(Mean_No_error ~ time, data=DATA.TOTAL_1,type="b", col="red", pch=16,cex=1.5,main="잠재성장모델 분석결과 예측된 평균변화 plot" )
  abline(v=0,lty=3,col=3,lwd=1)
  abline(v=1,lty=3,col=4,lwd=1)
  abline(v=2,lty=3,col=5,lwd=1)
  abline(v=3,lty=3,col=6,lwd=1)
  text(0.1,DATA.TOTAL_1[1,7],"1st")
  text(1.1,DATA.TOTAL_1[2,7],"2nd")
  text(2.1,DATA.TOTAL_1[3,7],"3rd")
  text(2.9,DATA.TOTAL_1[4,7],"4st")
}


# 2variable lgm #################################
#잠재성장 모델
lgm_2 <- function(x,NR_name){
  library(knitr)
  #estimate
  #parameterEstimates(x,standardized = T)

  #filter latent growth data
  #parameterEstimates(x,standardized = T) %>% filter(op=="~~"|op=="~1") %>% filter(est>0)
  #define data
  para.data1 <- parameterEstimates(x,standardized = T) %>% filter(op=="~~"|op=="~1") %>% filter(est>0)
  #para.data1
  #Extraction Slope
  SLOPTE<- para.data1 %>% filter(lhs=="Slope") %>% filter(op=="~1")%>% select(est)
  #SLOPTE
  #Extraction intercept
  INTERCEPT <- para.data1 %>% filter(lhs=="Inter")%>% filter(op=="~1")%>% select(est)
  #INTERCEPT
  #slope  intercept
  cov_T_S <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~~") %>% filter(rhs=="Slope") %>% select(est)
  #cov_T_S

  #slope  intercept
  cov_T_S.1 <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~~") %>% filter(rhs=="Slope") %>% select(est,std.all, pvalue)
  #cov_T_S.1


  cov_G_N<- parameterEstimates(x,standardized = T) %>% filter(lhs=="Gender") %>% filter(op=="~~") %>% filter(rhs==NR_name) %>% select(est)
  #cov_G_N

  cov_G_N.1<- parameterEstimates(x,standardized = T) %>% filter(lhs=="Gender") %>% filter(op=="~~") %>% filter(rhs==NR_name) %>% select(est, std.all, pvalue)
  #cov_G_N.1

  #resute
  data_1 <- cbind(INTERCEPT,SLOPTE,cov_T_S.1, cov_G_N.1)
  colnames(data_1)=c("Intercepts", "Slope", " Inter <-> Slope","cor(I-S)","p-value(I_S)","Gender <-> NR","cor(G-N)","I_S.p-value(G_N)")
  data_2 <- data_1 %>% kable(digits = 3, format = "pandoc", caption="Slope & intercept & Covariance")

  #data_2  #결과
  #extraction 0,1,2,3 time ####
  TIME <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Slope") %>% filter(op=="=~") %>% select(est)
  #TIME


  #1~4 ; error
  ERROR <- para.data %>% slice(1:4) %>% filter(op=="~~")%>% select(est)
  varname <- para.data %>% slice(1:4) %>% filter(op=="~~")%>% select(lhs)
  #ERROR
  #varname

  #
  #parameterEstimates(x,standardized = T)

  #regression
  Inter_Gen <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~") %>% filter(rhs=="Gender") %>% select(est)
  Slope_Gen <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Slope") %>% filter(op=="~") %>% filter(rhs=="Gender")%>% select(est)

  Inter_NR <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Inter") %>% filter(op=="~") %>% filter(rhs==NR_name)%>% select(est)
  Slope_NR <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Slope") %>% filter(op=="~") %>% filter(rhs==NR_name)%>% select(est)



  #Inter_Gen
  #Slope_Gen#
  #Inter_NR
  #Slope_NR

  Gen_intercept <- parameterEstimates(x,standardized = T) %>% filter(lhs=="Gender") %>% filter(op=="~1")%>% select(est)
  NR_intercept <- parameterEstimates(x,standardized = T) %>% filter(lhs==NR_name) %>% filter(op=="~1")%>% select(est)

  Gen_intercept_0 <- 0
  Gen_intercept_1 <- 1
  #Gen_intercept
  #NR_intercept

  DATA.TOTAL.grp0 <- cbind(varname, INTERCEPT ,TIME ,SLOPTE,  ERROR, Inter_Gen,Slope_Gen,Inter_NR,Slope_NR,Gen_intercept_0,NR_intercept )

  DATA.TOTAL.grp1 <- cbind(varname, INTERCEPT ,TIME ,SLOPTE,  ERROR, Inter_Gen,Slope_Gen,Inter_NR,Slope_NR,Gen_intercept_1,NR_intercept )



  NR_intercept_0 <- 0
  NR_intercept_1 <- 1


  #cbind data
  DATA.TOTAL.NR0 <- cbind(varname, INTERCEPT ,TIME ,SLOPTE,  ERROR, Inter_Gen,Slope_Gen,Inter_NR,Slope_NR,Gen_intercept,NR_intercept_0  )

  #cbind data
  DATA.TOTAL.NR1 <- cbind(varname, INTERCEPT ,TIME ,SLOPTE,  ERROR, Inter_Gen,Slope_Gen,Inter_NR,Slope_NR,Gen_intercept,NR_intercept_1  )


  #cbind data
  DATA.TOTAL.m <- cbind(varname, INTERCEPT ,TIME ,SLOPTE,  ERROR, Inter_Gen,Slope_Gen,Inter_NR,Slope_NR,Gen_intercept,NR_intercept )

  #DATA.TOTAL.m

  #DATA.TOTAL.m<- as.data.frame(DATA.TOTAL.m)
  #DATA.TOTAL.m<-DATA.TOTAL.m[,1:5]

  #names
  colnames(DATA.TOTAL.m)=c("name", "Intercepts", "time","Slope","error","Gen_inter","Gen_slope","NR_Inter","NR_Slope","Gen_intercept","NR_intercept" )

  colnames(DATA.TOTAL.grp0)=c("name", "Intercepts", "time","Slope","error","Gen_inter","Gen_slope","NR_Inter","NR_Slope","Gen_intercept","NR_intercept" )
  colnames(DATA.TOTAL.grp1)=c("name", "Intercepts", "time","Slope","error","Gen_inter","Gen_slope","NR_Inter","NR_Slope","Gen_intercept","NR_intercept" )


  colnames(DATA.TOTAL.NR0)=c("name", "Intercepts", "time","Slope","error","Gen_inter","Gen_slope","NR_Inter","NR_Slope","Gen_intercept","NR_intercept" )
  colnames(DATA.TOTAL.NR1)=c("name", "Intercepts", "time","Slope","error","Gen_inter","Gen_slope","NR_Inter","NR_Slope","Gen_intercept","NR_intercept" )


  #DATA.TOTAL.m

  #plot ####
  #dt_1 <- DATA.TOTAL.m %>% mutate(Resuilt_Mean = (Intercepts + Gen_inter*Gen_intercept + NR_Inter*NR_intercept) + time*(Slope+ Gen_slope*Gen_intercept +  NR_Slope*NR_intercept) +  error)
  #plot(Resuilt_Mean ~ time, dt_1, type="b", pch=16,cex=1.5,  col="red", main="평균변화")
  #dt_1
  DATA.TOTAL_grp0 <- DATA.TOTAL.grp0 %>% mutate(grp_0_Mean = (Intercepts + Gen_inter*Gen_intercept + NR_Inter*NR_intercept) + time*(Slope+ Gen_slope*Gen_intercept +  NR_Slope*NR_intercept) +  error) %>% kable(digits=3, format="pandoc", caption="latent growth data of Gender ==0 (female)")


  DATA.TOTAL_grp1 <- DATA.TOTAL.grp1 %>% mutate(grp_1_Mean = (Intercepts + Gen_inter*Gen_intercept + NR_Inter*NR_intercept) + time*(Slope+ Gen_slope*Gen_intercept +  NR_Slope*NR_intercept) +  error) %>% kable(digits=3, format="pandoc", caption="latent growth data of Gender == 1 (male)")




  DATA.TOTAL_NR0 <- DATA.TOTAL.NR0 %>% mutate(NR_0_Mean = (Intercepts + Gen_inter*Gen_intercept + NR_Inter*NR_intercept) + time*(Slope+ Gen_slope*Gen_intercept +  NR_Slope*NR_intercept) +  error) %>% kable(digits=3, format="pandoc", caption="latent growth data of NR ==0 (CASE-1")


  DATA.TOTAL_NR1 <- DATA.TOTAL.NR1 %>% mutate(NR_1_Mean = (Intercepts + Gen_inter*Gen_intercept + NR_Inter*NR_intercept) + time*(Slope+ Gen_slope*Gen_intercept +  NR_Slope*NR_intercept) +  error) %>% kable(digits=3, format="pandoc", caption="latent growth data of NR == 1 (CASE-2)")






  #DATA.TOTAL_1 ####
  DATA.TOTAL_1 <- DATA.TOTAL.m %>% mutate(Result_Mean = (Intercepts + Gen_inter*Gen_intercept + NR_Inter*NR_intercept) + time*(Slope+ Gen_slope*Gen_intercept +  NR_Slope*NR_intercept) +  error) %>% kable(digits=3, format="pandoc", caption="latent growth data of Gender and NR(New ratio)")

  result <- list(data_2,  DATA.TOTAL_grp0,  DATA.TOTAL_grp1,DATA.TOTAL_NR0 ,DATA.TOTAL_NR1, DATA.TOTAL_1)
  result
}

#정준상관분석####


#정준상관분석 ####
cancor.plot <- function(data,r,c){


  X=data[r]
  X1 <- scale(X, scale=F)
  #class(X)

  Y=data[c]
  Y1 <- scale(Y, scale=F)


  # cancor()
  # X,Y는 표준화한 데이터 사용 X1,Y1
  cancor.head <- cancor(X1,Y1)



  #row and column coordinates
  Rx <- X1%*%cancor.head$xcoef
  Cx <- cancor.head$xcoef
  Ry <- Y1%*%cancor.head$ycoef
  Cy <- cancor.head$ycoef
  colnames(Cx) <- colnames(X1)
  colnames(Cy) <- colnames(Y1)

  x11()
  par(mfrow=c(1,2))

  #first son
  biplot(Rx,Cx,cex=1, pch=16,
         xlim=c(-0.5,0.5), ylim = c(-0.5,0.5),
         xlab = "1st Dimension", ylab = "2nd Dimension", main = "(a) first Variable")
  abline(v=0,h=0,lty=2, col=4)

  #second son
  biplot(Ry,Cy,cex=1, pch=16,
         xlim=c(-0.5,0.5), ylim = c(-0.5,0.5),
         xlab = "1st Dimension", ylab = "2nd Dimension", main = "(b) second Variable")
  abline(v=0,h=0,lty=2, col=4)
  par(mfrow=c(1,1))



  return( cancor.head )
}
#실행방법
#cancor.plot(data[,1:2], data[,3:4])

#library(MVT)
#data(examScor)
#examScor[,1:2]
#examScor[,3:5]
#02
#cancor.plot(examScor[,1:2], examScor[,3:5])
#CCA TEST significant
cca_test <- function(data,r,c){
  library(knitr)
  # Sets of Variables : X, Y
  X1=data[,r] # Closed books
  X=scale(X1, scale=T)
  Y1=data[,c] # Opened books
  Y=scale(Y1, scale=T)
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)

  #[Step 2] Covariance Matrix S(or Correlation Matix R)
  R=round(cor(data),3)
  #R

  Rxx=R[r, r]
  Ryy=R[c, c]
  Rxy=R[r, c]
  Ryx=t(Rxy)

  # Rx ^ Ry
  Exx <- eigen(Rxx)
  Eyy <- eigen(Ryy)
  Rx <- Exx$vectors %*% diag(sqrt(Exx$values)) %*% t(Exx$vectors)
  Ry <- Eyy$vectors %*% diag(sqrt(Eyy$values)) %*% t(Eyy$vectors)

  #[Step 3] Spectral Decomposition : M=PDP
  M=solve(Rx)%*%Rxy%*%solve(Ryy)%*%Ryx%*%solve(Rx)
  #round(M, 3)
  eigen.M=eigen(M)
  eig=eigen.M$values
  round(eig, 3) # Eigenvalues
  rho=round(sqrt(eig), 3) #Canonical Correlation


  #[data analysis ]
  P=eigen.M$vectors
  U=round(solve(Rx)%*%P, 3) # Eigenvectors
  rownames(U)<-colnames(X)
  #U  #제1 정준변수쌍 axis 1

  V=solve(Ryy)%*%Ryx%*%U%*%diag(1/sqrt(eig))
  #V   #제2 정준변수쌍 axis 2

  #[Step 4] Canonical Varables Scores
  Zx=X%*%U
  #Zx
  Zy=Y%*%V
  #Zy



  #[Step 4] Testing Significant Canonical Correlations
  ev <- (1 - eigen.M$values)
  s <- min(p, q)
  #s
  w=n - 3/2 - (p + q)/2
  Lambda <- rev(cumprod(rev(ev)))
  # initialize
  d1 <- d2 <- f <- vector("numeric", s)
  for (i in 1:s) {
    t <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
    ti <- 1/t
    d1[i] <- p * q
    d2[i] <-w * t - p * q/2 + 1
    r <- (1 - Lambda[i]^ti)/Lambda[i]^ti
    f[i] <- r * d2[i]/d1[i]
    p <- p - 1
    q <- q - 1
  }

  p_value <- pf(f, d1, d2, lower.tail = FALSE)

  cancor.s <-cancor(X,Y)
  cancor.cor <- cancor.s$cor

  dmat <- cbind(Eigen_Value=eig,Canonical_Corr = cancor.cor,WilksLambda = Lambda, F = f, df_1 = d1, df_2 = d2, p_value = p_value)
  rownames(dmat)=c("1st_axis", "2nd_axis")

  # scienct numeric

  damat.re <-kable(dmat, format = "pandoc", digits=3, caption = "CCA TEST F- significant")


  cca_result<- rbind(U,V,Eigen_Value=eig,Canonical_Correlation=rho,WilksLambda = Lambda, F = f, df_1 = d1, df_2 = d2, p_value = p_value)
  colnames(cca_result)=c("1st pair", "2nd pair")
  secondResult<- kable(cca_result, format = "pandoc", digits = 3, caption = "Canonical Variables & Correlations & sigification")





  Re<- list(damat.re, secondResult)
  Re
}
cca_Matirx_test <- function(data,r,c,sample_N){
  library(knitr)
  # Sets of Variables : X, Y
  #   X=data[,r] # Closed books
  #  X=scale(X, scale=T)
  # Y=data[,c] # Opened books
  #  Y=scale(Y, scale=T)
  n=sample_N
  p=length(r)
  q=length(c)

  #[Step 2] Covariance Matrix S(or Correlation Matix R)
  R=round(data,3)
  #R

  Rxx=R[r, r]
  Ryy=R[c, c]
  Rxy=R[r, c]
  Ryx=t(Rxy)

  # Rx ^ Ry
  Exx <- eigen(Rxx)
  Eyy <- eigen(Ryy)
  Rx <- Exx$vectors %*% diag(sqrt(Exx$values)) %*% t(Exx$vectors)
  Ry <- Eyy$vectors %*% diag(sqrt(Eyy$values)) %*% t(Eyy$vectors)

  #[Step 3] Spectral Decomposition : M=PDP
  M=solve(Rx)%*%Rxy%*%solve(Ryy)%*%Ryx%*%solve(Rx)
  #round(M, 3)
  eigen.M=eigen(M)
  eig=eigen.M$values
  round(eig, 3) # Eigenvalues
  rho=round(sqrt(eig), 3) #Canonical Correlation


  #[data analysis ]#######################################
  P=eigen.M$vectors

  U=round(solve(Rx)%*%P, 3) # Eigenvectors
  rownames(U)<-colnames(X)
  #U  #제1 정준변수쌍 axis 1

  V=solve(Ryy)%*%Ryx%*%U%*%diag(1/sqrt(eig))
  #V   #제2 정준변수쌍 axis 2

  #[Step 4] Canonical Varables Scores
  Zx=X%*%U
  #Zx
  Zy=Y%*%V
  #Zy





  #[Step 4] Testing Significant Canonical Correlations
  ev <- (1 - eigen.M$values)
  s <- min(p, q)
  #s
  w=n - 3/2 - (p + q)/2
  Lambda <- rev(cumprod(rev(ev)))
  # initialize
  d1 <- d2 <- f <- vector("numeric", s)
  for (i in 1:s) {
    t <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
    ti <- 1/t
    d1[i] <- p * q
    d2[i] <-w * t - p * q/2 + 1
    r <- (1 - Lambda[i]^ti)/Lambda[i]^ti
    f[i] <- r * d2[i]/d1[i]
    p <- p - 1
    q <- q - 1
  }

  p_value <- pf(f, d1, d2, lower.tail = FALSE)

  cancor.s <-cancor(X,Y)
  cancor.cor <- cancor.s$cor

  dmat <- cbind(WilksLambda = Lambda, F = f, df_1 = d1, df_2 = d2, p_value = p_value)
  rownames(dmat)=c("1st_axis", "2nd_axis")

  # scienct numeric

  damat.re <-kable(dmat, format = "pandoc", caption = "CCA TEST F- significant")


  # second paper data
  cca_result<- rbind(U,V,Eigen_Value=eig,Canonical_Correlation=rho,WilksLambda = Lambda, F = f, df_1 = d1, df_2 = d2, p_value = p_value)
  colnames(cca_result)=c("1st pair", "2nd pair")
  secondResult<- kable(cca_result, format = "pandoc", digits = 3, caption = "Canonical Variables & Correlations & sigification ")


  Re<- list(damat.re, secondResult)
  Re
}


#평균중심변환 함수2020.11.17 ####
#변수 Centering 평균중심변환 함수

#변수 Centering 평균중심변환 함수
Mean_center1 <-function(dataset,variable){
  if(!is.numeric(variable)){return("연속형변수가 아닙니다. 데이터, `데이터$변수명`을 정확히 입력해주세요 ")}

  library(dplyr)
  data<-dataset %>% as.data.frame()
  data$x<-variable
  cen.data <- data %>% mutate(x_m=mean(x), x_Center = x - mean(x),SD_x=sd(x),sd_low=x_Center+sd(x),sd_high=x_Center-sd(x)) %>% round(2)
  cn<- cen.data[,c(1,2,3,7,8)]
  #aa<-list(centered=cn)#,full_data=cen.data)
  cn
}

Mean_center <-function(dataset,x){
  #  attach(dataset)
  library(dplyr)
  if(!is.numeric(x)){return("연속형변수가 아닙니다. 데이터, `데이터$변수명`을 정확히 입력해주세요 ")}

  dataset<-dataset %>% as.data.frame()
  #dataset$x <- x
  row_name<-rownames(dataset)
  cen.data <- dataset %>% mutate(x_m=mean(x), x_Center = x - mean(x),SD_x=sd(x),sd_low=x_Center+sd(x),sd_high=x_Center-sd(x)) %>% round(2)
  rownames(cen.data)=row_name
  res<-cen.data[,c(1,2,3,5,7,8)]
  res
}

#Mean_center(women,women$weight)




#simple 기술통계 aggregate####
#평균과 표준편차를 빠르게 구하기 aggregate####
dataMSD <- function(data,form){
  cbind(Mean=aggregate(form ,data, mean),SD=aggregate(form,data, sd)[,-1]) %>% round(2)
}

#dataMSD(mtcars, mpg~am)
#dataMSD(mtcars, mpg~am+vs)
# count(mtcars,am,vs)


#평균과 표준편차, 샘플수 구하기 데이터, 형식(종속변수~그룹변수, 그룹분류(데이터프레임))
dataMSDN <- function(data,form, ...){
  #counta <-table(data$x)
  library(dplyr)
  library(knitr)

  n<-count(data, ...)
  agg_mean <-aggregate(form ,data, mean)
  casemean <- agg_mean[,ncol(n)]
  agg_sd <- aggregate(form ,data, sd)
  a<-cbind(Case=agg_mean[,-ncol(n)], M=casemean,SD=agg_sd[,ncol(n)],N=n[,ncol(n)]) %>% round(2)
  #s=list(a,n);s
  a<-a %>% as.data.frame() %>% kable("pandoc",caption = "그룹별 기술통계")
  a
}
#실행
 #dataMSDN(mtcars, mpg~am,am)
 #dataMSDN(mtcars, mpg~am+vs,am,vs)
 #dataMSDN(mtcars, mpg~am+vs+cyl,am,vs,cyl)

#form에서 나타난 만큼 독립변수를  ,로 구분하여 입력 해야 함.



#interCal:::interaction계산기 #######################################
#데이터 생성  2x2
# 2by 2 생성
grp22<-function(xx,...){
  gr <-rbind(xx,...)
  colnames(gr)=c("M","SD","N")
  rownames(gr)=c("grp1_ap","grp2_aq","grp3_bp","grp4_bq")
  gr<-gr
  gr
}

# 데이터가 많은 경우 생성
grpdata<-function(xx,...){
  gr <-rbind(xx,...)
  colnames(gr)=c("M","SD","N")
  rownames(gr)=paste0(rep("grp_",nrow(gr)),1:nrow(gr))
  gr<-gr
}

#gr<-grp22(c(24.75,7.44831,12),c(23.5833,7.85,12),c(23.3333,7.36495,12),c(12.9167,6.40253,12))
#gr %>% class
#grpdata(c(24.89,2.37,9),c(21.44,2.79,9),c(14.56,2.46,9),c(12,2.69,9),c(1,2,3),c(34,3.2,1))

# MSDN으로 interaction 계산######
#데이터 생성  2x2
# 2by 2 생성
grp22<-function(xx,...){
  gr <-rbind(xx,...)
  colnames(gr)=c("M","SD","N")
  rownames(gr)=c("grp1_ap","grp2_aq","grp3_bp","grp4_bq")
  gr<-gr
  gr
}
#interCal #####
interCal <- function(MSDN, text="",sel=1){
  library(dplyr, warn.conflicts = FALSE)
  # Suppress summarise info
  options(dplyr.summarise.inform = FALSE)

  if(!is.character(MSDN)){if(!is.data.frame(MSDN)){
    if(!is.matrix(MSDN)){return("정확한 데이터인 dataframe 이나 (2X2)matrix를 입력해주세요.grp22를 이용하여 데이터를 생성하세요")}}
  }else{gdata <- MSDN}

 #데이터 입력
  gdata <- MSDN

  #구분변수 입력 파트
  if(is.character(text)){
    gl <-strsplit(text, split= "\\,") %>% unlist()
    grp1.ap <-paste("m1 =", gl[1], gl[3])
    grp2.aq <-paste("m2 =", gl[1], gl[4])
    grp3.bp <-paste("m3 =", gl[2], gl[3])
    grp4.bq <-paste("m4 =", gl[2], gl[4])
    ggt <- rbind(grp1.ap,grp2.aq ,grp3.bp,grp4.bq)
    gdata.total <-cbind(ggt,gdata)}
  else(return("문자변수로 `,로 구분해서입력하세요. 문자변수가 없으면 sel=1 or sel=2로 입력하세요"))

   library(knitr)
  gdata.total<-gdata.total %>% kable("pandoc", caption = "입력된 기술통계 데이터 M,SD,N ")


#자유도 및 합동분산 계산
  M<-gdata[,1]
  SD<-gdata[,2]
  N<-gdata[,3]

  df1<-N-1
  df1
  var=SD^2
  var
  spsq <-sum(df1*var)/sum((df1))
  #gdata$con<-c("m1","m2","m3","m4")
  #합동분산 pooled variance


  #contrast setting
  # ap aq bp bq
  contr1= cbind("inter:ap/aq/bp/bq:m1-m2-m3+m4=0"= c(1,-1,-1,1), #interaction
                "me_of_ab:(m1+m2)-(m3+m4)=0"    =c(1,1,-1,-1), # (m1+m2)-(m3+m4)=0
                "sme.of.ab_at.p:m1-m3=0"        =c(1,0,-1,0),  # m1-m3=0
                "sme.of.ab_at.q:m2-m4=0"        =c(0,1,0,-1),  # m2-m4=0

                "me_of_pq:(m1+m3)-(m2+m4)=0"=c(1,-1,1,-1),   #(m1+m3)-(m2-m4) =0
                "sme.of.pq_at.a:m1-m2=0"    =c(1,-1,0,0),    # m1-m2=0
                "sme.of.pq_at.b:m3-m4=0"    =c(0,0,1,-1))    # m3-m4=0
  rownames(contr1)=c("m1","m2","m3","m4")

  contr2= cbind("inter:ap/aq/bp/bq:m1-m2-m3+m4=0"= c(1,-1,-1,1),
                #"me_of_ab:(m1+m2)-(m3+m4)=0"=c(1,1,-1,-1), # (m1+m2)-(m3+m4)=0
                "sme.of.ab_at.p:m1-m3=0"     =c(1,0,-1,0),  # m1-m3=0
                "sme.of.ab_at.q:m2-m4=0"     =c(0,1,0,-1),  # m2-m4=0

                #"me_of_pq:(m1+m3)-(m2+m4)=0"=c(1,-1,1,-1),   #(m1+m3)-(m2-m4) =0
                "sme.of.pq_at.a:m1-m2=0"   =c(1,-1,0,0),    # m1-m2=0
                "sme.of.pq_at.b:m3-m4=0"   =c(0,0,1,-1))    # m3-m4=0
  rownames(contr2)=c("m1","m2","m3","m4")


  if(sel==1){contr<- contr1} #total 출력
  else{if(sel==2){contr <-contr2}   #simple main effect만 출력
    else{return("selcon = 1과 selcon = 2중에서 선택하세요")}}



  #행렬 계수*평균
  cim<- contr * M
  #cim

  #평균합
  cim_sum <-apply(cim, 2, sum)
  #cim_sum

  #등분산인 경우
  ci2sd2<- contr^2*spsq  #계수제곱* pooled variance
  ci2sd2_N<- ci2sd2/N    #N을 나눈 것
  SE <- sqrt(apply(ci2sd2_N,2,sum))
  #SE  #standard Erroe

  #등분산이 아닌 경우
  #{ci2sd2_notvar <- (contr1^2 *gdata[,2]^2)/N
  #SE.notvar<- sqrt(apply(ci2sd2_notvar,2,sum))
  #SE.notvar}


  df <- sum(N)-nrow(gdata)

  t <- cim_sum/SE
  cri.t <- abs(qt(0.05/2,df))
  p.val<- 2*(1-pt(abs(t),df))
  #p.val
  LowerCI <- cim_sum - SE * cri.t
  UpperCI <- cim_sum + SE * cri.t
  star <- ifelse(p.val<0.001,"***",ifelse(p.val<0.01,"**",ifelse(p.val<0.05,"*","NotSig")))
  options(scipen = 100)
  table <- cbind(cim_sum, SE, t, cri.t,df,LowerCI,UpperCI,p.val) %>% round(3)
  #pooled.var=spsq,
  table1<- cbind(table, star)
  tab <- table1 %>% as.data.frame


  res<-list(gdata.total,"대비검정"=tab,"대비계수"=t(contr), "합동분산"=spsq)
  res#결과

}

#사용방법

#실행하기gl

#a1 =c(24.75,   7.45 ,  12)
     #a2 =c(23.58,   7.86,   12)
     #a3=c( 23.33,   7.36 ,  12)
     #a4=c(12.92 ,  6.40,   12)

     #gr<-grp22(a1,a2,a3,a4)


     #interCal(grp22(a1,a2,a3,a4))
   #interCal(grp22(c(24.75,   7.45 ,  12),
     #               c(23.58,   7.86,   12),
     #               c( 23.33,   7.36 ,  12),
     #               c(12.92 ,  6.40,   12)  ))
     #interCal(grp22(a1,a2,a3,a4),"aa,bb,pp,qq",2)
     #interCal(grp22(a1,a2,a3,a4),1)
     #interCal(grp22(a1,a2,a3,a4),1,222)


     #interCal(gr) #maine
     # interCal(gr,sel=2)
     # interCal(gr,"")
     #interCal(gr,2) #sme
     #interCal(gr,sel=1) #maine
# interCal(gr,3) #오류 메시지 출력
# interCal(c(12,3,2),c(1,2,3))  # 오류메시지 확인
# interCal(c("al dadas",1))  #오류 메지시 확인


#################################

#t test를 SPSS양식으로 ####
ttestspss <- function(formula,data){
  library(dplyr)
  library(knitr)
  library(car)
  levent<-leveneTest (formula, data ,center=mean)
  ff<-levent$`F value`[1]
  pp<-levent$`Pr(>F)`[1]


  #student ttest
  mtc_vs_ttest <-t.test(formula, data, var.equal=T)

  t_val<- mtc_vs_ttest$statistic #t value
  t_df<-mtc_vs_ttest$parameter #df
  t_pval<- mtc_vs_ttest$p.value
  t_mean<- mtc_vs_ttest$estimate #mean
  t_SE<-mtc_vs_ttest$stderr
  t_CI<- t(mtc_vs_ttest$conf.int) #CI

  #Whechl t test
  ttest <-t.test(formula, data, var.equal=F)

  t_val1<- ttest$statistic #t value
  t_df1<-ttest$parameter #df
  t_pval1<- ttest$p.value
  t_mean1<- ttest$estimate #mean
  t_SE1<-ttest$stderr
  t_CI1<- t(ttest$conf.int) #CI

  #Levene적용
  result1<-cbind(Levene_F=ff,Leven_p=pp,t=t_val,df=t_df,p=t_pval,SE=t_SE,CI=t_CI) %>% round(3)
  result2<- cbind(Levene_F=ff,Leven_p=pp,t=t_val1,df=t_df1,p=t_pval1,SE=t_SE1,CI=t_CI1) %>% round(3)


  #result1<-cbind(t=t_val,df=t_df,p=t_pval,SE=t_SE,CI=t_CI) %>% round(4)
  #result2<- cbind(t=t_val1,df=t_df1,p=t_pval1,SE=t_SE1,CI=t_CI1) %>% round(4)
  result<-rbind(result1, result2)
  #colnames(result)=c("t","df","p","Std.Error","95%CI lower","95%CI upper" )
  colnames(result)=c("Levene's F","Sig","t","df","p","Std.Error","95%CI lower","95%CI upper" )
  rownames(result)=c("Equal var. assumed","Equal var. Not assumed")
  result %>% kable("pandoc", caption = "Independent t-test of SPSS type table By park Joonghee(2020)")
}
#독립변수는 factor로 할 것
# ttestspss(mpg~factor(am),mtcars)
#ttestspss(mpg~am,mtcars)  #error


#유의성 별 만들기 #####
star_make<-function(data,p.value){
  ndata<- data %>% as.data.frame() %>%
    mutate(stars=ifelse(p.value < 0.001, "***",
                        ifelse(p.value < 0.01, "**",
                               ifelse(p.value < 0.05, "*", ""))))
  ndata
}

#사용
#star_make(mmm,v2)

#v1=c(1,2,3,4)
#v2=c(0.1,0.04,0.003,0.00001)
#mmm<-cbind(v1,v2)
#mmm
#star_make(mmm,v2)





#카이제곱 검정::오류메시지 해결하는 함수 제작 ####
#The code for Monte Carlo simulation is a C translation of the Fortran algorithm of Patefield (1981).
#Hope, A. C. A. (1968). A simplified Monte Carlo significance test procedure. Journal of the Royal Statistical Society Series B, 30, 582–598. http://www.jstor.org/stable/2984263.
#Agresti, A. (2007). An Introduction to Categorical Data Analysis, 2nd ed. New York: John Wiley & Sons. Page 38.
advanced.chisq.test<-function(cri, option=TRUE){
  cc <- cri[rowSums(cri)>0,]
  if(option == TRUE){
    chisq.test(cc, simulate.p.value = TRUE)}
  else{chisq.test(cc)}
}

#실행
#advanced.chisq.test(online_freq)














#sobel test my function 박중희 #####
#sobel test my function 박중희 #####
sobel.test <- function(x,m,y){
  options(scipen=10)
  library(multilevel)
  sob.model<- sobel(pred=x,med=m,out=y)
  p_value<-  2*(1-pnorm(abs(sob.model$z.value)))
  library(bda)
  bda <-bda::mediation.test(m,x,y)
  z.value<-bda$Sobel[1]
  p.value<-bda$Sobel[2]
  sig <-ifelse(p_value < 0.05,paste0("간접효과(indirect effect=ab)는 통계적으로 유의미하다,"," z = ",
                                     round( z.value,2),", p = ",
                                     round(p.value,4),"."),
               "간접효과(indirect effect)는통계적으로 유의미하지 않다")
  res<-list(c(Sobel_z=sob.model$z.value,p.value=p_value,ind.effect=sob.model$Indirect.Effect),sig, bda)
  res
}

#sobel test my function 박중희 #####
BKsobel.test <- function(x,m,y){
  options(scipen=10)
  library(multilevel)
  sob.model<- sobel(pred=x,med=m,out=y)
  p_value<-  2*(1-pnorm(abs(sob.model$z.value)))
  library(bda)
  bda <-bda::mediation.test(m,x,y)
  z.value<-bda$Sobel[1]
  p.value<-bda$Sobel[2]
  #step1 total effect
  c.lm <- lm(y ~ x)
  a.lm <- lm(m ~ x)
  cp.b.lm <-lm(y~x+m)

  a<- coef(summary(a.lm))[2,1]
  SEa<- coef(summary(a.lm))[2,2]
  a_p.value <- coef(summary(a.lm))[2,4]
  b <- coef(summary(cp.b.lm))[3,1]
  SEb <-coef(summary(cp.b.lm))[3,2]
  b_p.value <- coef(summary(cp.b.lm))[3,4]
  cp<-coef(summary(cp.b.lm))[2,1]
  SEcp<-coef(summary(cp.b.lm))[2,2]
  cp_p.value<-coef(summary(cp.b.lm))[2,4]

  c<-coef(summary(c.lm))[2,1]
  SEc<-coef(summary(c.lm))[2,2]
  c_p.value<-coef(summary(c.lm))[2,4]

  Direct_effect.cp<-cp
  Total_effect.c<-c
  Indirect_effect.ab=a*b

  #BK step causal step
  coeff.c<-cbind(c,SEc,c_p.value ) #1
  coeff.a<-cbind(a,SEa,a_p.value) #2
  coeff.b<-cbind(b,SEb, b_p.value)  #3
  coeff.cp<-cbind(cp,SEcp,cp_p.value ) #3
  coeff=rbind(coeff.c,coeff.a,coeff.b,coeff.cp)
  name=c("BK:step1_X->Y(c)","BK:step2_X->M(a)" ,"BK:step3_M->Y(b)", "BK:step3_X->Y(c')")
  cname=c("B","SE","p.value","stars")

  #coeff<-data.frame(coeff)
  library(dplyr)
  coeff <- coeff %>% as.data.frame() %>%
    mutate(stars=ifelse(p.value < 0.001, "***",
                        ifelse(p.value < 0.01, "**",
                               ifelse(p.value < 0.05, "*", ""))))
  rownames(coeff)=name
  colnames(coeff)=cname

  effect<- cbind(Total_effect.c, Direct_effect.cp, Indirect_effect.ab)

  sig <-ifelse(p_value < 0.05,paste0("간접효과(indirect effect=ab)는 통계적으로 유의미하다,"," z = ",
                                     round( z.value,2),", p = ",
                                     round(p.value,4),"."),
               "간접효과(indirect effect)는통계적으로 유의미하지 않다")
  res<-list(coeff,effect,
            c(Sobel_z=sob.model$z.value,p.value=p_value,ind.effect=sob.model$Indirect.Effect),sig, bda)
  res
}

#각각의 변수를 넣어야 하는 함수X,M,Y
sobel.test2(tm7$facultysupport,tm7$achievementmotivation, tm7$publication)



#mtcars

BKsobel.test(mtcars$disp, mtcars$wt, mtcars$mpg)
sobel.test(mtcars$disp, mtcars$wt, mtcars$mpg)

#사용법:: 각의 변수를 넣어야 하는 함수
# sobel.test(tm7$achievementmotivation,tm7$facultysupport, tm7$publication)

#직접 값을 넣는 함수
SOBEL<-function(a,sa,b,sb){
  options("scipen"=10)
  AB=a*b
  A=a^2
  SA=Sa^2
  B=b^2
  SB=Sb^2

  rootSSC= sqrt(A*SB+B*SA)
  z=AB/rootSSC  #z value
  pval <- 2*(1-pnorm(abs(z))) #p value
  sig <-ifelse(pval < 0.05,"간접효과(indirect effect = ab)는 통계적으로 유의미하다","간접효과(indirect effect)는통계적으로 유의미하지 않다")
  list(c(z.value= z,p.value= pval), sig)
}


#a= 0.6694448
#Sa=0.09998841
#b=.8506020
#Sb=0.1437364

SOBEL(a,sa,b,sb)




#다차원척도법 plot####
#함수그리기 자동화

mdsplot<- function(point){
  #par(mar=c(4,4,4,4))
  plot(point[,1],point[,2], type = "p",
       xlim=c(-3,3),ylim=c(-3,3),pch=21,bg=c(1:ncol(point)),cex=2,
       main="다차원척도 2차원 plot by JH Park")
  grid()
  text(point+.2,cex=1)
  abline(v=0,lty=2,col="red")
  abline(h=0,lty=2,col="red")
}


mdsplot1<- function(point,variable,main_Title){
  #par(mar=c(4,4,4,4))
  plot(point[,1],point[,2], type = "p",
       xlim=c(-3,3),ylim=c(-3,3),pch=21,bg=c(1:ncol(point)),cex=2,
       main=main_Title)
  #     xlab="편리성",ylab="mobile")
  grid()
  text(point+.2,rownames(variable),cex=1)
  abline(v=0,lty=2,col="red")
  abline(h=0,lty=2,col="red")
}


mdsplot_edit<- function(point,variable,main_Title,xname,yname){
  par(mar=c(5,5,5,5))
  plot(point[,1],point[,2], type = "p",
       xlim=c(-3,3),ylim=c(-3,3),pch=21,bg=c(1:ncol(point)),cex=2,
       main=main_Title,
       xlab=xname,ylab=yname)
  #     xlab="편리성",ylab="mobile")
  grid()
  text(point+.2,rownames(variable),cex=1)
  abline(v=0,lty=2,col="red")
  abline(h=0,lty=2,col="red")
}



#MDS에 대한 데이터 isoMDS- data가 fulldata인 경우
#0을 포함한 자료인 경우
mdsDisMat<-function(data){
  library(dplyr)
  library(psych)
  library(lavaan)
  library(MASS)
#  data<-mds_onl
  des.data <- data %>% as.data.frame() %>% describe()
  mat <- des.data[,3] %>% lav_matrix_upper2full() #행렬
  dist<- mat %>% scale()%>% dist() #표준화 및 거리 함수
  iso<- dist %>% isoMDS(k=2) #MDS분석
  isoplot <- iso$points #MDS결과
  colnames(isoplot )=c("x","y")
  #shepard plot
  dist_sh <- Shepard(mat[lower.tri(mat)], iso$points)
  cdist_sh <-cbind(dist_sh$x, dist_sh$y, dist_sh$yf)
  colnames(cdist_sh)=c("x","y","yf")
  result<-list(Matrix=mat,Distance=dist,isoMDS_Result=iso,Plot_data=isoplot,Shepard_data=cdist_sh)
  result
}
#실행
#mdsDisMat(mds_onl)


#6개에 맞추어 둠 행렬제작 ####
#MDS데이터에 0이 안들어간 경우 만들어주기

#6개에 맞추어 둠 행렬제작 ####
#MDS데이터에 0이 안들어간 경우 만들어주기
mdsDataMaker <- function(data){
  library(psych)
  library(lavaan)
  library(dplyr)

  mds_data <- data %>% as.data.frame()
  mds0 =rep(0, nrow(mds_data)) %>% as.data.frame()

  dms_raw <-cbind(mds0,mds_data[,1:5],mds0,mds_data[,6:9],mds0,mds_data[,10:12],mds0,mds_data[,13:14],mds0,mds_data[,15],mds0)

  des.data <- describe(dms_raw)
  des.data <- des.data %>% as.data.frame()
  result.matrix <-  lav_matrix_upper2full(des.data[,3])
  colnames(result.matrix)=c("V1","V2","V3","V4","V5","V6")
  rownames(result.matrix)=c("V1","V2","V3","V4","V5","V6")
  #result.matrix
  result.matrix
}


#행렬의 열과 행의 이름 넣기####
matrix_name<- function(Matrix, name){
  colnames(Matrix)=name
  rownames(Matrix)=name
  Matrix
}




#SDS:scale -dist-shepard  plot####
#행렬을 입력하여 사용


mds_SDS_plot <-function(mat_data){
  library(dplyr)
  library(psych)
  library(lavaan)
  library(MASS)

  #des.data <- data %>% describe(data)
  #mat <- des.data[,3] #행렬
  mat<- mat_data %>% scale() #표준화
  dist<-mat %>% dist() #거리함수
  iso <-isoMDS(dist, k=2) #MDS분석
  isoplot <- iso$points #MDS결과plot자료
  stress_value <-iso$stress

  cmd<-cmdscale(dist,k=2,eig = T) #metric MDS
  cmdplot<-cmd$points
  cmd_Gof<-cmd$GOF

  #shepard plot
  dist_sh <- Shepard(dist, iso$points) #shepard plot 데이터
  cdist_sh=cbind(dist_sh$x, dist_sh$y, dist_sh$yf)
  colnames(cdist_sh)=c("원거리", "FitDist", "FitDATA")

  res=list(dist, cmdscale= cmdplot,cmdGOF = cmd_Gof, isoMDS=iso, Shepard=cdist_sh)
  res
}

#행렬만 가지고 분석
#mds_SDS_plot(on_mat)

#전체분석 도구 : :image를 순차적으로 분리함 #####
mds_SDS_plot_step <-function(mat_data,opt=FALSE){
  library(dplyr)
  library(psych)
  library(lavaan)
  library(MASS)

  #des.data <- data %>% describe(data)
  #mat <- des.data[,3] #행렬


  mat<- mat_data %>% scale() #표준화
  dist<-mat %>% dist() #거리함수

  #isoMDS data
  iso <-isoMDS(dist, k=2) #MDS분석
  isoplot <- iso$points #MDS결과plot자료
  stress_value <-iso$stress

  #cmdscale MDS data
  cmd<-cmdscale(dist,k=2,eig = T) #metric MDS
  cmdplot<-cmd$points
  cmd_Gof<-cmd$GOF

  #shepard plot data
  dist_sh <- Shepard(dist, iso$points) #shepard plot 데이터
  cdist_sh=cbind(dist_sh$x, dist_sh$y, dist_sh$yf)
  colnames(cdist_sh)=c("원거리", "FitDist", "FitDATA")


  #par(mfrow=c(2,2))

  #cmdscale plot
  #바탕만들기

  x1<-cmdplot[,1]
  y1<-cmdplot[,2]
  limx1<-c(-max(abs(x1)),max(abs(x1)))
  limy1<-c(-max(abs(y1)),max(abs(y1)))

  plot(cmdplot[,1],cmdplot[,2],type = "n",#pch=21,bg=c(1:6),
       xlim=limx1, ylim = limy1,
       xlab="",ylab="",axes=FALSE)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], border = FALSE, col="gray80")
  #par("usr")
  #abline(h=-3:3,col="white")
  #abline(v=-3:3, col="white")
  par(new=TRUE)

  plot(cmdplot[,1],cmdplot[,2], type = "p",
       xlim=limx1, ylim = limy1,
       pch=21,bg="white",#c(1:ncol(cmdplot)),
       cex=2,
       main="Metric cmdscale  Perceptual Map ",
       xlab= paste("Goodnes Of Fit =",round(cmd_Gof[1]*100,2),collapse=""),ylab="y variable ")
  grid(col = "white")
  text(cmdplot +c(.1,-.1),
       rownames(cmdplot),cex=1,col = "black")
  abline(v=0,lty=2,col="gray40")
  abline(h=0,lty=2,col="gray40")
  mtext("Yong & Householder(1938)")



  #isoMDS plot
  #colnames(isoplot)=c("X-Rename Variable set","Y-Rename Variable set")
  #바탕을 회색으로 만들기
  x2<-isoplot[,1]
  y2<-isoplot[,2]
  limx2<-c(-max(abs(x2)),max(abs(x2)))
  limy2<-c(-max(abs(y2)),max(abs(y2)))

  plot(isoplot[,1],isoplot[,2],type = "n",pch=21,bg=c(1:6),
       xlim=limx2, ylim = limy2,
       xlab="",ylab="",axes=FALSE)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], border = FALSE, col="gray80")
  #par("usr")
  #abline(h=-3:3,col="white")
  #abline(v=-3:3, col="white")
  par(new=TRUE)

  plot(isoplot[,1],isoplot[,2], type = "p",
       xlim=limx2, ylim = limy2,
       pch=21,bg="white",#c(1:ncol(isoplot)),
       cex=2,
       main="Nonmetric isoMDS Perceptual Map ",
       xlab= paste("stress =",round(stress_value,5),collapse="") ,ylab="y variable")
  grid(col = "white")
  text(isoplot+c(.1,-.1),
       rownames(cmdplot),cex=1,col = "black")
  abline(v=0,lty=2,col="gray40")
  abline(h=0,lty=2,col="gray40")
  mtext("Shepard(1962a,b), Kruskal(19641,b)")

  #shepard plot
  plot(cdist_sh[,1],cdist_sh[,3], pch=".",
       xlim=range(cdist_sh[,1]), ylim = range(cdist_sh[,3]),
       xlab="Dissimilarity", ylab="Distance",main="MDS_Guttman(1968) Shepard Diagram  ")#,
  #xlim=range(cdist_sh[,1]), ylim = range(cdist_sh[,1]))
  grid()
  lines(cdist_sh[,1],cdist_sh[,3], type="p",pch=1,col="red")
  lines(cdist_sh[,1],cdist_sh[,3], type="S", col="red")
  abline(lm(cdist_sh[,3]~cdist_sh[,1]),lty=2,col="gray70")#regression line
  #abline(0,1,lty=2,col="gray70")
  mtext("45도에 가까운 선형이 나타나면 적합")

  #image plot
  plot(cdist_sh[,2],cdist_sh[,3],xlim=range(cdist_sh[,2]), ylim = range(cdist_sh[,3]),
       xlab="Fit.Dissimilarity", ylab="Distance",main="MDS_Residual plot(image plot)")
  grid()
  lines(cdist_sh[,2],cdist_sh[,3], type="p",pch=21, bg=c(1:10),cex=1.3)
  lines(cdist_sh[,2],cdist_sh[,3], type="S",lty=2,col="gray70")

  abline(lm(cdist_sh[,3]~cdist_sh[,2]),lty=2,col="red")#regression line
  #par(mfrow=c(1,1))


  res=list(dist, cmdscale=cmdplot,cmd_Gof, isoMDS=iso, Shepard=cdist_sh)
  if(opt==TRUE)  {print(res)}else{return("결과를 보려면 opt=TRUE or T를 입력하세요 ")}
}
#mds_SDS_plot_step(on_mat,T)

#mds_SDS_plot_step(on_mat,F)

#연구결과 제출시 활용할 데이터 ####
mds_SDS_plot2by2 <-function(mat_data,opt=FALSE){
  library(dplyr)
  library(psych)
  library(lavaan)
  library(MASS)

  #des.data <- data %>% describe(data)
  #mat <- des.data[,3] #행렬


  mat<- mat_data %>% scale() #표준화
  dist<-mat %>% dist() #거리함수

  #isoMDS data
  iso <-isoMDS(dist, k=2, trace = FALSE) #MDS분석
  isoplot <- iso$points #MDS결과plot자료
  stress_value <-iso$stress

  #cmdscale MDS data
  cmd<-cmdscale(dist,k=2,eig = T) #metric MDS
  cmdplot<-cmd$points
  cmd_Gof<-cmd$GOF

  #shepard plot data
  dist_sh <- Shepard(dist, iso$points) #shepard plot 데이터
  cdist_sh=cbind(dist_sh$x, dist_sh$y, dist_sh$yf)
  colnames(cdist_sh)=c("원거리", "FitDist", "FitDATA")


  par(mfrow=c(2,2))

  #cmdscale plot
  #바탕만들기

  x1<-cmdplot[,1]
  y1<-cmdplot[,2]
  limx1<-c(-max(abs(x1)),max(abs(x1)))
  limy1<-c(-max(abs(y1)),max(abs(y1)))

  plot(cmdplot[,1],cmdplot[,2],type = "n",#pch=21,bg=c(1:6),
       xlim=limx1, ylim = limy1,
       xlab="",ylab="",axes=FALSE)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], border = FALSE, col="gray80")
  #par("usr")
  #abline(h=-3:3,col="white")
  #abline(v=-3:3, col="white")
  par(new=TRUE)

  plot(cmdplot[,1],cmdplot[,2], type = "p",
       xlim=limx1, ylim = limy1,
       pch=21,bg="white",#c(1:ncol(cmdplot)),
       cex=2,
       main="Metric cmdscale  Perceptual Map ",
       xlab= paste("Goodnes Of Fit =",round(cmd_Gof[1]*100,2),collapse=""),ylab="y variable ")
  grid(col = "white")
  text(cmdplot +c(.1,-.1),
       rownames(cmdplot),cex=1,col = "black")
  abline(v=0,lty=2,col="gray40")
  abline(h=0,lty=2,col="gray40")
  mtext("Yong & Householder(1938)")



  #isoMDS plot
  #colnames(isoplot)=c("X-Rename Variable set","Y-Rename Variable set")
  #바탕을 회색으로 만들기
  x2<-isoplot[,1]
  y2<-isoplot[,2]
  limx2<-c(-max(abs(x2)),max(abs(x2)))
  limy2<-c(-max(abs(y2)),max(abs(y2)))

  plot(isoplot[,1],isoplot[,2],type = "n",pch=21,bg=c(1:6),
       xlim=limx2, ylim = limy2,
       xlab="",ylab="",axes=FALSE)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], border = FALSE, col="gray80")
  #par("usr")
  #abline(h=-3:3,col="white")
  #abline(v=-3:3, col="white")
  par(new=TRUE)

  plot(isoplot[,1],isoplot[,2], type = "p",
       xlim=limx2, ylim = limy2,
       pch=21,bg="white",#c(1:ncol(isoplot)),
       cex=2,
       main="Nonmetric isoMDS Perceptual Map ",
       xlab= paste("stress =",round(stress_value,5),collapse="") ,ylab="y variable")
  grid(col = "white")
  text(isoplot+c(.1,-.1),
       rownames(cmdplot),cex=1,col = "black")
  abline(v=0,lty=2,col="gray40")
  abline(h=0,lty=2,col="gray40")
  mtext("Shepard(1962a,b), Kruskal(19641,b)")

  #shepard plot
  plot(cdist_sh[,1],cdist_sh[,3], pch=".",
       xlim=range(cdist_sh[,1]), ylim = range(cdist_sh[,3]),
       xlab="Dissimilarity", ylab="Distance",main="MDS_Guttman(1968) Shepard Diagram  ")#,
  #xlim=range(cdist_sh[,1]), ylim = range(cdist_sh[,1]))
  grid()
  lines(cdist_sh[,1],cdist_sh[,3], type="p",pch=1,col="red")
  lines(cdist_sh[,1],cdist_sh[,3], type="S", col="red")
  abline(lm(cdist_sh[,3]~cdist_sh[,1]),lty=2,col="gray70")#regression line
  #abline(0,1,lty=2,col="gray70")
  mtext("45도에 가까운 선형이 나타나면 적합")

  #image plot
  plot(cdist_sh[,2],cdist_sh[,3],xlim=range(cdist_sh[,2]), ylim = range(cdist_sh[,3]),
       xlab="Fit.Dissimilarity", ylab="Distance",main="MDS_#Residual plot(image plot)")
  grid()
  lines(cdist_sh[,2],cdist_sh[,3], type="p",pch=21, bg=c(1:10),cex=1.3)
  lines(cdist_sh[,2],cdist_sh[,3], type="S",lty=2,col="gray70")

  abline(lm(cdist_sh[,3]~cdist_sh[,2]),lty=2,col="red")#regression line
  par(mfrow=c(1,1))


  res<- list(dist, cmdscale=cmdplot,cmd_Gof, isoMDS=iso, Shepard=cdist_sh)
  if(opt==TRUE)print(res)
}

#Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis (Gower, 1966).

#Usage
#mds_SDS_plot2by2(on_mat, opt=TRUE)
#mtext("전체데이터")
#mds_SDS_plot2by2(gachon,TRUE)
#mtext("가천대 대학원생 ")
#mds_SDS_plot2by2(yonsei,F)
#mtext("연세대 대학원생")
#mds_SDS_plot2by2(academy)
#mtext("평생교육원 수강생")


#cmdcale and isoMDS plot ####
cmdscale_plot <- function(mds_data, title="",line=1:6,...) {
  if(is.numeric(line)==FALSE){
    stopifnot("입력란이 비었습니다. 숫자로다시입력하세요 ")
  }
  library(shapes)

  cmdplot_d <- mds_data
  cmdplot   <- cmdplot_d$cmdscale
  cmd_Gof   <- cmdplot_d$cmdGOF

  mds_data_d <- mds_data
  mds_data1  <- mds_data_d$cmdscale
  mds_stress <- mds_data_d$cmdGOF

  par(mfrow=c(1,2))
  #order
  plot(cmdplot[,1],cmdplot[,2], type = "b",
       xlim=c(-2.5,2.5),
       ylim=c(-2.5,2.5),lty=2,col=1,
       pch=21,bg=c(1:nrow(cmdplot)),cex=0.9,
       main= title ,
       xlab= paste("GOF =",round(cmd_Gof[1],3),collapse=""),ylab="y variable - Must Rename"
  )
  grid()
  text(cmdplot+c(-.2, .3),rownames(cmdplot),cex=.9, col = c(1,1,1,2,2,2))
  abline(v=0,lty=2,col="gray50")
  abline(h=0,lty=2,col="gray60")


  #shape 1,4,2,5,3,6
  plotshapes(mds_data1, joinline=c(line,...),color="black",symbol = 1)#mteric MDS 형상
  abline(v=0,h=0,col="gray50")
  points(mds_data1,pch=21,bg=c(3,3,3,6,6,6),cex=1)
  grid()
  text(mds_data1+c(-.2, .3),rownames(mds_data1),cex=0.9, col = c(rep("blue",3), rep("red",3)))
  title(paste("NonMetric MDS Shape of", title),paste("GOF(%) =",round(mds_stress*100,2)))
  par(mfrow=c(1,1))

}
#실행
# cmdscale_plot(mds_online)
# cmdscale_plot(mds_online,"Total MDS")
# cmdscale_plot(mds_online,"Total MDS", 1,6,3,2,4,5,1)
# cmdscale_plot(mds_gachon,"gachon student",1,6,3,2,4,5,1)
# cmdscale_plot(mds_yonsei,"yonsei student",1,6,3,2,4,5,1)
# cmdscale_plot(mds_academy,"Lecture",1,6,3,2,4,5,1)



#compare metric and nonmetric MDS
#계량형과 비계량형을 비교
MDS_plot <- function(mds_data,title="",line=c(1:6,1),...) {
  if(is.character(title)==FALSE){
    stopifnot("input title ")
  }
  library(shape)

  cmdplot_d <- mds_data
  cmdplot   <- cmdplot_d$cmdscale
  cmd_Gof   <- cmdplot_d$cmdGOF

  mds_data_d <- mds_data
  mds_data1  <- mds_data_d$isoMDS$points
  mds_stress <- mds_data_d$isoMDS$stress

  #par(bg = "gray90")
  par(mfrow=c(1,2)) #1x2 figuare
  #Metric MDS shape
  plotshapes(cmdplot, joinline=c(line,...),color="red",symbol = 1)#mteric MDS 형상
  abline(v=0,h=0)
  points(cmdplot,pch=21,bg=c(1,1,1,6,6,6),cex=1)
  grid()
  text(cmdplot+c(-.2, .3),rownames(cmdplot),cex=0.9, col = c(rep("blue",3), rep("red",3)))
  title(paste("Metric MDS of", title),paste("GOF(%) =",round(cmd_Gof*100),2))

  #Nonmetric shape coordinated
  plotshapes(mds_data1, joinline=c(line,...),color="red",symbol = 1)#mteric MDS 형상
  abline(v=0,h=0)
  points(mds_data1,pch=21,bg=c(1,1,1,6,6,6),cex=1)
  grid()
  text(mds_data1+c(-.2, .3),rownames(mds_data1),cex=0.9, col = c(rep("blue",3), rep("red",3)))
  title(paste("NonMetric MDS  of", title),paste("stress =",round(mds_stress/100,4)))

  par(mfrow=c(1,1))


}

#실행
#par(bg = "white")
# MDS_plot(mds_online)
# MDS_plot(mds_online,"total")
# MDS_plot(mds_online,"total", 1,6,2,3,4,5,1)
# MDS_plot(mds_gachon,"Gachon",1,5,4,2,3,6,1)
# MDS_plot(mds_yonsei,"yonsei",1,5,4,3,6,2,1)
# MDS_plot(mds_academy,"lecture",1,6,3,2,4,5,1)



#각각 하나만을 그려야 하는 경우
cmd_single_plot <- function(mds_data,title="",line=c(1:6,1),...) {
  if(is.character(title)==FALSE){
    stopifnot("input title ")
  }
  library(shape)

  cmdplot_d <- mds_data
  cmdplot   <- cmdplot_d$cmdscale
  cmd_Gof   <- cmdplot_d$cmdGOF

  mds_data_d <- mds_data
  mds_data1  <- mds_data_d$isoMDS$points
  mds_stress <- mds_data_d$isoMDS$stress

  #par(bg = "gray90")
  #par(mfrow=c(1,2)) #1x2 figuare
  #Metric MDS shape
  plotshapes(cmdplot, joinline=c(line,...),color="red",symbol = 1)#mteric MDS 형상
  abline(v=0,h=0)
  points(cmdplot,pch=21,bg=c(1,1,1,6,6,6),cex=1)
  grid()
  text(cmdplot+c(-.2, .3),rownames(cmdplot),cex=0.9, col = c(rep("blue",3), rep("red",3)))
  title(paste("Metric MDS of", title),paste("GOF(%) =",round(cmd_Gof*100),2))

  par(mfrow=c(1,1))

}

iso_single__plot <- function(mds_data,title="",line=c(1:6,1),...) {
  if(is.character(title)==FALSE){
    stopifnot("input title ")
  }
  library(shape)

  cmdplot_d <- mds_data
  cmdplot   <- cmdplot_d$cmdscale
  cmd_Gof   <- cmdplot_d$cmdGOF

  mds_data_d <- mds_data
  mds_data1  <- mds_data_d$isoMDS$points
  mds_stress <- mds_data_d$isoMDS$stress



  #Nonmetric shape coordinated
  plotshapes(mds_data1, joinline=c(line,...),color="red",symbol = 1)#mteric MDS 형상
  abline(v=0,h=0)
  points(mds_data1,pch=21,bg=c(1,1,1,6,6,6),cex=1)
  grid()
  text(mds_data1+c(-.2, .3),rownames(mds_data1),cex=0.9, col = c(rep("blue",3), rep("red",3)))
  title(paste("NonMetric MDS  of", title),paste("stress =",round(mds_stress/100,4)))

  par(mfrow=c(1,1))
}

#형상비교를 위한 plot함수
plotShapes<- function (A, B = 0, joinline = c(1, 1), orthproj = c(1, 2), color = 1,
                       symbol = 1) {
  CHECKOK <- TRUE
  if (is.array(A) == FALSE) {
    if (is.matrix(A) == FALSE) {
      cat("Error !! argument should be an array or matrix \n")
      CHECKOK <- FALSE
    }
  }
  if (CHECKOK) {
    k <- dim(A)[1]
    m <- dim(A)[2]
    kk <- k
    if (k >= 15) {
      kk <- 1
    }

    par(pty = "s")
    if (length(c(B)) != 1) {
      par(mfrow = c(1, 2))
    }
    if (length(dim(A)) == 3) {
      A <- A[, orthproj, ]
    }
    if (is.matrix(A) == TRUE) {
      a <- array(0, c(k, 2, 1))
      a[, , 1] <- A[, orthproj]
      A <- a
    }

    out <- defplotsize2(A)
    width <- out$width
    if (length(c(B)) != 1) {
      if (length(dim(B)) == 3) {
        B <- B[, orthproj, ]
      }
      if (is.matrix(B) == TRUE) {
        a <- array(0, c(k, 2, 1))
        a[, , 1] <- B[, orthproj]
        B <- a
      }
      ans <- defplotsize2(B)
      width <- max(out$width, ans$width)
    }
    n <- dim(A)[3]
    lc <- length(color)
    lt <- k * m * n/lc
    color <- rep(color, times = lt)
    lc <- length(symbol)
    lt <- k * m * n/lc
    symbol <- rep(symbol, times = lt)

    plot(A[, , 1], xlim = c(out$xl, out$xl + width), ylim = c(out$yl,
                                                              out$yl + width), type = "n", xlab = " ",
         ylab = " ",pch=21,bg=(1:ncol(A)))
    grid()
    abline(v=0,h=0,col="gray50",lty=2)


    for (i in 1:n) {
      select <- ((i - 1) * k * m + 1):(i * k * m)
      points(A[, , i], pch = symbol[select], col = color[select])
      lines(A[joinline, , i])
      # text(A + c(-.2, .3),rownames(A),cex=0.9)
    }
    if (length(c(B)) != 1) {
      A <- B
      if (is.matrix(A) == TRUE) {
        a <- array(0, c(k, 2, 1))
        a[, , 1] <- A
        A <- a
      }
      out <- defplotsize2(A)
      n <- dim(A)[3]
      plot(A[, , 1], xlim = c(ans$xl, ans$xl + width),
           ylim = c(ans$yl, ans$yl + width), type = "n",
           xlab = " ", ylab = " ",pch=21,bg=(1:nrow(A)))
      grid()
      abline(v=0,h=0,col="gray50",lty=2)



      for (i in 1:n) {
        points(A[, , i], pch = symbol[select], col = color[select])
        lines(A[joinline, , i])
        #text(A+c(-.2, .3),rownames(A),cex=0.9, col = c(1:nrow(A)))
      }
    }
  }
}
#사용자 함수
# plotShapes(mds_online$cmdscale,mds_gachon$cmdscale,color=2,symbol = 2,joinline=c(1,4,5,2,3,6,1))

#원래함수
# library(shapes)
# plotshapes(mds_online$cmdscale,mds_gachon$cmdscale,color=3,symbol = 21,joinline=c(1,4,5,2,3,6,1))





#형상 분석
#nonmetric- metric shape analysis ~ procrustes analysis
#이름을 지정해주어야 함.
vv=c("ZOOM","Webex","googleMEET","Youtube","naveBand","prism")
OPA_Plot <- function(metricD,nonmetricD,title="",row_name,line=c(1:6,1),...) {
  if(is.character(title)==FALSE){
    stopifnot("input title ")
  }
  library(shape)
  #metric(뒤) onto nonmetirc(앞)
  pnm<- procOPA(nonmetricD, metricD)
  #nonmetric on metirc
  pmn<- procOPA(metricD, nonmetricD)


  OPA_data_nm <- pnm
  OPA_data1  <- OPA_data_nm$Bhat
  OPA_OSS_nm    <- OPA_data_nm$OSS
  OPA_rmsd_nm   <- OPA_data_nm$rmsd

  OPA_data_mn <- pmn
  OPA_data2  <- OPA_data_mn$Bhat
  OPA_OSS_mn   <- OPA_data_mn$OSS
  OPA_rmsd_mn  <- OPA_data_mn$rmsd

  rownames(OPA_data1 )= row_name
  rownames(OPA_data2 )= row_name


  #metric shape coordinated
  par(mfrow = c(1,2))
  plotshapes( OPA_data1, joinline=c(line,...),color="red",symbol = 1)#mtericOPA 형상
  abline(v=0,h=0)
  points(OPA_data1,pch=21,bg=c(line,...),cex=1.5)
  grid()
  text(OPA_data1+c(-.2, .1),rownames(OPA_data1),cex=0.9, col = c(rep("blue",3), rep("red",3)))
  title(paste("shape analysis::", title),paste("OPSS =", round(OPA_OSS_nm),2),
        paste("RMSD = ",round(OPA_rmsd_nm,2)))


  #Non-metric shape coordinated
  plotshapes( OPA_data2, joinline=c(line,...),color="red",symbol = 1)#mtericOPA 형상
  abline(v=0,h=0)
  points(OPA_data2,pch=21,bg=c(line,...),cex=1.5)
  grid()
  text(OPA_data2+c(-.2, .1),rownames(OPA_data2),cex=0.9, col = c(rep("blue",3), rep("red",3)))
  title(paste("shape analysis::", title),paste("OPSS =", round(OPA_OSS_mn),2),
        paste("RMSD = ",round(OPA_rmsd_mn,2)))





  par(mfrow=c(1,1))

  tabnm <-cbind(OPA_OSS_nm,OPA_rmsd_nm)
  tabmn <-cbind(OPA_OSS_mn,OPA_rmsd_mn)
  tablefull<-rbind(tabnm, tabmn)
  rownames(tablefull)=c("n -> m", "m -> n")
  colnames(tablefull)=c("OPSS","RMSD")

  res<-list(PNM=OPA_data1,PMN=OPA_data2,OPSS_RMSD=tablefull)
  res
}

#metric MDS와 Nonmetric MDS형상 비교 ######
#자료를 활용할 것 MDS분석
#mds_SDS_plot(on_mat) 를 이용하여 MDS시행
#match A onto B
# OPA_Plot(mds_online$cmdscale, mds_online$isoMDS$points,"Total OPA",vv,1,2,5,4,6,3,1)
# OPA_Plot(mds_gachon$cmdscale, mds_gachon$isoMDS$points,"gachon OPA",vv,1,2,5,4,6,3,1)
# OPA_Plot(mds_yonsei$cmdscale, mds_yonsei$isoMDS$points,"yonsei OPA",vv,1,2,5,4,6,3,1)
# OPA_Plot(mds_academy$cmdscale, mds_academy$isoMDS$points,"academy OPA",vv,1,2,5,4,6,3,1)
# OPA_Plot(mds_grad$cmdscale, mds_grad$isoMDS$points,"grad OPA",vv,1,2,5,4,6,3,1)


#카이제곱 오류해결######
#The code for Monte Carlo simulation is a C translation of the Fortran algorithm of Patefield (1981).
#Hope, A. C. A. (1968). A simplified Monte Carlo significance test procedure. Journal of the Royal Statistical Society Series B, 30, 582–598. http://www.jstor.org/stable/2984263.
#Agresti, A. (2007). An Introduction to Categorical Data Analysis, 2nd ed. New York: John Wiley & Sons. Page 38.

advanced.chisq.test<-function(cri, option=TRUE, ...){
  cc <- cri[rowSums(cri)>0,]
  if(option == TRUE){
    cat(paste0("정확검정을 위해 Monte Carlo simulation이 ",...,"회 실행되었습니다(박중희, 2020).","\n"))
    chisq.test(cc, simulate.p.value = TRUE, ...)
  }
  else{chisq.test(cc)
    cat(" Monte Carlo simulation, option=TRUE, B=회수로 실행하세요. \n")}
}
#advanced.chisq.test(online_freq,F, B=5000)


#Dplyr해결 #####
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)



# unload tidyverse and ggplot2 as an example first
detach("package:tidyverse", unload = TRUE)
detach("package:ggplot2", unload = TRUE)
# magic option
options(tidyverse.quiet = TRUE)
library(tidyverse)

#text mining####

#text를 띄어쓰기로 분리하기
word_split0<- function(x){
  library(stringr)
  wordsplit <-unlist(str_extract_all(x, boundary("word")))
  return(wordsplit)
}



#자연어 분석 1- 명사추출+ggplot, wordcloud
word_split<- function(x){
  library(stringr)
  library(KoNLP)
  library(tidyverse)
  library(ggplot2)
  library(knitr)
  library(wordcloud2)
  x <-str_remove_all(x,"\\,")
  xx <- extractNoun(x)
  x1 <- xx %>% unlist() %>% table() %>%  as.data.frame(stringsAsFactors=F)
  #x1 %>% rename(word=".",n=Freq)
  colnames(x1)=c("word","n")  #뱐수명
  x2<-filter(x1,nchar(word)>=2)  #두글자 이상
  x2<-x2 %>% arrange(desc(n))  #빈도수별로 정렬
  #table.x2 <-x2 %>% kable("pandoc",caption = "명사별 빈도수")

  x3 <- x2 %>% mutate(word = fct_reorder(word,n,"mean")) #순서 정렬


  a1 <-ggplot(data= x3[1:20,], aes(x=word, y=n),fill=word)+geom_bar(stat="identity")+coord_flip()+
    labs(y="빈도수", x= "주요 키워드 ")

  a2<-wordcloud2(x2,size = 2,
                 color = "random-light",backgroundColor = "black", shape = 'circle')


  result <-list(Noun_extract=xx, Noun_count= x2,Bar_Plot=a1,WordCloud=a2)
  result
}


#자연어 분석 2- 형태소로 추출 ggplot, wordcloud
word_split2<- function(x){
  library(KoNLP)
  #useNIADic()
  #useSejongDic()   #사전 설치
  library(tidyverse)           #install.packages("dplyr")
  library(stringr)           #install.packages("stringr")
  library(wordcloud2)        #install.packages("wordcloud22")
  library(reshape2)          #install.packages("reshape2")


  xxx <-x
  txt_df <-  xxx%>% SimplePos09 %>%melt %>%as_tibble
  txt_s_df<-txt_df %>%  select(3,1)


  txt_count <-txt_s_df %>%
    mutate(noun=str_match(value, '([가-힣]+)/N')[,2]) %>%
    na.omit %>%
    filter(str_length(noun)>=2) %>%
    count(noun, sort=TRUE)

  txt_count <- txt_count %>% mutate(noun = fct_reorder(noun,n,"mean")) #순서 정렬
  table<-ggplot(data= txt_count[1:20,], aes(x=noun, y=n),fill=noun)+geom_bar(stat="identity")+coord_flip()+
    labs(y="빈도수", x="주요 키워드 ")
  wdc<-wordcloud2(txt_count,size = 2,
                  color = "random-light",backgroundColor = "black", shape = 'circle')

  list(part_of_Speech= txt_df,Count_data=txt_count,BarPlot=table,WodrCloud=wdc)
}


#형태소로 분석3 데이터만
word_split3<- function(x){
  library(KoNLP)
  #useNIADic()
  #useSejongDic()   #사전 설치
  library(tidyverse)           #install.packages("dplyr")
  library(stringr)           #install.packages("stringr")
  library(wordcloud2)        #install.packages("wordcloud22")
  library(reshape2)          #install.packages("reshape2")



  xxx <-x
  txt_df <-  xxx%>% SimplePos09 %>%melt %>%as_tibble
  txt_s_df<-txt_df %>%  select(3,1)


  txt_count <-txt_s_df %>%
    mutate(noun=str_match(value, '([가-힣]+)/N')[,2]) %>%
    na.omit %>%
    filter(str_length(noun)>=2) %>%
    count(noun, sort=TRUE)

  txt_count <- txt_count %>% mutate(noun = fct_reorder(noun,n,"mean")) #순서 정렬
  #table<-ggplot(data= txt_count[1:20,], aes(x=noun, y=n),fill=noun)+geom_bar(stat="identity")+coord_flip()+
  #  labs(y="빈도수", x="주요 키워드 ")
  #wdc<-wordcloud2(txt_count,size = 2,
  #                color = "random-light",backgroundColor = "black", shape = 'circle')

  result<- list(part_of_Speech= txt_df,Count_data=txt_count)#,BarPlot=table,WodrCloud=wdc)
  result
}



#사전에 단어를 추가 - 분석이 잘 안될 경우 명사를 추가하는 기능
dicAdd<- function(add_Keyword){
  library(KoNLP)
  #useNIADic()
  #useSejongDic()   #사전 설치
  buildDictionary(user_dic = data.frame(add_Keyword,rep("ncn",length(add_Keyword))),replace_usr_dic = T)
}





##텍스트 마이닝 #####
#tokens화 #####

Token_gen <- function(data){
  library(tidyverse)
  library(tidytext)
  data1 <-data %>% tibble(line=1:length(data))%>%
    unnest_tokens(word,text) %>%
    count(word,sort = T) %>% filter(str_length(word)>=2)
  data1
}

##word_kor--token화  한글단어추출 함수 ####
word_kor<- function(data){
  library(KoNLP)
  library(stringr)
  colnames("text")
  #SimplePos09() 형태소 구분
  data %>% mutate(word=str_match(text, '([가-힣]+)')[,2]) %>%
    na.omit %>%
    filter(str_length(word)>=2) %>%  count(word, sort=TRUE)
}

#불용어 제거함수: 데이터와 불용어를 넣어서 처리
#데이터, 컬럼, 제거할 단어
#자체네 tibble처리후에 변수를 reorder
word_remover <- function(data,remove_words){
  library(tidyverse)
  library(tidytext)
  colnames(data)=c("word","n")
  word_remove<-tibble(word=remove_words)
  data <-data %>% anti_join(word_remove)
  data<-data %>% mutate(word=reorder(word,n))
  data
}

#빈도수 그래프 ####
freq_bar<- function(data){
  library(ggplot2)
  library(dplyr)
  colnames(data)=c("word","n")
  data1 <- data %>% mutate(word=reorder(word,n)) #순서 정렬
  ggplot(data=data1[1:30,], aes(x=word, y=n))+geom_bar(stat="identity") + geom_col()+coord_flip()+
    labs(y="빈도수", x="주요 키워드 ")
}
