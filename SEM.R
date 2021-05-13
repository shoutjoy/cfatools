
library(tidyverse)
library(psych)
library(lavaan)
library(semPlot)
library(semTools)
library(readr)
library(PerformanceAnalytics) #chart.Correlation


#기술통계
Summarise <- function(x,type="col",digit=2){

  #col:세로로 나타내기 (변수가 세로)
  #row:가로로 나타내기 (변수가 가로)
  s_col <-summarise(x, across(
  .cols=c(1:ncol(x)),
  .fns = function(x){str_c(round(mean(x),digit),"( ", round(sd(x),digit)," )")}
  )) %>% t() %>% as.data.frame()
  colnames(s_col)="mean(sd)"
  s_row<- s_col %>% t() %>% data.frame()
  switch(type,
        col=s_col,
        row=s_row)
  }

#SEM process
SEM <- function(x, type="cfa"){

  library(dplyr)
  library(knitr)
  library(lavaan)
  library(semTools)
  library(tibble)
  library(semPlot)

  tryCatch({
  switch(type,
         cfa=cfa2(x),
         CR=CR(x),
         AVE=AVE(x),
         semPaths=SemPaths(x),
         fit=CompareFit(x),
         loadings=factor_loadings(x),
         convergent=Convergent_Validity(x),
         discriminant=Discriminant_Validity(x),
         effect=effect1(x),
         med=med_effect(x),
         label=effect2(x),
         regress=effect3(x),
         define =effect4(x),
         group=effect5(x)
         )
  },error=function(e) return("입력 오류입니다.cfa,CR,AVE,semPaths, fit, loadings,convergent,  discriminant,  effect, med,  label, regress, define,  group ")
)
}




#cfa과정에 대한 정리
cfa2 <- function(x, format="pandoc"){

  library(dplyr)
  library(knitr)
  library(lavaan)
  library(semTools)
  library(tibble)
  library(semPlot)

  # tryCatch({

  #01 fit table
  options(scipen = 100)

  fit.indices=c("chisq","pvalue", "df","rmsea",
                "gfi","agfi","srmr","cfi","tli","nfi","aic","bic")
  fitMeasures <- round(fitMeasures(x,fit.indices),3)
  fitMeasures_s <- round(fitMeasures(x,fit.indices),3)


  # fitMeasures <- as.data.frame(fitMeasures(x,fit.indices))
  fitMeasures <- as.data.frame(fitMeasures)
   fitMeasures$critera <- c("",
                           "*p.value >= 0.05",
                           "_chisq/df <= 3(<5(ok)",
                           "*RMSEA< 0.05",
                           "*GFI >= 0.95",
                           "_AGFI>= 0.90",
                           "*SRMR < 0.08",
                           "*CFI >= 0.95",
                           "_TLI >= 0.90",
                           "_NFI >= 0.90",
                           "_lower",
                           "_lower")
  fitMeasures$Ref <-c("-",
                      "-",
                      "Wheaton et al.(1977)",
                      "Kline(2011)",
                      "Kline(2011)",
                      "Tanaka & Huba(1985)",
                      "Kline(2011)",
                      "Kline(2011)",
                      "Bentler & Bonett(1980)",
                      "Bollen(1989)",
                      "Akaike(1973)",
                      "-")
  fitMeasures$chiq_df <- c("","",
                           round(fitMeasures[1,1]/fitMeasures[3,1],2),
                           "","","","","","","","","")
  # fitMeasures$fit_chek  <- c("absolute fit","",
  #                            "absolute fit",
  #                            "absolute fit ",
  #                            "absolute fit ",
  #                            "absolute fit ",
  #                            "absolute fit ",
  #                            "incremental fit",
  #                            "incremental fit",
  #                            "incremental fit",
  #                            "parsimonious fit",
  #                            "parsimonious fit")
  fit <- fitMeasures  %>%
    kable(digits=3, format=format,
    caption="FitMeasure and criterian
          (*)동시에 만족해야하는 것 표시 By kline(2011)")



  #modelfit
  fitdata <- fitMeasures(x,c("chisq","df","pvalue",
                             "rmsea","rmsea.ci.lower",
                             "rmsea.ci.upper","rmsea.pvalue",
                             "gfi","cfi","srmr"))

  criteria_data = c("chisq","df","pvalue>0.05",
                    "rmsea<0.08","90%CI.lower","90%CI.upper","p value",
                    "gfi>0.95", "cfi>0.95","srmr<0.08")

  modelfitdata <-cbind("criterian"=criteria_data,
                       "Value"=round(fitdata,3))

  fitMeasures_s1 <- modelfitdata %>%
    kable (format=format,
           caption = "01 Model fit information")




  #04 factor loading
  options(knitr.kable.NA="")
  factorloading <- parameterEstimates(x, standardized=TRUE) %>%
    filter(op=="=~") %>%
    mutate(stars=ifelse(pvalue < 0.001, "***",
                        ifelse(pvalue < 0.01, "**",
                               ifelse(pvalue < 0.05, "*", "")))) %>%
    mutate(label=ifelse(std.all>0.7,"Yes(Good)",
                        ifelse(std.all>0.5,"Yes(fair)","No"))) %>%
    dplyr::select("Latent"=lhs, Item=rhs, Est=est,S.E.=se,
                  Z=z, Sig.=stars, "p"=pvalue,
                  std=std.all, beta_Accept=label) %>%
    kable(digits=3, format=format,
          caption="02 Indicator Validity(1)-Factor Loadings::
          (1) c.r(=Estimate/S.E) p<0.05,
          (2) std.damda >= 0.5(Bagozzi & Yi(1988)")




  #Cronbach alpha
  alpha.1 <- reliability(x,return.total = T) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::select( "Cronbach's alpha"=alpha) %>%
    mutate(sig=ifelse(`Cronbach's alpha`>0.7,"Accept(>0.7) *",
                      ifelse(`Cronbach's alpha`>0.6,"Yes(poor) *", "Reject")))

  #," average variance extracted(AVE)"=avevar)

  #05 Reprort cronbach, AVE, C.R
  FL.1 <- cbind(alpha.1)
  FL<-FL.1%>%kable(digits=3, format=format,
                   caption="03-1. Internal consistency
                   (Cronbach's Alpha, 1951)")




  #CR,AVE

  AVE <- reliability(x,return.total = F) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::select( "AVE"=avevar)
  sqrt.AVE <- sqrt(AVE)
  colnames(sqrt.AVE)="sqrt.AVE"

  #correlations Matrix
  rho <- inspect(x,"std")$beta


  {#reliablity Construct Reliability & Average Variance Extracted#
  #AVE cal
#
#   l.matrix<- inspect(x,"std")$lambda
#   l.matrix[l.matrix==0]<-NA
#   AVE<-apply(l.matrix^2,2,mean,na.rm=T)
#
#   #sqrt.AVE
#   sqrt.AVE <- sqrt(AVE)
#   #latent correlation
  # rho <- inspect(x,"std")$psi
#   #Composit reliability
#   t.matrix<- inspect(x,"std")$theta
#   t.matrix[t.matrix==0]<-NA
#   t.matrix
#   cr1<-apply(l.matrix,2,sum,na.rm=T)
#   d.sum <-apply(1-l.matrix^2,2,sum,na.rm=T)
#   C.R<- cr1^2/(cr1^2+d.sum)
}


# Convergent validity
  alpha_AVE_CR <-  reliability(x,return.total = T) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::select("Cronbach"=alpha, "CR" = omega3, "AVE"=avevar) %>%
    mutate(sqrt.AVE=sqrt(AVE))%>%
     mutate(CR_check=ifelse(CR>0.7,"Accept(>0.7) *","Reject")) %>%
    mutate(AVE_check=ifelse(AVE>0.5,"Accept(>0.5) *","Reject"))%>%
    kable(digits = 3,format = format,
          caption = "03 Convergent validity
          AVE(>0.5) & CR(>0.7): Fornell & Lacker(1981)")




  #check data####


  #06 discriminant validity
  betaa <-inspect(x, "std")$beta

   if(is.null(betaa)){

    psi <-inspect(x, "std")$psi
    psi[lower.tri(psi)==FALSE]<-0

    rho1<- psi %>% as.data.frame()
    rho1$max<- apply(rho1,1,max)
    diff<- cbind(rho1$max,sqrt.AVE)
    diff$delta<-diff[,2]- diff[,1]
    diff$sig<-ifelse(diff$delta>0,"*","")

  FornellNacker <-cbind(psi, sqrt.AVE,sig=diff[,4]) %>% as.data.frame()


  validity <- FornellNacker %>%
      kable(digits=3, format=format,
            caption="04 Discriminant Validity:
          rho < Square Root of(AVE)
           By Fornell & Lacker(1981)")

  }else{

  rho1<-inspect(x, "std")$beta %>% as.data.frame()
  rho1$max<- apply(rho1,1,max)
  diff<- cbind(rho1$max,sqrt.AVE)
  diff$delta<-diff[,2]- diff[,1]
  diff$sig<-ifelse(diff$delta>0,"*","")

  FornellNacker <-cbind(rho, sqrt.AVE,sig=diff[,4]) %>% as.data.frame()


  validity <- FornellNacker %>%
    kable(digits=3, format=format,
          caption="04 Discriminant Validity:
          rho < Square Root of(AVE)
           By Fornell & Lacker(1981)")

  }



  lv.cor <-lavInspect(x, what="cor.lv")

  lv.cor.sig <- parameterEstimates(mmodel_sem1, standardized = T) %>%
    filter(op=="~") %>% select(lhs,op,rhs, std.all, pvalue) %>%
    mutate(sig=ifelse(pvalue < 0.001, "***",
                        ifelse(pvalue < 0.01, "**",
                               ifelse(pvalue < 0.05, "*", "")))) %>%
    kable(digits=3, format=format,
          caption="05 latent correlation Significant Check")


  all.reuslt <-list(fit_criterian=fit,
                    model_fit=fitMeasures_s1,
                    factorloadings=factorloading,
                    Cronbach=FL,
                    Convegent=alpha_AVE_CR,
                    Discriminant=validity,
                    betaMatrix=lv.cor,
                    betaMat_sig=lv.cor.sig )
  all.reuslt

  }



SemPaths<- function(x,lat=9,lat2=6,man=7,man2=3,
                         edgelabelcex=1,label.cex=1,
                    WhatLabels="par",
                         rotate=1,
                         layout="tree2",
                         Curve=1.2,
                         curvePivot=T,
                         style="lisrel",
                         edge_label_pos=.6,
                         resid=5,intTF=T,intercepts=F,
                         strut_T_F=F,
                         bg="gray84",
                         nDigits=3){
  library(semPlot)

  semPaths(x, whatLabels = WhatLabels, nCharNodes = 8, #text number
           rotation = rotate,intercepts = intercepts,
           style = style,residScale = resid,
           curvePivot=curvePivot,curve = Curve,
           layout = layout, #shape
           layoutSplit = F ,subScale = 1,
           subScale2 = 1, # manifest 1 row, 2column
           subRes = 4, #Default=4
           sizeLat = lat,sizeLat2 = lat2 ,
           sizeMan = man,sizeMan2 = man2,
           color = list(lat="skyblue", man="Gold", int="gray80")
           ,label.cex=label.cex,
           edge.label.cex = edgelabelcex, edge.color = "gray25",
           edge.label.position= edge_label_pos,
           edge.width= 1.6, asize=1.9,
           mar=c(6,6,6,6), bg=bg,
           structural = strut_T_F,
           nDigits = nDigits)
}


CompareFit <- function(..., option=1,format="pandoc") {
  library(magrittr)
  library(stargazer)
  library(tibble)
  library(knitr)
  library(dplyr)


  m <- list(...)


  if(option==1){
  result <-sapply(m, fitMeasures) %>%
    set_colnames(paste0("Model-", 1:length(m))) %>%
    as.data.frame() %>%
    rownames_to_column("Fit_Measures") %>%
    slice(match(c("chisq", "df", "pvalue",
                  "rmsea","rmsea.ci.upper","rmsea.ci.lower",
                  "gfi", "srmr","cfi","tli","aic","bic"),
                Fit_Measures)) %>%
    mutate(Fit_Measures=c("Chi-square", "df", "p-value",
                          "RMSEA(<0.05)","rmsea.ci.upper","rmsea.ci.lower",
                          "GFI(>0.95)","SRMR(<0.08)",
                          "CFI(>0.95)","TLI(NNFI)(>0.9)","AIC","BIC"))%>%
    kable(digits=3, format=format,
          caption="model comparison")
  }else{
  if(option==2){
    result <-sapply(m, fitMeasures) %>%
      set_colnames(paste0("Model-", 1:length(m))) %>%
      as.data.frame() %>%
      rownames_to_column("Fit_Measures") %>%
      slice(match(c("chisq", "df", "pvalue",
                    "rmsea","rmsea.ci.upper","rmsea.ci.lower",
                    "gfi", "srmr","cfi","tli","aic","bic"),
                  Fit_Measures)) %>%
      mutate(Fit_Measures=c("Chi-square", "df", "p-value",
                            "RMSEA(<0.05)","rmsea.ci.upper","rmsea.ci.lower",
                            "GFI(>0.95)","SRMR(<0.08)",
                            "CFI(>0.95)","TLI(NNFI)(>0.9)",
                            "AIC","BIC"))


  }
    else{ return("option = 1, 2중 선택")}
  }
  return(result)

}


CompareFit_2 <- function(..., format="pandoc") {
  library(magrittr)
  library(stargazer)
  library(tibble)
  library(knitr)
  library(dplyr)


  m <- list(...)



    result <-sapply(m, fitMeasures) %>%
      set_colnames(paste0("Model-", 1:length(m))) %>%
      as.data.frame() %>%
      rownames_to_column("Fit_Measures") %>%
      slice(match(c("chisq", "df", "pvalue",
                    "rmsea","rmsea.ci.upper","rmsea.ci.lower",
                    "gfi", "srmr","cfi","tli","aic","bic"),
                  Fit_Measures)) %>%
      mutate(Fit_Measures=c("Chi-square", "df", "p-value",
                            "RMSEA(<0.05)","rmsea.ci.upper","rmsea.ci.lower",
                            "GFI(>0.95)","SRMR(<0.08)",
                            "CFI(>0.95)","TLI(NNFI)(>0.9)","AIC","BIC"))
    result

}




#01 CFA:factor loading####
factor_loadings <- function(x,format="pandoc"){

  #factor loading

  options(knitr.kable.NA="")

  factorloading <- parameterEstimates(x, standardized=TRUE) %>%
    filter(op=="=~") %>%
    mutate(stars=ifelse(pvalue < 0.001, "***",
                        ifelse(pvalue < 0.01, "**",
                               ifelse(pvalue < 0.05, "*", "")))) %>%
    mutate(Accept=ifelse(std.all>0.7,"Yes(**)",ifelse(std.all>0.5,"Yes( *)","No"))) %>%
    dplyr::select("Latent"=lhs,op, Item=rhs, Est=est, SE=se,
                  "cr(t)"=z, Sig.=stars, #"p-value"=pvalue,
                  std=std.all,
                  Accept) %>%
    kable(digits=3, format=format, caption="Factor Loadings::
          (1) c.r(=Estimate/S.E) p<0.05,
          (2) std.damda >= 0.5")

  factorloading

}

# fit.math2 %>% factor_loadings()

#CR-수기계산 ################
CR <- function(x,digit=3){
  library(knitr)
  library(dplyr)
  inspect(x,"std")
  #lambda
  inspect(x,"std")$lambda
  l.matrix<- inspect(x,"std")$lambda
  l.matrix[l.matrix==0]<-NA
  #theta: error coefficient
  t.matrix<- inspect(x,"std")$theta
  t.matrix[t.matrix==0]<-NA

  #t.matrix
  cr1<-apply(l.matrix,2, sum, na.rm=T)
  d.sum <-apply(1-l.matrix^2,2,sum,na.rm=T)
  CR.c<- cr1^2/(cr1^2 + d.sum)
  CR <-CR.c %>% kable(digits=digit, format="pandoc", caption="Construct Reliability > 0.7(Bagozzi and Yi 1988)")
  return(CR)}

# obj.mycfa1 %>% CR(5)
# obj.mycfa2 %>% CR(5)

#AVE수기계산#####
AVE<-function(x, digit=3){
  library(knitr)
  l.matrix <- inspect(x,"std")$lambda
  l.matrix[l.matrix==0]<-NA   #불필요한 것은 삭제
  AVE.1<-apply(l.matrix^2, 2, mean, na.rm=T)
  AVE<- AVE.1 %>% kable(digits=digit, format="pandoc",
                        caption="AVE(average variance extracted)>0.5")
  return(AVE)
}


#02 CFA Convergent_Validity####
Convergent_Validity  <- function(x,format="pandoc", digit=3){
  library(dplyr)
  library(knitr)
  library(semTools)

    rel.1<- reliability(x) %>% t() %>%
      data.frame() %>%
      select(Cronbach=alpha, CR=omega3, AVE=avevar) #%>% print(digits=3)

    rel.2<-rel.1 %>% mutate(CR_check=ifelse(CR>0.7,"Accept(>0.7) *","Reject")) %>%
      mutate(AVE_check=ifelse(AVE>0.5,"Accept(>0.5) *","Reject"))

    rel<- rel.2 %>% dplyr::select(CR,CR_check,AVE,AVE_check) %>%
      kable(digits=digit, format=format,
            caption="집중타당도(Convergent Validity)- Hair(2009):
            CR(>0.7) and AVE(>0.5)")
    rels<- list(Reliability= round(rel.1,digit),
                Convergent_Validity=rel)
    rels
}

# fit %>% Convergent_Validity()

#

# #03CFA Discriminant####
# Discriminant_Validity<-function(x, format="pandoc"){
#   library(knitr)
#   library(psych)
#
#   library(semTools)
#   rel.1<- reliability(x) %>% t() %>%
#     data.frame() %>%
#     select(Cronbach=alpha, CR=omega3, AVE=avevar)
#
#   #Internal consistency
#   FornellNacker1 <- cbind(rel.1, sqrt.AVE)%>%
#     kable(digits=3, format=format)
#
#
#   AVE <- reliability(x,return.total = F) %>%
#     t() %>%
#     as.data.frame() %>%
#     dplyr::select( "AVE"=avevar)
#   sqrt.AVE <- sqrt(AVE)
#   colnames(sqrt.AVE)="sqrt.AVE"
#   # colnmaes(sqrt.AVE)=c("sqrtAVE")
#   #latent correlation
#   # rho <- inspect(x,"std")$beta #+inspect(x,"std")$psi
#   # lavInspect(mmodel_sem1, what="cor.lv")
#
#   rho1<-inspect(x, "std")$beta %>% as.data.frame()
#   rho1$max<- apply(rho1,1,max)
#   diff<- cbind(max_rho=rho1[,"max"],sqrt.AVE)
#   diff$delta<-diff[,2]- diff[,1]
#   diff$sig<-ifelse(diff$delta>0,"*","")
#
#   FornellNacker <-cbind(rho, sqrt.AVE,sig=diff[,4]) %>% as.data.frame()
#
#
#   validity <- FornellNacker %>%
#     kable(digits=3, format=format,
#           caption="04 Discriminant Validity:
#           rho < Square Root of(AVE)
#            By Fornell & Lacker(1981)")
#
#
#
#   result.1 <- list(CR_AVE=FornellNacker1, discriminant_validity=validity)
#   result.1
#
# }
#
#
#

Discriminant_Validity<-function(x, format="pandoc"){
  library(knitr)
  library(psych)

  library(semTools)

  rel.1<- reliability(x) %>% t() %>%
    data.frame() %>%
    select(Cronbach=alpha, CR=omega3, AVE=avevar)

  #Internal consistency
  FornellNacker1 <- cbind(rel.1, sqrt.AVE)%>%
    kable(digits=3, format=format)


  AVE <- reliability(x,return.total = F) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::select( "AVE"=avevar)
  sqrt.AVE <- sqrt(AVE)
  colnames(sqrt.AVE)="sqrt.AVE"
  # colnmaes(sqrt.AVE)=c("sqrtAVE")
  #latent correlation
  # rho <- inspect(x,"std")$beta #+inspect(x,"std")$psi
  # lavInspect(mmodel_sem1, what="cor.lv")



  #06 discriminant validity
  betaa <-inspect(x, "std")$beta

  if(is.null(betaa)){

    psi <-inspect(x, "std")$psi
    psi[lower.tri(psi)==FALSE]<-0

    rho1<- psi %>% as.data.frame()
    rho1$max<- apply(rho1,1,max)
    diff<- cbind(rho1$max,sqrt.AVE)
    diff$delta<-diff[,2]- diff[,1]
    diff$sig<-ifelse(diff$delta>0,"*","")

    FornellNacker <-cbind(psi, sqrt.AVE,sig=diff[,4]) %>% as.data.frame()


    validity <- FornellNacker %>%
      kable(digits=3, format=format,
            caption="04 Discriminant Validity:
          rho < Square Root of(AVE)
           By Fornell & Lacker(1981)")

  }else{

    rho1<-inspect(x, "std")$beta %>% as.data.frame()
    rho1$max<- apply(rho1,1,max)
    diff<- cbind(rho1$max,sqrt.AVE)
    diff$delta<-diff[,2]- diff[,1]
    diff$sig<-ifelse(diff$delta>0,"*","")

    FornellNacker <-cbind(rho, sqrt.AVE,sig=diff[,4]) %>% as.data.frame()


    validity <- FornellNacker %>%
      kable(digits=3, format=format,
            caption="04 Discriminant Validity:
          rho < Square Root of(AVE)
           By Fornell & Lacker(1981)")
  }

  result.1 <- list(Convegent=FornellNacker1, Discriminant_validity=validity)
  result.1


}









effect<- function(x,type="basic"){

  library(dplyr)
  library(stargazer)

tryCatch({
  switch(type,
         basic=effect1(x),
         label=effect2(x),
         regress=effect3(x),
         effect=effect4(x),
         group=effect5(x))
},error=function(e) {
  return("입력 오류입니다.
         정확한 명칭을 입력해주세요(basic, label, regress, effect, group) ")})

}

effect1 <-function(x){
  library(dplyr)
  library(stargazer)
  parameterEstimates(x, standardized = T, rsquare = T) %>%
    filter(op=="~"|op==":=") %>%
    mutate(stars=ifelse(pvalue<0.001,"***",
                        ifelse(pvalue<0.01,"**",
                               ifelse(pvalue<0.05,"*","")))) %>%
    mutate(op=ifelse(op=="~","<--",
                     ifelse(op==":=","effect",""))) %>%
    dplyr::select(Dependent=lhs,Path=op, Independent=rhs,est,
                  "Coefficient"= std.all, c.r=z,
                  Sig.=stars, "p-value"=pvalue) %>%
    stargazer(type="text", title="Regression & Total Effect .",
              summary = FALSE,
              digits = 3, digits.extra = 0, rownames = FALSE)

  }

#label표시 ####
effect2 <-function(x){
  library(dplyr)
  library(stargazer)
  parameterEstimates(x, standardized = T, rsquare = T) %>%
    filter(op=="~"|op==":=") %>%
    mutate(stars=ifelse(pvalue<0.001,"***",
                        ifelse(pvalue<0.01,"**",
                               ifelse(pvalue<0.05,"*","")))) %>%
    mutate(op=ifelse(op=="~","<--",
                     ifelse(op==":=","effect",""))) %>%
    dplyr::select(Dependent=lhs,Path=op, Independent=rhs,
                  label,est,"Coefficient"= std.all, c.r=z,
                  Sig.=stars, "p-value"=pvalue) %>%
    stargazer(type="text", title="Regression & Total Effect (label) .", summary = FALSE,
              digits = 3, digits.extra = 0, rownames = FALSE)

  }


#effect3 regress만 표시 ####
effect3 <-function(x){
  library(dplyr)
  library(stargazer)
  parameterEstimates(x, standardized = T, rsquare = T) %>%
    filter(op=="~") %>%
    mutate(stars=ifelse(pvalue<0.001,"***",
                        ifelse(pvalue<0.01,"**",
                               ifelse(pvalue<0.05,"*","")))) %>%
    mutate(op=ifelse(op=="~","<--","" )) %>%
    dplyr::select(Dependent=lhs,Path=op,
                  Independent=rhs,
                  est,
                  "Coefficient"= std.all,
                  c.r=z,
                  Sig.=stars,
                  "p-value"=pvalue) %>%
    stargazer(type="text", title="Regression", summary = FALSE,
              digits = 3, digits.extra = 0, rownames = FALSE)

}

#:= total effect####
effect4 <-function(x){
  library(dplyr)
  library(stargazer)
  parameterEstimates(x, standardized = T, rsquare = T) %>%
    filter(op==":=") %>%
    mutate(stars=ifelse(pvalue<0.001,"***",
                        ifelse(pvalue<0.01,"**",
                               ifelse(pvalue<0.05,"*","")))) %>%
    mutate(op=ifelse(op==":=","effect","")) %>%
    dplyr::select(Dependent=lhs,Path=op, Independent=rhs,label,est,
                  "Coefficient"= std.all, c.r=z,
                  Sig.=stars, "p-value"=pvalue) %>%
    stargazer(type="text", title="Total Effect(Direct, Indirect)",
              summary = FALSE,
              digits = 3, digits.extra = 0, rownames = FALSE)

}


#effect5 group 표시 #####
effect5 <-function(x){
  library(dplyr)
  library(stargazer)
  parameterEstimates(x, standardized = T, rsquare = T) %>%
    filter(op=="~"|op==":=") %>%
    mutate(stars=ifelse(pvalue<0.001,"***",
                        ifelse(pvalue<0.01,"**",
                               ifelse(pvalue<0.05,"*","")))) %>%
    mutate(op=ifelse(op=="~","<--",
                     ifelse(op==":=","effect",""))) %>%
    dplyr::select(Dependent=lhs,Path=op, Independent=rhs,group,
                  "Est"=est,c.r=z,
                  Sig.=stars, "p-value"=pvalue) %>%
    stargazer(type="text", title="Regression & Total Effect (label) .",
              summary = FALSE,
              digits = 3, digits.extra = 0, rownames = FALSE)

}



# effect<- function(x){
#
#   library(dplyr)
#   library(stargazer)
#   standardizedsolution(x) %>%
#     filter(op=="~"|op==":=") %>%
#     mutate(stars=ifelse(pvalue<0.001,"***",
#                         ifelse(pvalue<0.01,"**",
#                                ifelse(pvalue<0.05,"*","")))) %>%
#     mutate(op=ifelse(op=="~","<--",
#                      ifelse(op==":=","effect",""))) %>%
#     dplyr::select(Dependent=lhs,Path=op, Independent=rhs, "Coefficient"= est.std, c.r=z,
#                   Sig.=stars, "p-value"=pvalue) %>%
#     stargazer(type="text", title="Regression & Total Effect .", summary = FALSE,
#               digits = 3, digits.extra = 0, rownames = FALSE)
# }



#effect함수
med_effect <- function(x, effect="",
                       caption="",
                       option=1,
                       effect2="",
                       format="markdown"
                       ){
  library(knitr)
  library(dplyr)
  if(option==1){
    parameterEstimates(x,standardized = T ) %>% filter(op==":=") %>%
      filter(str_detect(lhs,effect) & str_detect(lhs,effect2)) %>%
      mutate(sig=cut(pvalue,c(-Inf,0.001,0.01,0.05,1),
                     labels=c("***","**","*","Not Sig"))) %>%
      select(lhs, rhs, est,se, z, pvalue, sig) %>%
      kable(format, digits=3, caption = caption)
  }
  else
  {if(option==2){
    parameterEstimates(x,standardized = T ) %>% filter(op==":=") %>%
      filter(str_detect(lhs,effect) & str_detect(lhs,effect2)) %>%
      mutate(sig=cut(pvalue,c(-Inf,0.001,0.01,0.05,1),
                     labels=c("***","**","*","Not Sig"))) %>%
      select(lhs, rhs, est,se, z, pvalue,std.all, sig)
  }else{}
  }
}
#

# #모든결과 보기
# sats.sem2 %>% med_effect(effect = "E",1,effect2 = "")
#
# #indirect effect계산과정 없이 보기
# sats.sem2 %>% med_effect(effect = "E",1)
#
# # option ==1 table, option==2, data.frame()
# sats.sem2 %>% med_effect("DE", "직접효과 ",1,"")
# sats.sem2 %>% med_effect("IE", "간접효과 ")
# sats.sem2 %>% med_effect("TE", "총효과 ",2)





#두그룹을 분리했을 때 사용
Group_lodings_diff<- function(M1,M2,option=1,op1="=~",digit=3,
                              lhs1="",rhs1=""){
  library(knitr)
  library(dplyr)

  #simple compare
    para_0 <- parameterEstimates(M1) %>% select(lhs, op,rhs)
    para_a<-parameterEstimates(M1) %>% select(lhs, op,rhs,est,se)
    para_b<-parameterEstimates(M2) %>% select(lhs, op,rhs, est,se)
    para_compare<-data.frame(para_0,
                             round(para_a[,4:5],digits = digit),
                             round(para_b[,4:5],digits = digit))
    colnames(para_compare)=c("lhs","op","rhs","est.1","se.1","est.2","se.2")
    # para_compare
    para_compare$diff_est<-round(para_a[,4]-para_b[,4], digits = digit)
    para_compare %>% filter(op==op1)
      # kable("pandoc",3,caption = "그룹별 모수추정치 비교")
}

#mesurement invariance####

#loadings처럼 한번에 처리한 경우 사용

#option 1 basic
# opption 2  label &dataframe
# option 3 add pvalue

Group_invariance_diff_inter <- function(fit,
                                        option=1,
                                        op1="=~",
                                        block1=1,
                                        block2=2){
  library(dplyr)
 library(knitr)

  if(option==1){
    para_0<-parameterEstimates(fit) %>%
      filter(block==block1) %>%
      select(lhs,op,rhs)
    para_1<- parameterEstimates(fit) %>%
      filter(block==block1) %>%
      select(lhs,op,rhs,"Grp1_Est"=est,"G1_se"=se,pvalue)
    para_2 <- parameterEstimates(fit) %>%
      filter(block==block2) %>%
      select(lhs,op,rhs, "Grp2_Est"=est,"G2_se"=se,pvalue)
    para_compare<-data.frame(para_0, para_1[,4:6],para_2[,4:6])
    # para_compare
    para_compare$diff_est<-para_1[,4]-para_2[,4]

    a1<- para_compare%>%  filter(op == op1) %>%
      kable("markdown", 3, caption = "두 모수의 차이비교")


  }else{ #1

    if(option==2){
      para_0<-parameterEstimates(fit) %>%
        filter(block==block1) %>%
        select(lhs,op,rhs,label)
      para_1<- parameterEstimates(fit) %>%
        filter(block==block1) %>%
        select(lhs,op,rhs,"Grp1_Est"=est,"G1_se"=se)
      para_2 <- parameterEstimates(fit) %>%
        filter(block==block2) %>%
        select(lhs,op,rhs, "Grp2_Est"=est,"G2_se"=se)
      para_compare<-data.frame(para_0,round(para_1[,4:5],3),round(para_2[,4:5],3))
      # para_compare
      para_compare$diff_est<-round(para_1[,4]-para_2[,4],3)


      a1<-para_compare%>%
        filter(op== op1)

    }
    else{#2


      if(option==3){
        para_0<-parameterEstimates(fit) %>%
          filter(block==block1) %>%
          select(lhs,op,rhs)
        para_1<- parameterEstimates(fit) %>%
          filter(block==block1) %>%
          select(lhs,op,rhs,"Grp1_Est"=est,"G1_se"=se,pvalue)
        para_2 <- parameterEstimates(fit) %>%
          filter(block==block2) %>%
          select(lhs,op,rhs, "Grp2_Est"=est,"G2_se"=se, pvalue)
        para_compare<-data.frame(para_0,round(para_1[,4:6],3),round(para_2[,4:6],3))

        para_compare$diff_est<-round(para_1[,4]-para_2[,4],3)

        a1<-para_compare%>%
          filter(op== op1)
      }else{
        return("select! option=1 or 2(label) or 3(pvalue)")

      }


      }#else3
    }#else2
  # }#esle1


  res=list(group=fit,loadings=a1)
  res


}#function



