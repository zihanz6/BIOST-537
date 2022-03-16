library(data.table)
library(survival)
library(foreign)
library(msm)
library(flexsurv)

##load in dataset
dt <- fread('/ihme/homes/mwalte10/bmt.csv')



#####################################################
#Analysis
#####################################################
## Question 1
{
  
  ###cleaning variables-----
  df$X <- NULL
  #survival time
  df$ts 
  #death indicator
  df$deltas <- factor(df$deltas,
                      levels = c(0,1),
                      labels = c("Alive","Dead"))
  #time until relapse
  df$tdfs
  #relapse indicators
  df$deltar <- factor(df$deltar,
                      levels = c(0,1),
                      labels = c("No Relapse","Relapse"))
  
  ##disease free survival indicator
  df$deltadfs <- factor(df$deltadfs,
                        levels = c(0,1),
                        labels = c("Disease free survival","Dead or Relapsed"))
  ###time until graft disease
  df$ta
  ####acute graftversus host disease
  df$deltaa <- factor(df$deltaa,
                      levels = c(0,1),
                      labels = c("No aGVHD","aGVHD"))
  #platelet recovery time
  df$tp
  #platelet outcome
  df$deltap <- factor(df$deltap,
                      levels = c(0,1),
                      labels = c("No Platelet Recovery",
                                 "Platelet Recovery"))
  #disease group
  df$disgroup <- factor(df$disgroup,
                        levels = c(1,2,3),
                        labels = c("ALL",
                                   "Low Risk AML",
                                   "High Risk AML"))
  #Patient age
  df$age <- as.numeric(df$age)
  #male
  df$male <- factor(df$male,
                    levels = c(0,1),
                    labels = c(
                      "Female",
                      "Male"
                    ))
  #CMV
  df$cmv <- factor(df$cmv,
                   levels = c(0,1),
                   labels = c(
                     "CMV Positive",
                     "CMV Negative"
                   ))
  #donor age
  df$donorage <- as.numeric(df$donorage)
  #donormale
  df$donormale <- factor(df$donormale,
                         levels = c(0,1),
                         labels = c("Female Donor","Male Donor"))
  #donor CMV
  df$donorcmv <- factor(df$donorcmv,
                        levels = c(0,1),
                        labels = c(
                          "Donor CMV Negative",
                          "Donor CMV Positive"
                        ))
  #wait time
  df$waittime <- as.numeric(df$waittime)
  #FAB
  df$fab <- factor(df$fab,
                   levels = c(1,0),
                   labels = c(
                     "FAB grade 4 or 5 and AML",
                     "Other"
                     
                   ))
  #hospital
  df$hospital <- factor(df$hospital,
                        levels = c(1,2,3,4),
                        labels = c(
                          "The Ohio State University (Columbus, OH)",
                          "Alfred (Melbourne, Australia)",
                          "St. Vincent (Sydney, Australia)",
                          "Hahnemann (Philadelphia, PA)"
                        ))
  #mtx
  df$mtx <- factor(df$mtx,
                   levels = c(0,1),
                   labels = c("No Methotrexate",
                              "Yes Methotrexate"))
  
  
  ###making the table------
  table_1 <- tableone::CreateTableOne(
    vars = c(
      "age",
      "male",
      "cmv",
      "donorage",
      "donormale",
      "donorcmv",
      "waittime",
      "disgroup",
      "mtx",
      "hospital",
      "deltas",
      "deltar",
      "deltaa",
      "deltadfs",
      "deltap"
    ),
    strata = "fab",
    data = df, test = TRUE
  )
  
  tab1 <- print(table_1, varLabels = T, showAllLevels = T)
  
}

## Question 2
{
  library(survival)
  library(ggplot2)
  surv<-Surv(dt$tdfs, dt$deltadfs)
  km<-survfit(surv~1, conf.type="log-log")
  plot(km, conf.int=0.95, xlab="Time (days)", ylab="Survival Probability",
       main="KM Estimation of Disease-Free Survival Time")
  
  median<-quantile(km, probs=0.5)
  tab<-cbind(median$quantile, median$lower, median$upper)
  colnames(tab)<-c("est", "lower", "upper")
  tab
}

## Question 3
{
  #age
  cox_age <- coxph(surv~age,data=dt)
  summary(cox_age)
  
  #sex
  dt$male<-as.factor(dt$male)
  cox_male <- coxph(surv~male,data=dt)
  summary(cox_male)
  
  #donor age
  cox_donorage <- coxph(surv~donorage,data=dt)
  summary(cox_donorage)
  
  #donor sex
  dt$donormale<-as.factor(dt$donormale)
  cox_donormale <- coxph(surv~donormale,data=dt)
  summary(cox_donormale)
  
  #CMV
  dt$cmv<-as.factor(dt$cmv)
  cox_cmv <- coxph(surv~cmv+age+male,data=dt)
  summary(cox_cmv)
  
  #donor CMV
  dt$donorcmv<-as.factor(dt$donorcmv)
  cox_donorcmv <- coxph(surv~donorcmv+donorage+donormale,data=dt)
  summary(cox_donorcmv)
  
  #wait time from diagnosis to transplantation
  cox_wait <- coxph(surv~waittime+hospital,data=dt)
  summary(cox_wait)
  
  #disease group
  dt$disgroup<-as.factor(dt$disgroup)
  cox_disgroup <- coxph(surv~disgroup+cmv+age+male,data=dt)
  summary(cox_disgroup)
  
  #FAB
  dt$fab<-as.factor(dt$fab)
  cox_fab <- coxph(surv~fab+cmv+age+male,data=dt)
  summary(cox_fab)
  
  #prophylactic use of methotrexate
  dt$mtx<-as.factor(dt$mtx)
  cox_mtx <- coxph(surv~mtx+cmv+age+male,data=dt)
  summary(cox_mtx)
}

## Question 4
{
  dt$sex_match<-ifelse((dt$male==1 & dt$donormale==1)|(dt$male==0 & dt$donormale==0), 1, 0)
  
  surv2<-Surv(dt$tdfs, dt$deltar)
  cox_deltaa2 <- coxph(surv2~deltaa+age+sex_match+cmv+mtx,data=dt)
  summary(cox_deltaa2)
  
  cox_deltaa <- coxph(surv~deltaa+age+sex_match+cmv+mtx,data=dt)
  summary(cox_deltaa)
  
  #AFT models
  surv=with(dt, Surv(tdfs, deltadfs))
  
  gengamma_aft <- flexsurvreg(surv~deltaa+age+sex_match+cmv+mtx,data=dt, dist="gengamma")
  gengamma_aft
  
  gengamma_aft.res <- gengamma_aft$res
  gengamma_aft.wald <- gengamma_aft.res[,1]/gengamma_aft.res[,4]
  gengamma_aft.p <- 2*pnorm(-abs(gengamma_aft.wald))
  gengamma_aft.p
}

## Question 5
{
  aDVHD<-subset(dt, deltaa==1)
  surv3<-Surv(aDVHD$tdfs, aDVHD$deltadfs)
  
  #age
  cox_age <- coxph(surv3~age,data=aDVHD)
  summary(cox_age)
  
  #sex
  aDVHD$male<-as.factor(aDVHD$male)
  cox_male <- coxph(surv3~male,data=aDVHD)
  summary(cox_male)
  
  #donor age 
  cox_donorage <- coxph(surv3~donorage,data=aDVHD)
  summary(cox_donorage)
  
  #donor sex
  aDVHD$donormale<-as.factor(aDVHD$donormale)
  cox_donormale <- coxph(surv3~donormale,data=aDVHD)
  summary(cox_donormale)
  
  #CMV
  aDVHD$cmv<-as.factor(aDVHD$cmv)
  cox_cmv <- coxph(surv3~cmv+age+male,data=aDVHD)
  summary(cox_cmv)
  
  #donor CMV
  aDVHD$donorcmv<-as.factor(aDVHD$donorcmv)
  cox_donorcmv <- coxph(surv3~donorcmv+donorage+donormale,data=aDVHD)
  summary(cox_donorcmv)
  
  #wait time from diagnosis to transplantation
  cox_wait <- coxph(surv3~waittime+hospital,data=aDVHD)
  summary(cox_wait)
  
  #disease group
  aDVHD$disgroup<-as.factor(aDVHD$disgroup)
  cox_disgroup1 <- coxph(surv3~cmv+age+male,data=aDVHD)
  cox_disgroup2 <- coxph(surv3~disgroup+cmv+age+male,data=aDVHD)
  anova(cox_disgroup1, cox_disgroup2)
  
  #FAB
  aDVHD$fab<-as.factor(aDVHD$fab)
  cox_fab <- coxph(surv3~fab+cmv+age+male,data=aDVHD)
  summary(cox_fab)
  
  #prophylactic use of methotrexate
  aDVHD$mtx<-as.factor(aDVHD$mtx)
  cox_mtx <- coxph(surv3~mtx+cmv+age+male,data=aDVHD)
  summary(cox_mtx)
  
}

## Question 6
{
  s.dt = with(dt, Surv(ta, deltaa))
  cox <- coxph(s.dt ~ as.factor(mtx), data = dt)
  cox_confounders <- coxph(s.dt ~ as.factor(mtx) + as.factor(fab) + as.factor(disgroup) + age, data = dt)
  
  surv<-Surv(dt$ta, dt$deltaa)
  surv.nomtx = with(dt[mtx == 0,], Surv(ta, deltaa))
  surv.mtx = with(dt[mtx == 1,], Surv(ta, deltaa))
  
  km_nomtx<-survfit(surv.nomtx~1, conf.type="log-log", data = dt[mtx == 0])
  km_mtx<-survfit(surv.mtx~1, conf.type="log-log", data = dt[mtx == 1])
  plot(km_nomtx,  xlab="Time (days)", ylab="Probability of no aGVHD", xlim = c(0,100), main = 'KM Estimation of Time until aGVHD')
  lines(km_mtx, col = 'red')
  legend('bottomright', col = c('black', 'red'), legend = c('No methotrexate', 'Methotrexate'), lty = 1)
}

## Question 7
{
  coxfit = coxph(surv~deltap+age+male+disgroup+mtx,data=dt)
  summary(coxfit)
  
  gengamma_aft <- flexsurvreg(surv~deltap+age+male+disgroup+mtx,data=dt, dist="gengamma")
  gengamma_aft
  
  gengamma_aft.res <- gengamma_aft$res
  gengamma_aft.wald <- gengamma_aft.res[,1]/gengamma_aft.res[,4]
  gengamma_aft.p <- 2*pnorm(-abs(gengamma_aft.wald))
  gengamma_aft.p
  
  surv_relapse<-Surv(dt$tdfs, dt$deltar)
  coxfit2 = coxph(surv_relapse~deltap+age+male+disgroup+mtx,data=dt)
  summary(coxfit2)
  
  
}

#####################################################
#Tables and figures
#####################################################
## Table 1
{
  ###cleaning variables-----
  df$X <- NULL
  #survival time
  df$ts 
  #death indicator
  df$deltas <- factor(df$deltas,
                      levels = c(0,1),
                      labels = c("Alive","Dead"))
  #time until relapse
  df$tdfs
  #relapse indicators
  df$deltar <- factor(df$deltar,
                      levels = c(0,1),
                      labels = c("No Relapse","Relapse"))
  
  ##disease free survival indicator
  df$deltadfs <- factor(df$deltadfs,
                        levels = c(0,1),
                        labels = c("Disease free survival","Dead or Relapsed"))
  ###time until graft disease
  df$ta
  ####acute graftversus host disease
  df$deltaa <- factor(df$deltaa,
                      levels = c(0,1),
                      labels = c("No aGVHD","aGVHD"))
  #platelet recovery time
  df$tp
  #platelet outcome
  df$deltap <- factor(df$deltap,
                      levels = c(0,1),
                      labels = c("No Platelet Recovery",
                                 "Platelet Recovery"))
  #disease group
  df$disgroup <- factor(df$disgroup,
                        levels = c(1,2,3),
                        labels = c("ALL",
                                   "Low Risk AML",
                                   "High Risk AML"))
  #Patient age
  df$age <- as.numeric(df$age)
  #male
  df$male <- factor(df$male,
                    levels = c(0,1),
                    labels = c(
                      "Female",
                      "Male"
                    ))
  #CMV
  df$cmv <- factor(df$cmv,
                   levels = c(0,1),
                   labels = c(
                     "CMV Positive",
                     "CMV Negative"
                   ))
  #donor age
  df$donorage <- as.numeric(df$donorage)
  #donormale
  df$donormale <- factor(df$donormale,
                         levels = c(0,1),
                         labels = c("Female Donor","Male Donor"))
  #donor CMV
  df$donorcmv <- factor(df$donorcmv,
                        levels = c(0,1),
                        labels = c(
                          "Donor CMV Negative",
                          "Donor CMV Positive"
                        ))
  #wait time
  df$waittime <- as.numeric(df$waittime)
  #FAB
  df$fab <- factor(df$fab,
                   levels = c(1,0),
                   labels = c(
                     "FAB grade 4 or 5 and AML",
                     "Other"
                     
                   ))
  #hospital
  df$hospital <- factor(df$hospital,
                        levels = c(1,2,3,4),
                        labels = c(
                          "The Ohio State University (Columbus, OH)",
                          "Alfred (Melbourne, Australia)",
                          "St. Vincent (Sydney, Australia)",
                          "Hahnemann (Philadelphia, PA)"
                        ))
  #mtx
  df$mtx <- factor(df$mtx,
                   levels = c(0,1),
                   labels = c("No Methotrexate",
                              "Yes Methotrexate"))
  
  
  ###making the table------
  table_1 <- tableone::CreateTableOne(
    vars = c(
      "age",
      "male",
      "cmv",
      "donorage",
      "donormale",
      "donorcmv",
      "waittime",
      "disgroup",
      "mtx",
      "hospital",
      "deltas",
      "deltar",
      "deltaa",
      "deltadfs",
      "deltap"
    ),
    strata = "fab",
    data = df, test = TRUE
  )
  
  tab1 <- print(table_1, varLabels = T, showAllLevels = T)
  
}

## Figure 1
{
  surv<-Surv(dt$tdfs, dt$deltadfs)
  km<-survfit(surv~1, conf.type="log-log")
  plot(km, conf.int=0.95, xlab="Time (days)", ylab="Survival Probability",
       main="KM Estimation of Disease-Free Survival Time")
  
}

## Figure 2
{
  surv<-Surv(dt$ta, dt$deltaa)
  surv.nomtx = with(dt[mtx == 0,], Surv(ta, deltaa))
  surv.mtx = with(dt[mtx == 1,], Surv(ta, deltaa))
  
  km_nomtx<-survfit(surv.nomtx~1, conf.type="log-log", data = dt[mtx == 0])
  km_mtx<-survfit(surv.mtx~1, conf.type="log-log", data = dt[mtx == 1])
  plot(km_nomtx,  xlab="Time (days)", ylab="Probability of no aGVHD", xlim = c(0,100), main = 'KM Estimation of Time until aGVHD')
  lines(km_mtx, col = 'red')
  legend('bottomright', col = c('black', 'red'), legend = c('No methotrexate', 'Methotrexate'), lty = 1)
}