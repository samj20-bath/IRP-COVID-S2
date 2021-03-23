library(plotly)
library(readxl)

setwd("~/Desktop/SAMBa/Semester 1/IRP/Data")
data<-read.csv("Data IRP Covid19 Strand 2 - Hoja 3.csv",header=T,stringsAsFactors = F)
three_months_days<-data[data$YEAR=="2020"&data$MONTH%in%c(8,9,10),]$days.from.12.March
data$Admissions<-as.numeric(data$Admissions)

t1 <-11  # LOCKDOWN
t2 <-103 # 4 july, NO LOCKDOWN
t3 <- 238 # 6 nov, LOCKDOWN AGAIN
t4 <- 266 # 3 dec, NO LOCKDOWN AGAIN
t5 <- 298 # 4 ene, LOCKDOWN AGAIN

#***************************** P L O T S ********************************
quartz(height = 4,width = 10)
par(mfrow=c(1,2))
plot(data$days.from.12.March,data$newCasesByPublishDate,type="l",xlab="Days from March 12",
     ylab="Infected people",col="brown3",lwd=2)
grid(nx = NULL, ny = NULL, col = "azure2", lty = "dotted")
abline(v=t0,lty=4,lwd=1.5)
abline(v=t1,lty=4,lwd=1.5)
abline(v=t2,lty=4,lwd=1.5)
abline(v=t3,lty=4,lwd=1.5)
abline(v=t4,lty=4,lwd=1.5)
abline(v=271,lty=4,col="forestgreen",lwd=1.5)
abline(v=281,lty=4,col="blue",lwd=1.5)

text=c("Change in lockdown policy",
       "New strain exceeds 50%","Vaccination starts")
legend("topleft",legend=text,
       text.width = strwidth(text)[1]/2,col=c("black","blue","forestgreen"),lty=c(4,4,4),lwd=c(1.5,1.5,1.5),cex=0.6,bg="white")

plot(data$days.from.12.March,data$Admissions,type="l",xlab="Days from March 12",
     ylab="Hospitalisations",col="orange",lwd=2)
grid(nx = NULL, ny = NULL, col = "azure2", lty = "dotted")
abline(v=t0,lty=4,lwd=1.5)
abline(v=t1,lty=4,lwd=1.5)
abline(v=t2,lty=4,lwd=1.5)
abline(v=t3,lty=4,lwd=1.5)
abline(v=t4,lty=4,lwd=1.5)
abline(v=271,lty=4,col="forestgreen",lwd=1.5)
abline(v=281,lty=4,col="blue",lwd=1.5)
text=c("Change in lockdown policy",
       "New strain exceeds 50%","Vaccination starts")
legend("topleft",legend=text,
       text.width = strwidth(text)[1]/2,col=c("black","blue","forestgreen"),
       lty=c(4,4,4),lwd=c(1.5,1.5,1.5),cex=0.6,bg="white")
title("Cases by reported date and hospitalisations", line=-3 ,outer = TRUE)

#**************************************************************************************
data_google<-read_xlsx("google_activity_by_London_Borough.xlsx",sheet=2)

quartz(height = 5,width = 9)
plot(-20:366,data_google$`Retail and recreation`,type="l",col="purple",
     xlab="Days from March 12",ylab="Change in mobility index(%)",main="Google mobility index")
lines(-20:366,data_google$`Grocery and pharmacy`,type="l",col="darkturquoise")
abline(v=t0,lty=4,lwd=1.5)
abline(v=t1,lty=4,lwd=1.5)
abline(v=t2,lty=4,lwd=1.5)
abline(v=t3,lty=4,lwd=1.5)
abline(v=t4,lty=4,lwd=1.5)
legend("topright",legend=c("Retail and recreation","Grocery and pharmacy","Change in lockdown policy"),col=c("purple","darkturquoise","black"),
       lty=c(1,1,2),bg="white")


#******************* Estimation of alphas according to google ***********************************
gugl<-data.frame(time=-5:366,data_google$`Retail and recreation`[-c(1:15)])
gugl$rr<-gugl$data_google..Retail.and.recreation.-mean(gugl$data_google..Retail.and.recreation.)
gugl$rr<-gugl$rr+abs(min(gugl$rr))
gugl$rr<-gugl$rr/max(gugl$rr)
gugl$rr<-gugl$rr*(0.86-0.2034884)+0.2034884
mean(gugl[gugl$time>=t1&gugl$time<t2|(gugl$time>=t3&gugl$time<t4),]$rr)
mean(gugl[gugl$time>=t2&gugl$time<t3|(gugl$time>=t4&gugl$time<t4),]$rr)

#***********  Estimation of theta and p_h for hospitalisation data ********************

ll<-function(theta,h){
  sum<-0
  hospit<-data[data$days.from.12.March%in%c(three_months_days+h),]$Admissions
  newCases<-data[data$days.from.12.March%in%c(three_months_days),]$newCasesByPublishDate
  for(i in 1:length(hospit)){
  sum<-sum+dbinom(hospit[i],newCases[i],theta,log=T)
  }
  return(-log(-sum))
}

sum(data[data$days.from.12.March%in%c(three_months_days+5),]$Admissions)/sum(newCases<-data[data$days.from.12.March%in%c(three_months_days),]$newCasesByPublishDate)

loglik_h<-c()
for(i in 1:50){
  hospit<-data[data$days.from.12.March%in%c(three_months_days+i),]$Admissions
  newCases<-data[data$days.from.12.March%in%c(three_months_days),]$newCasesByPublishDate
  loglik_h[i]<-print(ll(sum(hospit)/sum(newCases),i))
}

plot(1:50,loglik_h,type="l") # h=5 is the value with greatest loglikelihood

h<-1:10
theta<-seq(from=0.03,to=0.1,length.out=100)
x<-c()
hache<-c()
for(i in 1:length(h)){
  for(j in 1:100){
    x<-c(x,ll(theta[j],h[i]))
    hache<-c(hache,h[i])
  }
}

df<-data.frame(theta=rep(theta,length(h)),ll=x,h=hache)
df$h<-factor(df$h)
rownames(df)<-paste0("a",df$h,1:nrow(df))


fig <- plot_ly(df, x = ~theta, y = ~h, z = ~ll, 
               type = "scatter3d", mode = "lines",color=~h,line=list(width=4)) 
fig <- fig %>% layout(
  title = "Log-likelihood* of (θ,pₕ)",
  scene = list(
    xaxis = list(title = "pₕ"),
    yaxis = list(title = "θ"),
    zaxis = list(title = "ll"))
  )

# fig

# ***PLOT**** Effective reproduction numberand L-S Estimation of alphas

R_data<-read.csv("R-and-growth-rate-time-series-26-Feb-2021.csv",sep=";",header=TRUE)
# 
# beta<-function(t,alpha1=3.5/15,alpha2=5.5/15){ #t0 is the date of implementation
#   ifelse(t < t0, beta_r,
#          ifelse(t >= t0 & t < t1, beta_r*alpha1,
#                 ifelse(t >= t1 & t < t2, beta_r*alpha2,
#                        ifelse(t >= t2 & t < t3 , beta_r*alpha1,
#                               ifelse(t >= t3 & t < t4 , beta_r*alpha2,
#                                      ifelse(t>=t4, beta_r*alpha1, beta_r))))))
# }
# plot(1:353,sapply(1:353,function(x) beta(x)/gamma),type="l",ylab="Rₒ(t)",xlab="Days from March 12",
#      main="Time dependent effective reproductive number",ylim=c(-.2,3.6))
# abline(v = 78)
# abline(v=281,col="green")
# points(R_data$days,R_data$LL,type="o",lty=2,col="blue")
# points(R_data$days,R_data$UP,type="o",lty=2,col="blue")
# legend("topright",legend=c("R estimate interval","New strain exceeds 50%"),col=c("blue","green"),lty=c(2,1),cex=0.7)

# Least squares estimation of alpha ***********************************************

alpha1_<-mean(R_data[R_data$days>=t1&R_data$days<t2|(R_data$days>=t3&R_data$days<t4),]$England)
alpha2_<-mean(R_data[R_data$days>=t2&R_data$days<t3|(R_data$days>=t4&R_data$days<t5),]$England)

beta_r <- 0.86# Constant transmission rate by day R_0=3.4444 (According to Liu) (R_0/avg.days)
rho <-0.33333 # 1/mean of incubation period
gamma<-0.25 # 1/mean of infectious period (4 days according to Britton)

alpha1<-(gamma*alpha1_)/beta_r

alpha2<-(gamma*alpha2_)/beta_r
print(alpha1)
print(alpha2)
# 0.2797965116 0.3136474908
