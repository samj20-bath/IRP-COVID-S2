# SEIR stochastic model (almost an "Standard SEIR model" see Britton p.11)
rm(list=ls())

setwd("~/Desktop/SAMBa/Semester 1/IRP/R Code/IRP-COVID-S2/Data")
data<-read.csv("Data IRP Covid19 Strand 2 - Hoja 3.csv",header=T,stringsAsFactors = F)
data$Admissions<-as.numeric(data$Admissions)


HOSP<-matrix(NA,nrow=1005,ncol=100)
SUS<-matrix(NA,nrow=1001,ncol=100)
INF<-matrix(NA,nrow=1001,ncol=100)
REC<-matrix(NA,nrow=1001,ncol=100)


a1vals<-seq(0.1,0.3,length.out = 20)
a2vals<-seq(0.1,0.5,length.out = 20)
  
# AMATRIX<-matrix(NA,ncol=20,nrow=20)
# 
# for(a1 in 1:length(a1vals)){
# for(a2 in a1:length(a2vals)){
#   
for(reps in 1:100){
N=56000000 # Population size
m<- 245 # Number of initial infectives (t0 = 12 march)

beta_r <- 0.86# Constant transmission rate by day R_0=3.4444 (According to Liu) (R_0/avg.days)
rho <-0.33333 # 1/mean of incubation period
gamma<-0.25 # 1/mean of infectious period (4 days according to Britton)

theta<-5 # Maximium likelihood estimator
p_H <- 0.06423 #  MLE Estimated by me

t1<-11 #  days to implement first lockdown
t2 <- 103 # 4 july, NO LOCKDOWN
t3 <- 238 # 6 nov, LOCKDOWN AGAIN
t4 <- 266 # 3 dec, NO LOCKDOWN AGAIN
t5 <- 298 # 4 ene, LOCKDOWN AGAIN

# We're assuming the estimated rates are in days units

# Values of impact of lockdown and postlockdown (Obtained by least squares)
#alpha1 <-0.2797965116 
#alpha2 <-0.3136474908

# alpha1<-a1vals[a1]
# alpha2<-a2vals[a2]
  
alpha1<-0.2443221
alpha2<-0.3206636

beta<-function(t){ #t0 is the date of implementation
    ifelse(t < t1, beta_r,
    ifelse(t >= t1 & t < t2, beta_r*alpha1,jajaja
    ifelse(t >= t2 & t < t3, beta_r*alpha2,
    ifelse(t >= t3 & t < t4 , beta_r*alpha1,
    ifelse(t >= t4 & t < t5 , beta_r*alpha2,
           ifelse(t>=t5, beta_r*alpha1, beta_r))))))
}

# Jumping probabilities
P_ <- function(t,I){
  1-exp(-beta(t)*I/N)
}
pc <- 1-exp(-rho)
pr <- 1-exp(-gamma)

# Series of compartments and jumps
S<-c()
B<-c();E<-c() # B is the jump; E the no. of exposed
C<-c();I<-c() 
D<-c();R<-c()
H<-c() # Hospitalized

# Initialization of compartments
S[1]<-N-m
E[1]<-0
I[1]<-m
R[1]<-0
H[1:theta]<-0


# Simulation
for(i in 1:1000){
  B[i]<-rbinom(1,S[i],P_(i,I[i]))
  C[i]<-rbinom(1,E[i],pc)
  D[i]<-rbinom(1,I[i],pr)
  H[i+theta]<-rbinom(1,I[i],p=p_H)
  #******************************
  S[i+1]<-S[i]-B[i]
  E[i+1]<-E[i]+B[i]-C[i]
  I[i+1]<-I[i]+C[i]-D[i]
  R[i+1]<-N-S[i+1]-E[i+1]-I[i+1]
}

I2=I
I=I/N
S=S/N
R=R/N
E=E/N


SUS[,reps]<-S
INF[,reps]<-I2
REC[,reps]<-R
HOSP[,reps]<-H
} # End of repetitions for simulating
#     AMATRIX[a1,a2]<-mean((apply(HOSP,1,mean)[4:330]-data$Admissions[4:330])^2)
#   }
# }
# ************** P L O T S ***************

# Infected people ----------------
quartz(width = 11,height=4.4)
par(mfrow=c(1,3))
plot(log(apply(INF,1,mean),10),type="l",col="red",lwd=2,xlim=c(0,500),
     main="Log of Infected people",ylab="log-10 of Infected",xlab="")
lines(log(apply(INF,1,mean)-2*apply(INF,1,sd),10),type="l",col="gray")
lines(log(apply(INF,1,mean)+2*apply(INF,1,sd),10),type="l",col="gray")
abline(v=11,lty=2)
abline(v=103,lty=2)
abline(v=238,lty=2)
abline(v=266,lty=2)
abline(h=0)

# Hospitalisations ----------------
plot(apply(HOSP,1,mean),type="l",main="Hospitalisations",ylab="Hospitalised people",xlab="Days from March 12",
     xlim=c(0,500),col="forestgreen",lwd=2)#,ylim=c(0,5000))
lines(apply(HOSP,1,mean)-2*apply(HOSP,1,sd),type="l",col="gray")
lines(apply(HOSP,1,mean)+2*apply(HOSP,1,sd),type="l",col="gray")
abline(v=11,lty=2)
abline(v=103,lty=2)
abline(v=238,lty=2)
abline(v=266,lty=2)
abline(h=0)
points(data$days.from.12.March,data$Admissions,type="l",col="red")
legend("topright",lty=c(1,1),lwd=c(2,2),col=c("forestgreen","red"),
       legend=c("Simulated","Observed"),cex = 0.8,bg = "white")

# print<-sum(I[1:339])
# ***********S I R series in one plot ---------
# quartz()
# par(mfrow=c(1,1))
IN<-INF/N

plot(apply(IN,1,mean),type="l",xlim=c(0,500),ylim=c(0,1),col="red",ylab="People",lwd=2,xlab="",
     main="COVID-19 Stochastic SEIR model")
lines(apply(IN,1,mean)-2*apply(IN,1,sd),lwd=1,col="gray")
lines(apply(IN,1,mean)+2*apply(IN,1,sd),lwd=1,col="gray")
lines(apply(SUS,1,mean),col="blue",lwd=2)
lines(apply(REC,1,mean),col="purple",lwd=2)
abline(v=11,lty=2)
abline(v=103,lty=2)
abline(v=238,lty=2)
abline(v=266,lty=2)
abline(h=0)
legend("right",lty=c(1,1,1),lwd=c(2,2,2),col=c("red","blue","purple","n"),
       legend=c("Infected","Susceptible","Recovered"),cex = 0.8,bg = "white")

sum(apply(INF,1,mean)[1:339])

#************** Assessment of simulation *********************

