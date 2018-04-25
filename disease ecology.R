library(deSolve)
library(ggplot2)
library(reshape2)
#Function which contains our model
SIRfunc = function(t, x, parameters){
  #x contains each of our population groups at time t:
  S1 = x[1]
  S2 = x[2]
  I1 = x[3]
  I2 = x[4]
  R1 = x[5]
  R2 = x[6]
  D  = x[7]
  
  #Calculate each of the changes in subgroup population, then return them as a list
  with(as.list(parameters),{
    N = S1+S2+I1+I2+R1+R2+D     
    
    dS1 = b*N -beta*(I1+I2)*S1/N - delta*S1 #healthy suspt change
    dS2 = b*N -beta*(I1+I2)*S2/N - delta*S2 #immuno suspt change
    
    dI1 = beta*(I1+I2)*S1/N - gamma1*I1 - alpha1*I1 #healthy infect change
    dI2 = beta*(I1+I2)*S2/N - gamma2*I2 - alpha2*I2 #immuno infect change
    
    dR1 = gamma1*I1 - delta*R1 #healthy recov change
    dR2 = gamma2*I2 - delta*R1 #immuno recov change
    
    dD = alpha1*I1+alpha2*I2 #death change
    
    out = c(dS1,dS2,dI1,dI2,dR1,dR2,dD)
    list(out)
  })
}
#Initial Parameters
time = c(0:150)
b =0.00024 #birth rate  #.00003426
beta = .41 #infectivity
delta = .0001589 #general death rate
alpha1 = .138/365 #healthy infect death rate
alpha2 = .26/365 #immuno infect death rate
gamma1 = .2959 #healthy recov rate
gamma2 = .2 #immuno recov rate

#Starting Values
N = 10000
percent = .5 #percent immunocomp Change value here
I1_0 = 1
I2_0 = 0
R1_0 = 0
R2_0 = 0
S1_0 = N-percent*N - I1_0
S2_0 = percent*N
D_0 = 0

parameters = c(alpha1 = alpha1, alpha2 = alpha2, beta = beta, gamma1 = gamma1, gamma2 = gamma2)
initial = c(S1 = S1_0, S2 = S2_0, I1 = I1_0,I2 = I2_0, R1 = R1_0, R2 = R2_0, D = D_0)
model.out = as.data.frame(lsoda(initial, time, SIRfunc, parameters))
model.out$Itotal <- model.out$I1 + model.out$I2
head(model.out)

modelToPlot <- melt(model.out[1:8], id.vars = c("time"), variable.name = "Category", value.name = "n")
ggplot(modelToPlot, aes(x = time, y = n, color = Category)) +theme(text = element_text(size=15)) + geom_line() + ggtitle("Disease Progression Over Time") + scale_color_discrete(labels = c("Susceptible", "Susceptible-IC", "Infected","Infected-IC", "Recovered","Recovered-IC", "Dead"))
final = model.out[150,]

print(paste("Susceptible: ", (final$S1/N*100), "%"))
print(paste("Susceptible I: ", (final$S2/N*100), "%"))
print(paste("Infected: ", (final$I1/N*100), "%"))
print(paste("Infected I: ", (final$I2/N*100), "%"))
print(paste("Recovered", (final$R1/N*100),"%"))
print(paste("Recovered", (final$R2/N*100),"%"))

print(paste("Max infected: ", round(max(model.out$I1+model.out$I2)/N*100),"%"))

maxVal <- max(model.out$Itotal)/N*100

index <- which(model.out$Itotal == max(model.out$Itotal))
print(paste("Time Max infected: ", model.out$time[[53]]," days"))

print(paste("Dead: ", (final$D/N*100),"%"))