####################### Gems package: Vignette ############################

#Link to vignette: https://cran.r-project.org/web/packages/gems/vignettes/gems.pdf

#Install gems package
install.packages("gems")
library(gems)

###################### Example in vignette (p. 14) ########################

#The tavi data set contains data on kidney injuries, bleeding complications and the combined
#endpoint of stroke or death for 194 patients. The variables kidney, bleeding, death are in-
#dicator variables that show if an event has occurred; the variables kidney.dur, bleeding.dur,
#death.dur are the times at which the events occurred or the patients were censored.

#Load data
data("tavi", package = "gems")
#Preview first few rows in data
head(tavi)


#DAG depicted in Figure 6 is assumed (pg. 15). Since no patients
#experience both kidney injury and bleeding complications, we assume these events to be
#mutually exclusive.

#Create Transition matrix using mstate package according to the DAG
  #Install mstate package
  install.packages("mstate")
  library(mstate)
  
  #Transition Matrix
  tmat <- transMat(x = list(c(2, 3, 4), c(4), c(4), c()),
                   names = c("TAVI", "Kidney Injury", "Bleeding", "Stroke/Death"))
  tmat
  
  
#To estimate the transition-specific hazard functions use "msprep" to 
#get the data into long format, and "split" to split the data
#according to the transition. 
  
  #msprep -> data to long format
   mstavi <- msprep(data = tavi, trans = tmat,
                   time = c(NA, "kidney.dur", "bleeding.dur", "death.dur"),
                   status = c(NA, "kidney", "bleeding", "death"))
   head(mstavi)
  
  
  #split -> split data according to transition
   mstavi$time[mstavi$time == 0] <- .Machine$double.eps
   msplit <- split(mstavi, mstavi$trans)
   head(msplit[[5]])
  

#Now, we fit an exponential distribution to all transition times. For each transition,
#we estimate the rate and the variance.
   
   exp.fit <- sapply(msplit, function(x) summary(survreg(Surv(time,
                                                              status) ~ 1, data = x, dist = "exponential")))
   exp.coef <- unlist(exp.fit["coefficients", ])
   exp.var <- unlist(exp.fit["var", ])

#Next, Simulate the model and compare simulated mortality to a Kaplan-Meier graph of Mortality.

   states <- 4 #specify 4 states/transitions in the model
   maxtime <- max(mstavi$time)
   ind <- which(!is.na(tmat), arr.ind = TRUE)
   
   #Generate Hazard Matrix
   hm <- generateHazardMatrix(states)
      for (i in 1:dim(ind)[1]) {
      hm[[ind[i, 1], ind[i, 2]]] <- "Weibull"
     }
   
   #Generate Parameter Matrix
   par <- generateParameterMatrix(hm)
       for (i in 1:dim(ind)[1]) {
       par[[ind[i, 1], ind[i, 2]]] <- list(shape = 1, scale = exp(exp.coef[i]))
       }
   
   #Generate Covariance Matrix
   cov <- generateParameterCovarianceMatrix(par)
       for (i in 1:dim(ind)[1]) {
       cov[[ind[i, 1], ind[i, 2]]] <- matrix(c(0, 0, 0, exp.var[i]), nrow = 2)
       }
   
   #Simulate Cohort: include matrices created above and cohort size
   ds <- simulateCohort(transitionFunctions = hm, parameters = par,
                        cohortSize = 100 * nrow(tavi), parameterCovariances = cov,
                        to = maxtime)
   cinc <- cumulativeIncidence(ds, 0:maxtime, colnames(tmat), M = 100)
     
  
#Split time into monthly intervals and calculate piecewise constant hazard 
#functions.
#Use pehaz function in muhaz package.
   
  install.packages("muhaz")   
  library(muhaz)

  timeStep <- 30
  pwexp <- sapply(msplit, function(x) pehaz(x$time, x$status,
                                             width = timeStep, min.time = 0, max.time = max(mstavi$time)))
  cuts <- pwexp["Cuts", ]
  pwhazard <- pwexp["Hazard", ]
  
  
#Re-parameterise the hazard functions with piecewise constant hazards and simulate again.

    #Generate Hazard Matrix
     hm2 <- generateHazardMatrix(states)
     
     for (i in 1:dim(ind)[1]) {
         hm2[[ind[i, 1], ind[i, 2]]] <- function(t, rates)
         rates[t / timeStep + 1]
          }
     
     #Generate Parameter Matrix
      par2 <- generateParameterMatrix(hm2)
      
      for (i in 1:dim(ind)[1]) {
       par2[[ind[i, 1], ind[i, 2]]] <- list(rates = pwhazard[[i]])
      }
      
      #Simulate Cohort
      ds2 <- simulateCohort(transitionFunctions = hm2, 
                            parameters = par2, cohortSize = 100 * nrow(tavi), to = maxtime)
      cinc2 <- cumulativeIncidence(ds2, 0:maxtime, colnames(tmat), M = 100)
      
      