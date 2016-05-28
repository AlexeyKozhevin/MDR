source("simulate.R")
source("MDR-EFE.R")
source("distribution.R")
source("powerplot.R")

#the following parameters are fixed for all experiments

r = 3
n = 100
significant_factors = c(2,3,5)
model_probabilities <- c(0.05,0.1,0.2)/2
ratio <- 0.1
MC <- 100

#the following parameters take different values. 
#For each combinations of them we should simulate MC = 100 experiments

costs <- (1:5)*100
models <- 1:3

#auxiliary variables for names in csv table

signif_name <- paste0("pred",1:r)
ind_name <- paste0("ind",1:r)
fix_name <- paste0("fix",1:r)

#TMR <- rep(0,5)

for (model in models)
{
  for (cost in costs)
  {
    
    #for each combination of variables model and cost create data.frame 
    #in order to save results of the experiments 
    
    results <- matrix(0,ncol=10,nrow = MC)
    results <- data.frame(results)
    colnames(results)<-c("model","total_cost","Nstr","tildeN1","tildeN2","powerindkn","powerindun",
                         "powerstr1","powerstr2kn","powerstr2un")
    
    #compute some parameter for MDR which depends on model and cost
    
    
    pr_of_case <- model_probabilities[model]
    N_fix <- sstr(cost,ratio,pr_of_case)[1]
    N_fix <- 2*(N_fix%/%2)

    #start experiments
    for (rep in 1:MC)
    {   
      #simulate dataset
      
      time <- Sys.time()
      data <- simulated_SNPs(n,N_fix,cost,significant_factors,model,type="comb")
      print(Sys.time() - time)
      
      #iMDR, probability of case is known
      
      time <- Sys.time() 
      indep.model.knownp <- mdr.efe(data[[1]],r,type = "indep",prob = pr_of_case,parallel = T)
      print(Sys.time() - time)
      
      #iMDR, probability of case is unknown
      
      time <- Sys.time() 
      indep.model.unknownp <- mdr.efe(data[[1]],r,type = "indep",parallel = T)
      print(Sys.time() - time)
      
      #sMDR, w = 0.1, probability of case is known
      
      time <- Sys.time() 
      fixed.model1 <- mdr.efe(data[[2]],r,type = "fix",prob = pr_of_case,parallel = T)
      print(Sys.time() - time)
      
      #sMDR, w = 0, probability of case is known
      
      time <- Sys.time() 
      fixed.model2.knownp <- mdr.efe(data[[3]],r,type = "fix",prob = pr_of_case,parallel = T)
      print(Sys.time() - time)
      
      #sMDR, w = 0, probability of case is unknown
      
      time <- Sys.time() 
      fixed.model2.unknownp <- mdr.efe(data[[3]],r,type = "fix",prob = data[[4]],parallel = T)
      print(Sys.time() - time)
      
      #compute results of the experiments
      
      powerindkn <- power(significant_factors,indep.model.knownp)
      powerindun <- power(significant_factors,indep.model.unknownp)
      powerstr1 <- power(significant_factors,fixed.model1)
      powerstr2kn <- power(significant_factors,fixed.model2.knownp)
      powerstr2un <- power(significant_factors,fixed.model2.unknownp)
      
      #log is a result which I should know for each experiment 
      
      log <- c(model,cost,N_fix,data[[5]],data[[6]],powerindkn,powerindun,powerstr1,powerstr2kn,powerstr2un)
      
      TMR <- TMR + c(powerindkn,powerindun,powerstr1,powerstr2kn,powerstr2un)
      
      results[rep,] <- log
      print(c(rep,log))
      print(TMR/rep)
      
      #each experiment save results as one string in .csv
      filename <- paste0("./csv/tmp/model",model,"cost",cost,"_",format(Sys.time(), "%m%d%H%M%S"),"_tmp.csv")
      write.csv(results[rep,],filename)
    }    
    #all experiments for fixed cost and model in one .csv
    filename <- paste0("./csv/model",model,"cost",cost,"_",format(Sys.time(), "%m%d%H%M%S"),"_final.csv")
    write.csv(results,filename)
  }
}
