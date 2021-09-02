#' Fits a NULL model for EmpSPA
#'
#' Fits a null logistic regression model for binary phenotype or a null Cox proportional hazards model for time-to-event phenotype and then calculates the raw residuals or martingale residuals and the empirical cumulant generation function (CGF) of the raw residuals or martingale residuals,
#' or use residuals from a null generalized linear model to calculate its the empirical CGF.
#' @param traits a character value corresponding to phenotype. It should be "binary" for binary phenotype and "survival" for time-to-event phenotype.
#' @param formula a formula to be passed to function glm() or coxph(). For more details, please refer to package survival.
#' @param data a data.frame in which to interpret the variables named in the formula
#' @param pIDs a character vector of subject IDs. NOTE: its order should be the same as the subjects order in the formula.
#' @param gIDs a character vector of subject IDs. NOTE: its order should be the same as the subjects order of the Geno.mtx (i.e. the input of the function EmpSPA()).
#' @param range a two-element numeric vector (default: c(-100,100)) to specify the domain of the empirical CGF.
#' @param length.out a positive integer (default: 9999) for empirical CGF. Larger length.out corresponds to longer calculation time and more accurate estimated empirical CGF.
#' @param Cova.mtx a covariate matrix including a column of 1s (only needed when the type of traits is neither survival nor binary).
#' @param resid a vector of residuals from a generalized linear regression model in which Cova.mtx is the covariate matrix (only needed when the type of traits is neither survival nor binary).
#' @param ... Other arguments passed to function glm() or coxph(). For more details, please refer to package survival.
#' @return an object with a class of "binary_Null_Model" for binary phenotype or "survival_Null_Model" for time-to-event phenotype.
#' @examples
#' Please check help(EmpSPA) for a simulated example.
#' @export
#' @import survival
EmpSPA_Null_Model = function(traits="survival/binary/others",
                             formula=NULL,
                             data=NULL,
                             pIDs=NULL,
                             gIDs=NULL,
                             range=c(-100,100),
                             length.out = 10000,
                             Cova.mtx=NULL,
                             resid=NULL,
                             ...)

{
  if(traits=="survival"){
    Call = match.call()
    ### Fit a Cox model
    obj.coxph = coxph(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.coxph, range)

    ### Get the covariate matrix to adjust for genotype
    mresid = obj.coxph$residuals
    Cova = obj.coxph$x

    X = cbind(1, Cova)
    X.invXX = X %*% solve(t(X)%*%X)
    tX = t(X)

    ### calculate empirical CGF for martingale residuals
    idx0 = qcauchy(1:length.out/(length.out+1))
    idx1 = idx0 * max(range) / max(idx0)

    cumul = NULL
    print("Start calculating empirical CGF for martingale residuals...")
    c = 0
    for(i in idx1){
      c = c+1
      t = i
      e_resid = exp(mresid*t)
      M0 = mean(e_resid)
      M1 = mean(mresid*e_resid)
      M2 = mean(mresid^2*e_resid)
      K0 = log(M0)
      K1 = M1/M0
      K2 = (M0*M2-M1^2)/M0^2
      cumul = rbind(cumul, c(t, K0, K1, K2))
      if(c %% 1000 == 0) print(paste0("Complete ",c,"/",length.out,"."))
    }

    K_org_emp = approxfun(cumul[,1], cumul[,2], rule=2)
    K_1_emp = approxfun(cumul[,1], cumul[,3], rule=2)
    K_2_emp = approxfun(cumul[,1], cumul[,4], rule=2)

    var.resid = var(mresid)

    re=list(resid=mresid,
            var.resid=var.resid,
            K_org_emp=K_org_emp,
            K_1_emp=K_1_emp,
            K_2_emp=K_2_emp,
            Call=Call,
            obj.coxph=obj.coxph,
            tX=tX,
            X.invXX=X.invXX,
            gIDs)

    class(re)<-"EmpSPA_Null_Model"
  }
  else if(traits=="binary"){
    Call = match.call()
    ### Fit a logistic model
    obj.logistic = glm(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.logistic, range)

    ### Get the covariate matrix to adjust for genotype
    mu = obj.logistic$fitted.values
    mresid = obj.logistic$y - mu
    Cova = obj.logistic$x

    X = Cova
    X.invXX = X %*% solve(t(X)%*%X)
    tX = t(X)

    ### calculate empirical CGF for martingale residuals
    idx0 = qcauchy(1:length.out/(length.out+1))
    idx1 = idx0 * max(range) / max(idx0)

    cumul = NULL
    print("Start calculating empirical CGF for martingale residuals...")
    c = 0
    for(i in idx1){
      c = c+1
      t = i
      e_resid = exp(mresid*t)
      M0 = mean(e_resid)
      M1 = mean(mresid*e_resid)
      M2 = mean(mresid^2*e_resid)
      K0 = log(M0)
      K1 = M1/M0
      K2 = (M0*M2-M1^2)/M0^2
      cumul = rbind(cumul, c(t, K0, K1, K2))
      if(c %% 1000 == 0) print(paste0("Complete ",c,"/",length.out,"."))
    }

    K_org_emp = approxfun(cumul[,1], cumul[,2], rule=2)
    K_1_emp = approxfun(cumul[,1], cumul[,3], rule=2)
    K_2_emp = approxfun(cumul[,1], cumul[,4], rule=2)

    var.resid = var(mresid)

    re=list(resid=mresid,
            var.resid=var.resid,
            K_org_emp=K_org_emp,
            K_1_emp=K_1_emp,
            K_2_emp=K_2_emp,
            Call=Call,
            obj.logistic=obj.logistic,
            tX=tX,
            X.invXX=X.invXX,
            gIDs)

    class(re)<-"EmpSPA_Null_Model"
  }
  else{
    #resid = resid
    #Cova.mtx = Cova.mtx
   re = EmpSPA_R_Null_Model(Cova.mtx, resid, pIDs = pIDs, gIDs = gIDs, range = range, length.out = length.out)
  }
  return(re)
}

#' Uses residuals from a null generalized linear model to calculate its empirical cumulant generation function (CGF)
#'
#' Calculate the empirical cumulant generation function (CGF) of the raw residuals from a generalized linear regression model
#' @param Cova.mtx a covariate matrix including a column of 1s
#' @param resid a vector of residuals from a genearlized linear regression model in which X is the covariate matrix
#' @param pIDs a character vector of subject IDs. NOTE: its order should be the same as the subjects order in the covariate matrix X.
#' @param gIDs a character vector of subject IDs. NOTE: its order should be the same as the subjects order of the Geno.mtx (i.e. the input of the function EmpSPA()).
#' @param range a two-element numeric vector (default: c(-100,100)) to specify the domain of the empirical CGF.
#' @param length.out a positive integer (default: 9999) for empirical CGF. Larger length.out corresponds to longer calculation time and more accurate estimated empirical CGF.
#' @examples
#'# example 1
#'# Simulation phenotype and genotype
#'N = 10000
#'nSNP = 1000
#'MAF = 0.3
#'X1 = rnorm(N)
#'X2 = rbinom(N,1,0.5)
#'mu = 0.1+0.5*X1+0.5*X2
#'Y = rnorm(N, mu, 1)
#'Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                      Cov1 = X1,
#'                      Cov2 = X2,
#'                      y = Y)
#'obj.linear = lm(Y~Cov1+Cov2, data=Phen.mtx)
#'R = obj.linear$residuals
#'X = cbind(1, X1, X2)
#'C = t(X) %*% R
#'Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)

#'# NOTE: The row and column names of genotype matrix are required.
#'rownames(Geno.mtx) = paste0("IID-",1:N)
#'colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#'Geno.mtx[1:10,1]=NA   # please use NA for missing genotype
#'obj.null = EmpSPA_R_Null_Model(X, R, pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx), range=c(-100,100), length.out = 10000)
#' @export

EmpSPA_R_Null_Model = function(Cova.mtx=NULL,
                               resid=NULL,
                               pIDs,
                               gIDs,
                               range,
                               length.out)
{
  X = Cova.mtx
  R = resid
  Call = match.call()
  p2g = check_input_R(pIDs, gIDs, R, range)

  ### Get the covariate matrix to adjust for genotype
  mresid = R
  X.invXX = X %*% solve(t(X)%*%X)
  tX = t(X)

  ### calculate empirical CGF for residuals
  idx0 = qcauchy(1:length.out/(length.out+1))
  idx1 = idx0 * max(range) / max(idx0)

  cumul = NULL
  print("Start calculating empirical CGF for martingale residuals...")
  c = 0
  for(i in idx1){
    c = c+1
    t = i
    e_resid = exp(mresid*t)
    M0 = mean(e_resid)
    M1 = mean(mresid*e_resid)
    M2 = mean(mresid^2*e_resid)
    K0 = log(M0)
    K1 = M1/M0
    K2 = (M0*M2-M1^2)/M0^2
    cumul = rbind(cumul, c(t, K0, K1, K2))
    if(c %% 1000 == 0) print(paste0("Complete ",c,"/",length.out,"."))
  }

  K_org_emp = approxfun(cumul[,1], cumul[,2], rule=2)
  K_1_emp = approxfun(cumul[,1], cumul[,3], rule=2)
  K_2_emp = approxfun(cumul[,1], cumul[,4], rule=2)

  var.resid = var(R)

  re=list(resid=R,
          var.resid=var.resid,
          K_org_emp=K_org_emp,
          K_1_emp=K_1_emp,
          K_2_emp=K_2_emp,
          Call=Call,
          #obj.coxph=obj.coxph,
          tX=tX,
          X.invXX=X.invXX,
          gIDs)

  class(re)<-"EmpSPA_Null_Model"
  return(re)
}

######################################################################################
#'Empirical SaddlePoint Approximation implementation of a case-control(binary phenotype), surival(time-to-event phenotype) or analysis of other traits for homogeneous case.
#'
#' A fast and accurate method for a genome-wide case-control study, survival analysis or analysis of other traits on a large-scale dataset.
#' @param obj.null an R object returned from function EmpSPA_Null_Model().
#' @param Geno.mtx a numeric genotype matrix with each row as an individual and each column as a genetic variant.
#'                 Column names of genetic variations and row names of subject IDs are required.
#'                 Missng genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param min.maf a numeric value (default: 0.0001) to specify the cutoff of the minimal MAF. Any SNP with MAF < cutoff will be excluded from the analysis.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#'                      "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e. 2p where p is MAF).
#' @param missing.cutoff a numeric value (default: 0.15) to specify the cutoff of the missing rates.
#'                       Any variant with missing rate higher than this cutoff will be excluded from the analysis.
#' @details To run EmpSPA.homo, the following two steps are required:
#' \itemize{
#'   \item Step 1: Use function EmpSPA_Null_Model() to fit a null logistic regression model for binary phenotype or a null Cox model for time-to-event phenotype.
#'   \item Step 2: Use function EmpSPA() to calculate p value for each genetic variant.
#' }
#'
#' EmpSPA calculate p values with both empirical saddlepoint approximation(EmpSPA) and empirical normal distribution approximation(Empnorm).
#' Generally speaking, empirical saddlepoint approximation(EmpSPA) is more accurate than, but a little slower than, the empirical normal approximation(Empnorm).
#' When the score statistic is very close to 0 (i.e. p-values are not small), empirical saddlepoint approximation(EmpSPA) may be numerically unstable .
#' To solve this problem, we replace the wrong p value from EmpSPA by the p value from Empnorm.
#'
#' To calibrate the score statistics for Cox PH model, EmpSPA uses martingale residuals which are calculated via R package survival.
#' All extentions (such as strata, ties, left-censoring) supported by package survival could also be used in EmpSPA.
#' Time-varying covariates are also supported by splitting each subject into several observations.
#' Simulation studies and real data analyses indicate that SPA works well if one subject corresponds to 2~3 observations.
#' While, if there are more than 4 observations for each subject, EmpSPA has not been fully evaluated and the results should be carefully intepreted.
#'
#' Sometimes, the order of subjects between phenotype data and genotype data are different, which could lead to some errors.
#' To avoid that, we ask users to specify the IDs of both phenotype data (pIDs) and genotype data (gIDs) when fitting the null model.
#' Users are responsible to check the consistency between pIDs and formula, and the consistency between gIDs and Geno.mtx.
#'
#' @return an R matrix with the following columns
#' \item{MAF}{Minor allele frequencies}
#' \item{missing.rate}{Missing rates}
#' \item{p.value.empspa}{p value (recommanded) from a empirical saddlepoint approximation.}
#' \item{p.value.empnorm}{p value from a empirical normal distribution approximation.}
#' \item{Stat}{score statistics}
#' \item{Var}{estimated variances of the score statistics}
#' \item{z}{z values corresponding to the score statistics}
#' @examples
#' # example 1  time-to-event phenotype (homogeneous)
#' Simulation phenotype and genotype
#' N = 10000
#' nSNP = 1000
#' MAF = 0.1
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       event=rbinom(N,1,0.5),
#'                       time=runif(N),
#'                       Cov1=rnorm(N),
#'                       Cov2=rbinom(N,1,0.5))
#' Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)
#'
#' # NOTE: The row and column names of genotype matrix are required.
#' rownames(Geno.mtx) = paste0("IID-",1:N)
#' colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#' Geno.mtx[1:10,1]=NA   # please use NA for missing genotype
#'
#' # Attach the survival package so that we can use its function Surv()
#' library(survival)
#' obj.null = EmpSPA_Null_Model("survival" ,Surv(time,event)~Cov1+Cov2, data=Phen.mtx,
#'                              pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
#' survival.res = EmpSPA.homo(obj.null, Geno.mtx)
#'
#' # we recommand using column of 'p.value.empspa' to associate genotype with time-to-event phenotypes
#' head(survival.res)
#'
#' # example 2  binary phenotype (homogeneous)
#' # Simulation phenotype and genotype
#' N = 10000
#' nSNP = 1000
#' MAF = 0.3
#' X1 = rnorm(N)
#' X2 = rbinom(N,1,0.5)
#' mu = exp(0.1+0.5*X1+0.5*X2) /(1+ exp(0.1+0.5*X1+0.5*X2))
#' Y= rbinom(N, 1, mu)
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       Cov1 = X1,
#'                       Cov2 = X2,
#'                       y = Y)
#' Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)
#'
#' # NOTE: The row and column names of genotype matrix are required.
#' rownames(Geno.mtx) = paste0("IID-",1:N)
#' colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#' Geno.mtx[1:10,1]=NA   # please use NA for missing genotype
#'
#' obj.null = EmpSPA_Null_Model("binary",y~Cov1+Cov2,family=binomial(link="logit"), data=Phen.mtx,
#'                              pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
#' binary.res = EmpSPA.homo(obj.null, Geno.mtx)
#'
#' # we recommand using column of 'p.value.empspa' to associate genotype with binary phenotypes
#' head(binary.res)
#'
#' # example 3 continuous phenotype (traits = "others") (homogeneous)
#  # Simulation phenotype and genotype
#'N = 10000
#'nSNP = 1000
#'MAF = 0.3
#'X1 = rnorm(N)
#'X2 = rbinom(N,1,0.5)
#'mu = 0.1+0.5*X1+0.5*X2
#'Y = rnorm(N, mu, 1)
#'Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                      Cov1 = X1,
#'                      Cov2 = X2,
#'                      y = Y)
#'obj.linear = lm(Y~Cov1+Cov2, data=Phen.mtx)
#'R = obj.linear$residuals
#'X = cbind(1, X1, X2)
#'C = t(X) %*% R
#'Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)
#'# NOTE: The row and column names of genotype matrix are required.
#'rownames(Geno.mtx) = paste0("IID-",1:N)
#'colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#'Geno.mtx[1:10,1]=NA   # please use NA for missing genotype
#'obj.null = EmpSPA_Null_Model ("others",resid=R, Cova.mtx=X, data=Phen.mtx, pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx), range=c(-100,100), length.out = 10000)
#'linear.res = EmpSPA.homo(obj.null, Geno.mtx)
#'head(linear.res)
#'@export

EmpSPA.homo = function(obj.null,
                       Geno.mtx,
                       impute.method = "fixed",
                       missing.cutoff = 0.15,
                       min.maf = 0.0001,
                       G.model = "Add")
{
  ## check input
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  impute.method=impute.method,
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf,
                  G.model=G.model)

  # check_input1(obj.null, Geno.mtx, par.list)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 11)
  colnames(output) = c("MAF","missing.rate","p.value.empspa","p.value.norm1","p.value.norm2", "p.value.empspa2", "Stat","Var1","z1","Var2","z2")
  rownames(output) = colnames(Geno.mtx)

  ### Start analysis
  print("Start Analyzing...")
  print(Sys.time())

  # Cycle for genotype matrix
  for(i in 1:n.Geno){

    g = Geno.mtx[,i]
    output.one.SNP = EmpSPA.one.SNP.homo(g,
                                         obj.null,
                                         #Cutoff,
                                         impute.method,
                                         missing.cutoff,
                                         min.maf,
                                         #CovAdj.cutoff,
                                         G.model)
    output[i,] = output.one.SNP
    output[i,3]<-ifelse(abs(output[i,4] - output[i,3]) - 0.5 > 0,
                        output[i,4],output[i,3])
    output[i,6]<-ifelse(abs(output[i,5] - output[i,6]) - 0.5 > 0,
                        output[i,5],output[i,6])
  }

  print("Analysis Complete.")
  print(Sys.time())
  return(output)
}

#'Empirical SaddlePoint Approximation implementation of a case-control(binary phenotype), surival(time-to-event phenotype) or analysis of other traits for heterogeneous case.
#'
#' A fast and accurate method for a genome-wide case-control study, survival analysis or analysis of other traits on a large-scale dataset.
#' @param obj.null an R object returned from function EmpSPA_Null_Model().
#' @param Geno.mtx a numeric genotype matrix with each row as an individual and each column as a genetic variant.
#'                 Column names of genetic variations and row names of subject IDs are required.
#'                 Missng genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param min.maf a numeric value (default: 0.0001) to specify the cutoff of the minimal MAF. Any SNP with MAF < cutoff will be excluded from the analysis.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#'                      "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e. 2p where p is MAF).
#' @param missing.cutoff a numeric value (default: 0.15) to specify the cutoff of the missing rates.
#'                       Any variant with missing rate higher than this cutoff will be excluded from the analysis.
#' @details To run EmpSPA.nonhomo, the following two steps are required:
#' \itemize{
#'   \item Step 1: Use function EmpSPA_Null_Model() to fit a null logistic regression model for binary phenotype or a null Cox model for time-to-event phenotype.
#'   \item Step 2: Use function EmpSPA() to calculate p value for each genetic variant.
#' }
#'
#' EmpSPA calculate p values with both empirical saddlepoint approximation(EmpSPA) and empirical normal distribution approximation(Empnorm).
#' Generally speaking, empirical saddlepoint approximation(EmpSPA) is more accurate than, but a little slower than, the empirical normal approximation(Empnorm).
#' When the score statistic is very close to 0 (i.e. p-values are not small), empirical saddlepoint approximation(EmpSPA) may be numerically unstable .
#' To solve this problem, we replace the wrong p value from EmpSPA by the p value from Empnorm.
#'
#' To calibrate the score statistics for Cox PH model, EmpSPA uses martingale residuals which are calculated via R package survival.
#' All extentions (such as strata, ties, left-censoring) supported by package survival could also be used in EmpSPA.
#' Time-varying covariates are also supported by splitting each subject into several observations.
#' Simulation studies and real data analyses indicate that SPA works well if one subject corresponds to 2~3 observations.
#' While, if there are more than 4 observations for each subject, EmpSPA has not been fully evaluated and the results should be carefully intepreted.
#'
#' Sometimes, the order of subjects between phenotype data and genotype data are different, which could lead to some errors.
#' To avoid that, we ask users to specify the IDs of both phenotype data (pIDs) and genotype data (gIDs) when fitting the null model.
#' Users are responsible to check the consistency between pIDs and formula, and the consistency between gIDs and Geno.mtx.
#'
#' @return an R matrix with the following columns
#' \item{MAF}{Minor allele frequencies}
#' \item{missing.rate}{Missing rates}
#' \item{p.value.empspa}{p value (recommanded) from a empirical saddlepoint approximation.}
#' \item{p.value.empnorm}{p value from a empirical normal distribution approximation.}
#' \item{Stat}{score statistics}
#' \item{Var}{estimated variances of the score statistics}
#' \item{z}{z values corresponding to the score statistics}
#' @examples
#'
#' # example 1  time-to-event phenotype (heterogeneous)
#' # Simulation phenotype and genotype
#' N = 10000
#' nSNP = 1000
#' MAF = 0.1
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       event=rbinom(N,1,0.5),
#'                       time=runif(N),
#'                       Cov1=rnorm(N),
#'                       Cov2=rbinom(N,1,0.5))
#' Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)
#'
#' # NOTE: The row and column names of genotype matrix are required.
#' rownames(Geno.mtx) = paste0("IID-",1:N)
#' colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#' Geno.mtx[1:10,1]=NA   # please use NA for missing genotype
#'
#' # Attach the survival package so that we can use its function Surv()
#' library(survival)
#' obj.null = EmpSPA_Null_Model("survival" ,Surv(time,event)~Cov1+Cov2, data=Phen.mtx,
#'                              pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
#' survival.res = EmpSPA.hete(obj.null, Geno.mtx)
#'
#' # we recommand using column of 'p.value.empspa' to associate genotype with time-to-event phenotypes
#' head(survival.res)
#'
#' # example 2  binary phenotype (heterogeneous)
#' # Simulation phenotype and genotype
#' N = 10000
#' nSNP = 1000
#' MAF = 0.1
#' X1 = rnorm(N)
#' X2 = rbinom(N,1,0.5)
#' mu = exp(0.1+0.5*X1+0.5*X2) /(1+ exp(0.1+0.5*X1+0.5*X2))
#' Y= rbinom(N, 1, mu)
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       Cov1 = X1,
#'                       Cov2 = X2,
#'                       y = Y)
#' Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)
#'
#' # NOTE: The row and column names of genotype matrix are required.
#' rownames(Geno.mtx) = paste0("IID-",1:N)
#' colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#' Geno.mtx[1:10,1]=NA   # please use NA for missing genotype
#'
#' obj.null = EmpSPA_Null_Model("binary",y~Cov1+Cov2,family=binomial(link="logit"), data=Phen.mtx,
#'                              pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
#' binary.res = EmpSPA.hete(obj.null, Geno.mtx)
#'
#' # we recommand using column of 'p.value.empspa' to associate genotype with binary phenotypes
#' head(binary.res)
#'
#'# example 3 continuous phenotype (traits = "others") (heterogeneous)
#'# Simulation phenotype and genotype
#'N = 10000
#'nSNP = 1000
#'MAF = 0.3
#'X1 = rnorm(N)
#'X2 = rbinom(N,1,0.5)
#'mu = 0.1+0.5*X1+0.5*X2
#'Y = rnorm(N, mu, 1)
#'Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                      Cov1 = X1,
#'                      Cov2 = X2,
#'                      y = Y)
#'obj.linear = lm(Y~Cov1+Cov2, data=Phen.mtx)
#'R = obj.linear$residuals
#'X = cbind(1, X1, X2)
#'C = t(X) %*% R
#'Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)

#'# NOTE: The row and column names of genotype matrix are required.
#'rownames(Geno.mtx) = paste0("IID-",1:N)
#'colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#'Geno.mtx[1:10,1]=NA   # please use NA for missing genotype
#'obj.null = EmpSPA_Null_Model ("others",resid=R, Cova.mtx=X, data=Phen.mtx, pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx), range=c(-100,100), length.out = 10000)
#'linear.res = EmpSPA.hete(obj.null, Geno.mtx)
#'head(linear.res)
#' @export

EmpSPA.hete = function(obj.null,
                          Geno.mtx,
                          impute.method = "fixed",
                          missing.cutoff = 0.15,
                          min.maf = 0.0001,
                          G.model = "Add")
{
  ## check input
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  impute.method=impute.method,
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf,
                  G.model=G.model)

  # check_input1(obj.null, Geno.mtx, par.list)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 11)
  colnames(output) = c("MAF","missing.rate","p.value.empspa","p.value.norm1","p.value.norm2", "p.value.empspa2", "Stat","Var1","z1","Var2","z2")
  rownames(output) = colnames(Geno.mtx)

  ### Start analysis
  print("Start Analyzing...")
  print(Sys.time())

  # Cycle for genotype matrix
  for(i in 1:n.Geno){

    g = Geno.mtx[,i]
    output.one.SNP = EmpSPA.one.SNP.hete(g,
                                         obj.null,
                                         #Cutoff,
                                         impute.method,
                                         missing.cutoff,
                                         min.maf,
                                         #CovAdj.cutoff,
                                         G.model)
    output[i,] = output.one.SNP
    output[i,3]<-ifelse(abs(output[i,4] - output[i,3]) - 0.5 > 0,
                        output[i,4],output[i,3])
    output[i,6]<-ifelse(abs(output[i,5] - output[i,6]) - 0.5 > 0,
                        output[i,5],output[i,6])
  }

  print("Analysis Complete.")
  print(Sys.time())
  return(output)
}

#' Empirical SaddlePoint Approximation implementation of logistic regression case-control analysis, Cox regression surival analysis or analysis of other traits (One-SNP-version for homogeneous case)
#'
#' One-SNP-version EmpSPA function. This function is to facilitate users that prefer reading and analyzing genotype line-by-line.
#' @param g a numeric genotype vector. Missing genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param others the same as function EmpSPA.homo. NOTE that we do not check subject order in this one-snp-version !!!
#' @return the same as function EmpSPA.
#' @export
EmpSPA.one.SNP.homo = function(g,
                               obj.null,
                               impute.method = "fixed",
                               missing.cutoff = 0.15,
                               min.maf = 0.0001,
                               G.model = "Add")
{
  ## calculate MAF and update genotype vector
  MAF = mean(g, na.rm=T)/2
  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N

  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }

  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }

  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)

  if(MAF < min.maf)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA))

  #if(!is.null(obj.null$p2g))
  # g = g[obj.null$p2g]

  ## Score statistic
  S = sum(g * obj.null$resid)

  ## estimated variance without adjusting for covariates
  G1 = g - 2*MAF   # centered genotype (such that mean=0)
  S.var1 = obj.null$var.resid * sum(G1^2)
  z1 = S/sqrt(S.var1)
  g.var.est = 2 * MAF * (1 - MAF)
  S.var2 = sum((obj.null$resid)^2 * g.var.est)
  z2 = S/sqrt(S.var2)

  N1set = which(g!=0)  # position of non-zero genotypes
  N0 = N-length(N1set)

  # G1norm = G1/sqrt(S.var1)  # normalized genotype (such that sd=1)

  G1N1 = G1[N1set]
  G1N0 = -2*MAF   # all subjects with g=0 share the same normlized genotype, this is to reduce computation time

  re1 = GetProb_SPA(obj.null, G1N1, G1N0, N1set, N0, abs(S), lower.tail = FALSE) # EmpSPA p value and a
  re2 = GetProb_SPA(obj.null, G1N1, G1N0, N1set, N0, -abs(S), lower.tail = TRUE) # EmpSPA p value and a
  a1 = re1[2]
  a2 = re2[2]

  pval1 = re1[1] # EmpSPA
  pval2 = re2[1] # EmpSPA
  pval3 = pnorm(abs(z1), lower.tail = FALSE) # Normal 1
  pval4 = pnorm(-abs(z1), lower.tail = TRUE) # Normal 1
  pval5 = pnorm(abs(z2), lower.tail = FALSE) # Normal 2
  pval6 = pnorm(-abs(z2), lower.tail = TRUE) # Normal 2
  pval7 = pnorm(a1 * sqrt(S.var1/S.var2), lower.tail = FALSE) # EmpSPA 2
  pval8 = pnorm(a2 * sqrt(S.var1/S.var2), lower.tail = TRUE) # EmpSPA 2
  pval.empspa = pval1 + pval2
  pval.norm1 = pval3 + pval4
  pval.norm2 = pval5 + pval6
  pval.empspa2 = pval7 + pval8

  pval = c(pval.empspa,  pval.norm1, pval.norm2, pval.empspa2) # 4 elements: element 1 is from empspa, element 2 is from Normal 1, element 3 is from Normal 2, element 4 is from

  return(c(MAF, missing.rate, pval, S, S.var1, z1, S.var2, z2))
}

#' Empirical SaddlePoint Approximation implementation of logistic regression case-control analysis, Cox regression surival analysis or analysis of other traits (One-SNP-version for heterogeneous case)
#'
#' One-SNP-version EmpSPA function. This function is to facilitate users that prefer reading and analyzing genotype line-by-line.
#' @param g a numeric genotype vector. Missing genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param others the same as function EmpSPA.homo. NOTE that we do not check subject order in this one-snp-version !!!
#' @return the same as function EmpSPA.nonhomo.
#' @export

EmpSPA.one.SNP.hete = function(g,
                               obj.null,
                               impute.method = "fixed",
                               missing.cutoff = 0.15,
                               min.maf = 0.0001,
                               G.model = "Add")
{
  ## calculate MAF and update genotype vector
  MAF = mean(g, na.rm=T)/2
  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N

  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }

  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }

  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)

  if(MAF < min.maf)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA))

  #if(!is.null(obj.null$p2g))
  # g = g[obj.null$p2g]

  ## Score statistic
  S = sum(g * obj.null$resid)

  ## estimated variance without adjusting for covariates
  N1set = 1:N
  N0 = 0
  G1 = g - obj.null$X.invXX %*% (obj.null$tX[,N1set,drop=F] %*% g[N1set])   # centered genotype (such that mean=0)
  S.var1 = obj.null$var.resid * sum(G1^2)
  z1 = S/sqrt(S.var1)
  MAF.est = 0.5 * obj.null$X.invXX %*% (obj.null$tX[,N1set,drop=F] %*% g[N1set])
  g.var.est = 2 * MAF.est * (1 - MAF.est)
  S.var2 = sum((obj.null$resid)^2 * g.var.est)
  z2 = S/sqrt(S.var2)

  # N1set = which(g!=0)  # position of non-zero genotypes
  # N0 = N-length(N1set)

  # G1norm = G1/sqrt(S.var1)  # normalized genotype (such that sd=1)
  # N1set = 1:N
  # N0 = 0
  G1N1 = G1
  G1N0 = 0   # since N0=0, this value actually does not matter
  # G1N1 = G1[N1set]
  # G1N0 = -2*MAF   # all subjects with g=0 share the same normlized genotype, this is to reduce computation time

  re1 = GetProb_SPA(obj.null, G1N1, G1N0, N1set, N0, abs(S), lower.tail = FALSE) # EmpSPA p value and a
  re2 = GetProb_SPA(obj.null, G1N1, G1N0, N1set, N0, -abs(S), lower.tail = TRUE) # EmpSPA p value and a
  a1 = re1[2]
  a2 = re2[2]

  pval1 = re1[1] # EmpSPA
  pval2 = re2[1] # EmpSPA
  pval3 = pnorm(abs(z1), lower.tail = FALSE) # Normal 1
  pval4 = pnorm(-abs(z1), lower.tail = TRUE) # Normal 1
  pval5 = pnorm(abs(z2), lower.tail = FALSE) # Normal 2
  pval6 = pnorm(-abs(z2), lower.tail = TRUE) # Normal 2
  pval7 = pnorm(a1 * sqrt(S.var1/S.var2), lower.tail = FALSE) # EmpSPA 2
  pval8 = pnorm(a2 * sqrt(S.var1/S.var2), lower.tail = TRUE) # EmpSPA 2
  pval.empspa = pval1 + pval2
  pval.norm1 = pval3 + pval4
  pval.norm2 = pval5 + pval6
  pval.empspa2 = pval7 + pval8

  pval = c(pval.empspa,  pval.norm1, pval.norm2, pval.empspa2) # 4 elements: element 1 is from empspa, element 2 is from Normal 1, element 3 is from Normal 2, element 4 is from

  return(c(MAF, missing.rate, pval, S, S.var1, z1, S.var2, z2))
}


GetProb_SPA = function(obj.null, G2NB, G2NA, NBset, N0, q2, lower.tail){

  out = uniroot(K1_adj, c(-20,20), extendInt = "upX",
                G2NB=G2NB, G2NA=G2NA, NBset=NBset,
                N0=N0, q2=q2, obj.null=obj.null)
  zeta = out$root

  k1 = K_org(zeta,  G2NB=G2NB, G2NA=G2NA, NBset=NBset, N0=N0, obj.null=obj.null)
  k2 = K2(zeta,  G2NB=G2NB, G2NA=G2NA, NBset=NBset, N0=N0, obj.null=obj.null)

  temp1 = zeta * q2 - k1

  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}
  a = w + 1/w * log(v/w)

  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  #pval.norm = pnorm(q2, lower.tail = lower.tail)

  re = c(pval, a)
  return(re)
}

K_org = function(t, G2NB, G2NA, NBset, N0, obj.null){

  n.t = length(t)
  out = rep(0,n.t)
  for(i in 1:n.t){
    t1 = t[i]
    t2NA = t1*G2NA
    t2NB = t1*G2NB
    out[i] = N0*obj.null$K_org_emp(t2NA) + sum(obj.null$K_org_emp(t2NB))
  }
  return(out)
}

K1_adj = function(t, G2NB, G2NA, NBset, N0, q2, obj.null)
{
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    t2NA = t1*G2NA
    t2NB = t1*G2NB
    out[i] = N0*G2NA*obj.null$K_1_emp(t2NA) + sum(G2NB*obj.null$K_1_emp(t2NB)) - q2
  }
  return(out)
}

K2 = function(t, G2NB, G2NA, NBset, N0, obj.null)
{
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    t2NA = t1*G2NA
    t2NB = t1*G2NB
    out[i] = N0*G2NA^2*obj.null$K_2_emp(t2NA) + sum(G2NB^2*obj.null$K_2_emp(t2NB))
  }
  return(out)
}

check_input = function(pIDs, gIDs, obj, range)
{
  if(is.null(pIDs) & is.null(gIDs)) stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  if(any(sort(unique(pIDs))!=sort(unique(gIDs)))) stop("unique(pIDs) should be the same as unique(gIDs).")
  if(anyDuplicated(gIDs)!=0) stop("Argument 'gIDs' should not have a duplicated element.")
  if(range[2]!=-1*range[1]) stop("range[2] should be -1*range[1]")
  mresid = obj$residuals
  if(length(mresid)!=length(pIDs)) stop("Argument 'pIDs' should be of the same length as input data.")

  if(all(pIDs == gIDs)) p2g = NULL
  else p2g = match(pIDs, gIDs)

  return(p2g)
}

check_input_R = function(pIDs, gIDs, R, range)
{
  if(is.null(pIDs) & is.null(gIDs)) stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  if(any(sort(unique(pIDs))!=sort(unique(gIDs)))) stop("unique(pIDs) should be the same as unique(gIDs).")
  if(anyDuplicated(gIDs)!=0) stop("Argument 'gIDs' should not have a duplicated element.")
  if(range[2]!=-1*range[1]) stop("range[2] should be -1*range[1]")
  mresid = R
  if(length(mresid)!=length(pIDs)) stop("Argument 'pIDs' should be of the same length as input data.")

  if(all(pIDs == gIDs)) p2g = NULL
  else p2g = match(pIDs, gIDs)

  return(p2g)
}


# G = Geno.mtx.PCA
# p.est0 = 0.5 * colMeans(G)
# p.est0.mtx = matrix(0, N, nSNP)
# for (i in 1:N) {
#   p.est0.mtx[i,] = p.est0
# }
# Z = (G - 2 * p.est0.mtx)/sqrt(2 * p.est0.mtx * (1 - p.est0.mtx))
# GRM = Z %*% t(Z) / nSNP
#
# pca.pr<-princomp(covmat = GRM)
# summary(pca.pr,loadings = TRUE)
# PC1 = pca.pr$loadings[,1]
# PC2 = pca.pr$loadings[,2]
# PC3 = pca.pr$loadings[,3]
# PC4 = pca.pr$loadings[,4]
# library(ggplot2)
# qplot(PC1, PC2)
# top10_PC = pca.pr$loadings[,1:10]

#'Calculate principle components from an empirical genomic relationship matrix (GRM) when all samples are independent.
#'
#' Use an empirical genomic relationship matrix (GRM) to calculate principle components when all samples are independent.
#' @param Geno.mtx a numeric genotype matrix with each row as an individual and each column as a genetic variant.
#' @return top 10 PCs
#' @examples
#'
#' # example 1  time-to-event phenotype (non-homogeneous)
#' # Simulation phenotype and genotype
#'N = 4000
#'N1 = N/2
#'N2 = N/2
#'nSNP = 10000
#'p.PCA = runif(nSNP, min = 0.2, max = 0.5)
#'p1.PCA = runif(nSNP, min = 0.2, max = 0.5)
#'p2.PCA = runif(nSNP, min = 0.2, max = 0.5)
#'P.PCA = cbind(p1.PCA, p2.PCA)
#'### i.i.d subject & i.i.d SNPS matrix
#'P1.mtx.PCA = matrix(0, N1, nSNP)
#'for (i in 1:N1) {
#'  P1.mtx.PCA[i,]= p1.PCA
#'}
#'P2.mtx.PCA = matrix(0, N2, nSNP)
#'for (i in 1:N2) {
#'  P2.mtx.PCA[i,]= p2.PCA
#'}
#'P.mtx.PCA = rbind(P1.mtx.PCA, P2.mtx.PCA) ### true MAF
#'
#'Geno.mtx1.PCA = matrix(0, N1, nSNP)
#'for (i in 1:N1) {
#'  Geno.mtx1.PCA[i,]= rbinom(nSNP, 2, p1.PCA)
#'}
#'Geno.mtx2.PCA = matrix(0, N2, nSNP)
#'for (i in 1:N2) {
#'  Geno.mtx2.PCA[i,]= rbinom(nSNP, 2, p2.PCA)
#'}
#'Geno.mtx.PCA = rbind(Geno.mtx1.PCA, Geno.mtx2.PCA)
#'
#'GRM_PCA(Geno.mtx.PCA)
#'
#' @export

GRM_PCA = function(Geno.mtx)
{
  p.est0 = (1/2) * colMeans(Geno.mtx)
  p.est0.mtx = matrix(0, N, nSNP)
  for (i in 1:N) {
    p.est0.mtx[i,] = p.est0
  }
  Z = (Geno.mtx - 2 * p.est0.mtx)/sqrt(2 * p.est0.mtx * (1 - p.est0.mtx))
  GRM = Z %*% t(Z) / nSNP
  PC = princomp(covmat = GRM)
  top10_PC = PC$loadings[,1:10]
  return(top10_PC)
}
