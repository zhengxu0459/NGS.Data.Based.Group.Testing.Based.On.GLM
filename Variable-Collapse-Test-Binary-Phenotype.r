# 
# R Script file to implement our variable collapse test for binary phenotype
# Version 1.2 (continuously updated)
# Author: Zheng Xu

variable.collapse.test.binary.phenotype = function(pheno,snp.glf,
  MAF=NULL,cov=NULL,weight=1,
  MAF.upper.bound=0.99, MAF.lower.bound = 0.01){ 

  # the size of snp.glf needs to be a multiple of 3.
  sub.n.marker=(dim(snp.glf)[2])/3 
  use = rep(TRUE,sub.n.marker)
  if (sub.n.marker>1){
  	for(w in 1:(sub.n.marker-1)){
		if (use[w]==TRUE){
			for (ww in (w+1):sub.n.marker){
				if (     all ( snp.glf[,((3*w-2):(3*w))] == 
						 snp.glf[,((3*ww-2):(3*ww))] )   ){
					use[ww]=FALSE
				}

			}
		}	
  	}
  }
  
  sub.n.marker = sum(use)
  snp.glf = snp.glf[,as.numeric(rbind(use,use,use))==1]		
  sub.n.person=length(pheno) 
  
  if (is.null(MAF)){
  	sub.estimated.MAF = rep(0.05,sub.n.marker)
  	sub.minus.log.likelihood = function(p,ThreeColumns){-mean(log((1-p)^2*ThreeColumns[,1]+(2*p*(1-p))*ThreeColumns[,2]+p^2*ThreeColumns[,3]))}
  	for (w in 1:sub.n.marker){
     		sub.estimated.MAF[w] = optimize(sub.minus.log.likelihood,interval=c(MAF.lower.bound,MAF.upper.bound),ThreeColumns=snp.glf[,(3:5)-5+w*3])$minimum
  	}
  }	else{
	sub.estimated.MAF = MAF[use]
  }

  if (is.null(cov)){
    extended.cov = matrix(1,sub.n.person,1)
  } else{
    extended.cov = cbind(1,cov)
  }

  if ( (length(weight)==1) && (weight==1) ){
     weight = matrix(1,1,sub.n.marker)
  } else{
     weight = matrix(weight,1,sub.n.marker)
  }

  sub.n.covariate = dim(extended.cov)[2]
  
  #### for binary phenotype
  
    #### MLE estimator
    if (is.null(cov)){
      temp = summary(glm(pheno~1, family = binomial(link="logit")))  
    } else{
      temp = summary(glm(pheno~cov, family = binomial(link="logit")))  
    }
    glm.result = temp$coefficients[,1] 
    alpha = glm.result   
    eta = extended.cov%*%alpha 
    fitted = exp(eta)/(1+exp(eta))
    deviation = pheno-fitted 
    Sequencing.Given.Genotype.tensor = 
		array(NA,dim = c(sub.n.person,sub.n.marker,3))
    Sequencing.Given.Genotype.tensor[,,1] = 
		as.matrix(snp.glf[,-5+5+c(seq(from=1,to=3*sub.n.marker,by=3))])
    Sequencing.Given.Genotype.tensor[,,2] = 
		as.matrix(snp.glf[,-5+5+c(seq(from=2,to=3*sub.n.marker,by=3))])
    Sequencing.Given.Genotype.tensor[,,3] = 
		as.matrix(snp.glf[,-5+5+c(seq(from=3,to=3*sub.n.marker,by=3))])
    MAF.matrix = matrix(sub.estimated.MAF[], 
		byrow = TRUE, sub.n.person, sub.n.marker)
    Genotype.Given.MAF.tensor = array(NA, dim = c(sub.n.person,sub.n.marker,3) )
    Genotype.Given.MAF.tensor[,,1] = 
		(1 - MAF.matrix)^2
    Genotype.Given.MAF.tensor[,,2] = 
		2*(1 - MAF.matrix)*MAF.matrix
    Genotype.Given.MAF.tensor[,,3] = 
		MAF.matrix^2
    Joint.Sequencing.Genotype.tensor =  
		Sequencing.Given.Genotype.tensor * 	Genotype.Given.MAF.tensor
    
    if (sub.n.marker==1){
      Posterior.Expected.Genotype.matrix = 
		(1*matrix(Joint.Sequencing.Genotype.tensor[,,2])+2*matrix(Joint.Sequencing.Genotype.tensor[,,3]))/(matrix(Joint.Sequencing.Genotype.tensor[,,1])+matrix(Joint.Sequencing.Genotype.tensor[,,2])+matrix(Joint.Sequencing.Genotype.tensor[,,3]))
    } else{
      Posterior.Expected.Genotype.matrix = 
		(1*Joint.Sequencing.Genotype.tensor[,,2]+2*Joint.Sequencing.Genotype.tensor[,,3])/(Joint.Sequencing.Genotype.tensor[,,1]+Joint.Sequencing.Genotype.tensor[,,2]+Joint.Sequencing.Genotype.tensor[,,3])
    }
   
    Posterior.Expected.Genotype.Genotype.tensor.Numerator=
		array(0,dim=c(sub.n.person,sub.n.marker,sub.n.marker))
    Posterior.Expected.Genotype.Genotype.tensor.Denominator=
		array(0,dim=c(sub.n.person,sub.n.marker,sub.n.marker))
    for (k in 1:sub.n.marker){
      for (j in 1:sub.n.marker){
         if (k!=j){
            for (g1 in 0:2){
              for (g2 in 0:2){
                 Posterior.Expected.Genotype.Genotype.tensor.Numerator[,k,j] = 	
					Posterior.Expected.Genotype.Genotype.tensor.Numerator[,k,j]+g1*g2*Joint.Sequencing.Genotype.tensor[,k,(g1+1)]*Joint.Sequencing.Genotype.tensor[,j,(g2+1)]
                 Posterior.Expected.Genotype.Genotype.tensor.Denominator[,k,j] = 
					Posterior.Expected.Genotype.Genotype.tensor.Denominator[,k,j]+Joint.Sequencing.Genotype.tensor[,k,(g1+1)]*Joint.Sequencing.Genotype.tensor[,j,(g2+1)]
              }
           }
        } else{
            for (g1 in 0:2){
              for (g2 in 0:2){
      	         if (g1==g2){ 
                      Posterior.Expected.Genotype.Genotype.tensor.Numerator[,k,j] = 
						Posterior.Expected.Genotype.Genotype.tensor.Numerator[,k,j]+g1*g2*Joint.Sequencing.Genotype.tensor[,k,(g1+1)]
                      Posterior.Expected.Genotype.Genotype.tensor.Denominator[,k,j] = 
						Posterior.Expected.Genotype.Genotype.tensor.Denominator[,k,j]+Joint.Sequencing.Genotype.tensor[,k,(g1+1)]
                     }  
              }
            } 
        }
      }
    }
    Posterior.Expected.Genotype.Genotype.tensor = 
		Posterior.Expected.Genotype.Genotype.tensor.Numerator/Posterior.Expected.Genotype.Genotype.tensor.Denominator
      
    # score function and Hessian information matrix
    Score.Beta.T = 
		matrix(apply( matrix(deviation,sub.n.person,sub.n.marker,byrow=FALSE)*Posterior.Expected.Genotype.matrix,2,sum),sub.n.marker,1) 
    Score.Alpha.T = matrix(0,sub.n.covariate,1)
    Score.Phi = 0
    Hessian.Alpha.T.Alpha = matrix(0,sub.n.covariate,sub.n.covariate)
    for (k in 1:sub.n.covariate)
      for (j in 1:sub.n.covariate)
        Hessian.Alpha.T.Alpha[k,j] =  
			sum(exp(eta) / (1+exp(eta))^2 * extended.cov[,k] * extended.cov[,j]) 
    
    Hessian.Alpha.T.Beta = matrix(0,sub.n.covariate,sub.n.marker)
    for (k in 1:sub.n.covariate)
      for (j in 1:sub.n.marker)
        Hessian.Alpha.T.Beta[k,j]=  
			sum ( exp(eta) / (1+exp(eta))^2 * extended.cov[,k] * Posterior.Expected.Genotype.matrix[,j] ) 
    Hessian.Alpha.T.Phi = matrix(0,sub.n.covariate,1)
    
    Hessian.Beta.T.Beta = matrix(0,sub.n.marker,sub.n.marker)
    for (k in 1:sub.n.marker){
      for (j in 1:sub.n.marker){
        Hessian.Beta.T.Beta[k,j] = 
			- sum( deviation^2 * (Posterior.Expected.Genotype.Genotype.tensor[,k,j]-Posterior.Expected.Genotype.matrix[,k]*Posterior.Expected.Genotype.matrix[,j])- exp(eta) / (1+exp(eta))^2  * Posterior.Expected.Genotype.Genotype.tensor[,k,j]    )  
      }
    }

    Hessian.Beta.T.Phi = matrix(0,sub.n.marker,1) 
    Hessian.Phi.Phi = 0 

    Score.vector = rbind(Score.Alpha.T,Score.Beta.T,Score.Phi)
    Hessian.matrix = 
      rbind( cbind(Hessian.Alpha.T.Alpha,Hessian.Alpha.T.Beta,Hessian.Alpha.T.Phi),
             cbind(t(Hessian.Alpha.T.Beta),Hessian.Beta.T.Beta,Hessian.Beta.T.Phi),
             cbind(t(Hessian.Alpha.T.Phi),t(Hessian.Beta.T.Phi),Hessian.Phi.Phi)) 
    Reduced.Score.vector = rbind(Score.Alpha.T,Score.Beta.T)
    Reduced.Hessian.matrix = 
      rbind( cbind(Hessian.Alpha.T.Alpha,Hessian.Alpha.T.Beta),
             cbind(t(Hessian.Alpha.T.Beta),Hessian.Beta.T.Beta) ) 
    Score = t(Reduced.Score.vector)%*%solve(Reduced.Hessian.matrix)%*%Reduced.Score.vector
   
	# Apply the chain rule for variable collapse test
	Score.Beta0 = weight %*% Score.Beta.T
	Hessian.Alpha.T.Beta0 = Hessian.Alpha.T.Beta %*% t(weight)
	Hessian.Beta0.Beta0 = weight %*% Hessian.Beta.T.Beta %*% t(weight)
	Hessian.Beta0.Phi = weight %*% Hessian.Beta.T.Phi
	Collapsed.Score.vector = rbind(Score.Alpha.T,Score.Beta0,Score.Phi)
	Collapsed.Hessian.matrix = 
	rbind( cbind(Hessian.Alpha.T.Alpha,Hessian.Alpha.T.Beta0,Hessian.Alpha.T.Phi),
		   cbind(t(Hessian.Alpha.T.Beta0),Hessian.Beta0.Beta0,Hessian.Beta0.Phi),
		   cbind(t(Hessian.Alpha.T.Phi),t(Hessian.Beta0.Phi),Hessian.Phi.Phi))
	Collapsed.Reduced.Score.vector = rbind(Score.Alpha.T,Score.Beta0)
	Collapsed.Reduced.Hessian.matrix = 
	rbind( cbind(Hessian.Alpha.T.Alpha,Hessian.Alpha.T.Beta0),
		   cbind(t(Hessian.Alpha.T.Beta0),Hessian.Beta0.Beta0) )
	Collapsed.Score = t(Collapsed.Reduced.Score.vector)%*%solve(Collapsed.Reduced.Hessian.matrix)%*%Collapsed.Reduced.Score.vector 
  
	# Score statistic and p value for variable collapse test for binary phenotype
	p.value=1-pchisq(Score,df=sub.n.marker)
	collapsed.p.value=1-pchisq(Collapsed.Score,df=1)
	
	return(list( df=sub.n.marker, score.statistic=Collapsed.Score, p.value=collapsed.p.value))
} 

