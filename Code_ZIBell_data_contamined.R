#######################################################################################
#
#    Programme de simulation MDPDE ZI-Bell (Methode  d'estimation de Power divergence)
#              Les données Y_c,i  contaminées
#######################################################################################

install.packages("bellreg")
install.packages("LambertW")
install.packages("lamW")


library(bellreg)
library(LambertW)
library(lamW)

rm(list=ls())  
      
data_c=function(n,b,g) {
  
  inter=rep(1,n)
  X1=rnorm(n,0,1)
  X=rbind(inter,X1)
  Z=inter
  w= exp(t(g)%*%Z)/(1+exp(t(g)%*%Z)) # zero-inflation probability
  mu=exp(t(b)%*%X)
  
  Y=c()
  S=rbinom(n,1,w) 
  for (i in 1:n)
  {    
    Y[i]=0*(S[i]==1)+rbell(1,W(mu[i]))*(S[i]!=1) 
  }
  list(Y=Y, X=X, Z=Z)
  
  P=rbinom(n,1,0.02 ) 
  Y_0=rpois(n,12) 
  
  Q=c()
  
  for (i in 1:n)
  {    
    
    Q[i]=P[i]*Y_0[i] 
  }
  list(Q=Q, P=P, Y_0=Y_0)
  
  Y_c=c()
  for (i in 1:n)
    
    if(Y[i]==0)
      Y_c[i]=Y[i]+Q[i]
  else
    Y_c[i]=Y[i]
  list(Y_c=Y_c,X=X, Z=Z)   
}
b=c(-0.5, 1.2); g=-1.1  ## 0.25
#b=c(-0.5, 1.2); g=.2 
Dc500=data_c(500,b,g)  


#############  Calcul de g(y, mu, w)########################

g.density=function(y, mu, w){
  
  Result=(w+(1-w)*dbell(y,W(mu)))*as.integer(y==0)+(1-w)*dbell(y,W(mu))*as.integer(y!=0)
  
  return(Result)
}

################## Calcul de Ell_alpha ##############

Ell_alpha= function(y_d, mu,w,alpha, N=20){
  
  M=matrix(0:N,N+1,1)
  
  G.D=function(y)
    (g.density(y, mu, w))^(1+alpha)
  
  S=apply(M,1,G.D)
  
  sum(S) -(1+1/alpha)*(g.density(y_d, mu, w))^(alpha)
  
}
###################  Calcul de H_alpha ##################

Fonc_H=function(theta, Y_c,X,Z, alpha){
  
  b=theta[1:2]; g=theta[3]
  
  w_i= exp(t(g)%*%Z)/(1+exp(t(g)%*%Z)) # zero-inflation probability
  mu_i=exp(t(b)%*%X)
  
  M.y=matrix(1:length(Y_c),length(Y_c), 1)
  
  ELL=function(i)
    Ell_alpha(Y_c[i], mu_i[i],w_i[i],alpha)
  
  sum(apply(M.y, 1, ELL))/length(Y_c)
  
  
}

################ Fonction d'estimation #########################


Estime_Rep =function(Y_c,X,Z, alpha)                       
{
  Theta_In=  c(.1,.1,.1) # Le vecteur initial pour d?marrer l'algorithme
  
  lik.opt=function(theta)
  { Fonc_H(theta, Y_c,X,Z, alpha) } 
  
  
  opt = optim(par=Theta_In, lik.opt,hessian = TRUE ,  method = "L-BFGS-B")
  
  
  return(opt$par )
}

Estime_Rep(Dc500$Y_c,Dc500$X,Dc500$Z, alpha=0.1) #### renvoie les paramètres estimés


########### Fonction pour enregistrer les données générer ##########

Simul_rep=function(n,nbrep, b, g){      # nbrep=nrow(My)
  
  My=matrix(0, nbrep, n)
  Mx=array(0,dim=c(2,n,nbrep)); Mx[1,,]=1
  
  for (i in 1:nbrep) {
    X.sim = data_c(n, b, g)
    My[i,]=X.sim$Y; Mx[,,i]= X.sim$X
  } 
  list(My=My, Mx=Mx)
}
CData500_ZI25= Simul_rep(500,500,b=c(-0.5, 1.2), g=-1.1    )
save(CData500_ZI25,file = "CData500_ZI25.rdata")


load("CData500_ZI25.rdata")
donnees=CData500_ZI25

####################### Monté Carlo #######################

MC = function(donnees,alpha)
{
  b=c(-0.5, 1.2); g=-1.1
  True_theta= c(b,g)
  My=donnees$My; Mx=donnees$Mx; Z=rep(1, ncol(My))
  
  Mat.Est= matrix(0,nrow(My), length(True_theta))
  
  for ( k in 1:nrow(My))             #nrow(My)=nbrep
  { print(k)
    Mat.Est[k,] = Estime_Rep(My[k,],Mx[,,k],Z, alpha ) 
  }
  
  Mean = apply(Mat.Est,2,mean)
  biais= apply(Mat.Est,2,mean)-True_theta
  Rmse= sqrt( apply((Mat.Est - matrix(True_theta,nrow(My),ncol=length(True_theta), byrow=T))^2,2,mean) )
  
  list( True_theta=True_theta, Mean.est =round(Mean,3) ,biais=round(biais,3),Rmse=round(Rmse,3))
  
}


###############  Resultats  ##########################
b=c(-0.5, 1.2); g=-1.1 #0.25
#b=c(-0.5, 1.2); g=.2 

T1<-Sys.time()
Res= MC( CData500_ZI25, alpha=0.3)  ###### Ici il faut changer alpha= 0.1
T2<-Sys.time() ; Temps = T2-T1; Temps


