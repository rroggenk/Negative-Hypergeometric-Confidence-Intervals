###############################################################
#-------------------------------------------------------------#
# Preliminary functions                                       #
#-------------------------------------------------------------#
###############################################################

###############################################################
#-------------------------------------------------------------#
# Returns probability sums: P(a<=X<=b)=sum_{k=a}^b P(X=k)     # 
# for vector x=a:b for HG, binom, and NB                      #
#-------------------------------------------------------------#
###############################################################

###################################
#---------------------------------#
#  Binomial(n,p); n unknown       #  
#---------------------------------#
###################################

binomn <- function(x, n, p){ 
  prob=NA 
  i=1
  while(i<=length(n)){ prob[i]=sum(dbinom(x,size=n[i],prob=p)); i=i+1 }
  return(prob) 
}

###################################
#---------------------------------#
#  NB(r,p); r unknown, X=# fails  #  
#---------------------------------#
###################################
# X~NB(r,p)represents the number of failures which occur 
# in a sequence of Bernoulli trials before a target number 
# of r successes is reached.

negbinomr <- function(x, r, p){ 
  prob=NA 
  i=1
  while(i<=length(r)){ prob[i]=sum(dnbinom(x,size=r[i],prob=p)); i=i+1 }
  return(prob) 
}

###################################
#---------------------------------#
# NB(r,p); r unknown, N=# trials  #  
#---------------------------------#
###################################
# N=Y=r+X represents the number of trials 
# in a sequence of Bernoulli trials before rth success.

negbinomrN <- function(y, r, p){ 
  prob=NA 
  i=1
  while(i<=length(r)){ prob[i]=sum(dnbinom(y-r[i],size=r[i],prob=p)); i=i+1 }
  return(prob) 
}


###############################################################
#-------------------------------------------------------------#
# Calculates the acceptance set of SCALAR theta for a         #
# given confidence procedure                                  # 
#-------------------------------------------------------------#
###############################################################

###################################
#---------------------------------#
#  Binomial(n,p); n unknown       #  
#---------------------------------#
###################################

accept.binomn <- function(n,p,LL,UL){
  left = floor(n*p)
  right = ceiling(n*p)
  while( LL[left+1] <=n & UL[left+1]>=n ){ left = left-1; if(left<0){break} }
  while( LL[right+1] <=n & UL[right+1]>=n ){ right = right+1; if(right>n){break} }
  return((left+1):(right-1))
}

###################################
#---------------------------------#
#  NB(r,p); r unknown, X=# fails  #  
#---------------------------------#
###################################

accept.negbinomr <- function(r,p,LL,UL){
  left = floor(r*(1-p)/p);left
  right = ceiling(r*(1-p)/p);right
  while( LL[left+1] <=r & UL[left+1]>=r ){ left = left-1; if(left<0){break} }
  while( LL[right+1] <=r & UL[right+1]>=r ){ right = right+1}#; if(right>r){break} }
  return((left+1):(right-1))
}

###################################
#---------------------------------#
# NB(r,p); r unknown, N=# trials  #  
#---------------------------------#
###################################
# N=Y=r+X represents the number of trials 

accept.negbinomrN <- function(r,p,LL,UL){
  left = floor(r/p);left
  right = ceiling(r/p);right
  while( LL[left] <=r & UL[left]>=r ){ left = left-1; if(left<r){break} }
  while( LL[right] <=r & UL[right]>=r ){ right = right+1}
  return((left+1):(right-1))
}


###############################################################
#-------------------------------------------------------------#
# Calculates the coverage of a VECTOR of the unknown          #
# parameter for a given confidence procedure                  # 
#-------------------------------------------------------------#
###############################################################

###################################
#---------------------------------#
#  Binomial(n,p); n unknown       #  
#---------------------------------#
###################################

cov.binomn <- function(n,p,LL,UL){
  cov=n
  card=n
  for(i in 1:length(n)){
    acc=accept.binomn(n=n[i],p=p,LL=LL,UL=UL)
    cov[i]=binomn(acc,n=n[i],p=p)
    card[i]=length(acc)
  }
  return(data.frame(cov=cov,card=card))
}

###################################
#---------------------------------#
#  NB(r,p); r unknown, X=#of fails#  
#---------------------------------#
###################################

cov.negbinomr <- function(r,p,LL,UL){
  cov=r
  card=r
  a=r
  b=r
  for(i in 1:length(r)){
    acc=accept.negbinomr(r=r[i],p=p,LL=LL,UL=UL)
    cov[i]=negbinomr(acc,r=r[i],p=p)
    card[i]=length(acc)
    a[i]=acc[1]
    b[i]=acc[length(acc)]
  }
  return(data.frame(r=r,a=a,b=b,cov=cov,card=card))
}

###################################
#---------------------------------#
# NB(r,p); r unknown, N=# trials  #  
#---------------------------------#
###################################
# N=Y=r+X represents the number of trials 

cov.negbinomrN <- function(r,p,LL,UL){
  cov=r
  card=r
  a=r
  b=r
  for(i in 1:length(r)){
    acc=accept.negbinomrN(r=r[i],p=p,LL=LL,UL=UL)
    cov[i]=negbinomrN(acc,r=r[i],p=p)
    card[i]=length(acc)
    a[i]=acc[1]
    b[i]=acc[length(acc)]
  }
  return(data.frame(r=r,a=a,b=b,cov=cov,card=card))
}


###############################################################
#-------------------------------------------------------------#
# Returns the value of r_max(a,b)                            #
# i.e. the location of the maximum of P(a<=X<=b|r)           #
#-------------------------------------------------------------#
###############################################################



###################################
#---------------------------------#
#  Binomial(n,p); n unknown       #  
#---------------------------------#
###################################

max.binomn <- function(a,b,p=p){
  n.max=b
  if(a>0){
    n.max=b
    while(binomn(x=a:b,n=n.max,p=p)<=binomn(x=a:b,n=n.max+1,p=p)){
        n.max=n.max+1
    }
  }
  return(n.max)
}

###################################
#---------------------------------#
#  NB(r,p); r unknown, X=#of fails#  
#---------------------------------#
###################################


#NB(k,mu)
max.negbinomr <- function(a,b,p=p){
  r.max=1
  if(a>0){
    while(negbinomr(x=a:b,r=r.max,p=p)<=negbinomr(x=a:b,r=r.max+1,p=p) & negbinomr(x=a:b,r=r.max,p=p)<1){
      r.max=r.max+1
    }
  }
  return(r.max)
}


###################################
#---------------------------------#
# NB(r,p); r unknown, N=# trials  #  
#---------------------------------#
###################################
# N=Y=r+X represents the number of trials 

max.negbinomrN <- function(a,b,p=p){
  r.max=1
  if(a>1){
    while(negbinomrN(a:b,r=r.max,p=p)<=negbinomrN(a:b,r=r.max+1,p=p)){
      r.max=r.max+1
    }
  }
  return(r.max)
}





