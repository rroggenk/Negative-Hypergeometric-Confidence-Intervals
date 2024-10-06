###############################################################
#-------------------------------------------------------------#
# Conditional Minimal Cardinality Procedure (CMC)             # 
# Schilling, Holladay, and Doi (submitted 2020)               #
# binom n unknown AND NB r unknown                            #
#-------------------------------------------------------------#
###############################################################

###################################
#---------------------------------#
# Source R Scripts                #
#---------------------------------#
###################################

# Preliminary functions                                       
source('PreliminaryFunctions.R', encoding = 'UTF-8')

###################################
#---------------------------------#
# CG, OC, BK; binom n unknown     #
#---------------------------------#
###################################


CMC.binomn <-function(p,x.max,conf.level=.90){
  LL=NA
  UL=NA
  
  ## Determine m(a), for all a<=x.max+1. Additionally, 
  # determine rbreak, the smallest value of r
  # such that P(a<=X<=b|r)>=conf.level
  m.vec=NA
  nbreak.vec=NA
  
  # Only acceptance set when n=0 is AS(n=0)={0} so 0-0 must
  # have probability 1 at n=0
  m.vec[1]=0
  nbreak.vec[1]=0
  
  a=1;b=1
  stop=0
  while(a<=x.max+10 | stop<10){
    # For each a increment b by 1 until AC(a,b) first exceeds the 
    # confidence level; then set m(a)=b.  For each a we can 
    # start the loop at b=m(a-1)+1 b/c m(a)>=m(a-1)+1
    while(binomn(x=a:b,n=max.binomn(a,b,p=p),p=p)<conf.level){
      b=b+1
    }
    m.vec[a+1]=b
    
    n.break=max.binomn(a,b,p=p)
    while(binomn(x=a:b, n.break-1, p)>=conf.level & n.break-1>=0){
      n.break=n.break-1
    }
    nbreak.vec[a+1]=n.break
    
    a=a+1
    
    b=max(b,a)  #b is currently at m(a) but when we have card 1 m(a) can be less a+1
    #print(c(a,b))
    #so in that case we need to start b at a
    #below is old
    #b=b-2 #I think we can start the search for m(a+1) at m(a) b/c m(a+1)>=m(a) (
    #BUT I HAVE NOT PROVEN THIS YET, so will start at m(a)-2 for now).
    
    # I added the stopping feature below in case RB(x.max+1) is skipped. Normally the upper endpoint for
    # x.max is the value of n when we first start using RB(x.max+1) (i.e. x.max not long used); however,
    # it might be possible that we never use RB(x.max+1) if RB(x.max+k) rises above confidence level at
    # same n or earlier n for any k. If same n then UL for x.max wouldn't change anyways, but if earlier 
    # the UL would be smaller.  Therefore we need to find n.break for subsequent rainbows.  We hope 
    # that continuing until 10 values of n.break are > than the nbreak for x.max+1 will be sufficient
    # Not too happy about arbitrarily choosing 10 but should be sufficient in most cases
    # may try to improve this part in future!
    if(a>x.max+1){
      if(nbreak.vec[a-1]>nbreak.vec[x.max+1]){
        stop=stop+1
      }
    }
    
  }
  
  
  ## Determine confidence interval endpoints for 0<=x<=x.max.  
  #store {a}_r and {b}_r sequences and cov_r at each r
  a.vec=NA
  b.vec=NA
  cov=NA
  
  n=0
  a=0
  b=m.vec[1]
  a.vec[n+1]=0; a.vec
  b.vec[n+1]=b; b.vec
  cov[n+1]=binomn(x=a:b, n, p)
  LL[(0:b)+1]=0
  
  ## Determine lower endpoints for 0<=x<=x.max.
  # Run until l(x.max) determined
  while(is.na(UL[x.max+1])){
    
    
    ## (1) increment r while the current curve P(a<=X<=b|r) remains
    # in play (i.e. P(a<=X<=b|r) remains above conf.level) and 
    # current rainbow still in play (haven't got to rbreak of subseqent
    # rainbow)
    
    # nextRB is the location where we begin to use subsequent rainbow
    # sometimes rainbows are skipped if the core of any subsequent rainbow 
    # rises above conf.level sooner or at the same location as current RB 
    nextRB=min(nbreak.vec[(a+2):length(nbreak.vec)])
    nexta=max(which(nbreak.vec==nextRB))-1  
    
    
    while(binomn(x=a:b, n+1, p)>=conf.level & n+1<nextRB){
      n=n+1;n
      a.vec[n+1]=a
      b.vec[n+1]=b
      cov[n+1]=binomn(x=a:b, n, p)
    }
    
    
    ## (2) If the current curve is no longer in play (and thus we have 
    # exited the loop above) and it is not yet time to transition 
    # to the next rainbow then we transition between curves from 
    # the same rainbow, AC(a,b) to AC(a,b+j), j>=1 when A(a,b) drops below 
    # the confidence.  This corresponds the lower endpoint for b+1,...,b+j
    if(n+1<nextRB){
      b.temp=b
      while(binomn(x=a:(b+1), n+1, p)<conf.level){
        b=b+1
      }
      
      
      b=b+1;b
      n=n+1;n
      LL[(b.temp+2):(b+1)]=n
      
      a.vec[n+1]=a
      b.vec[n+1]=b
      cov[n+1]=binomn(x=a:b, n, p)
    }
    
    ## (3) If its time to transition to "next" rainbows rainbow 
    # RB(a+j).  This transition occurs when the core of RB(a+j) first rises 
    # above the conf.level and this determines the lower endpoints for 
    # x=b+1,...,m(a+j)and the upper end point for a+1,...a+j
    # '+j' instead of '+1' b/c rainbows can be skipped if  a future 
    # rainbow (usually the following one) rises above the conf.level
    # sooner or at the same location as current RB 
    if(n+1==nextRB){
      
      n=n+1
      LL[(b+2):(m.vec[nexta+1]+1)]=nextRB
      UL[(a+1):nexta]=nextRB-1
      
      a=nexta
      b=m.vec[a+1]
      
      a.vec[n+1]=a
      b.vec[n+1]=b
      cov[n+1]=binomn(x=a:b, n, p)
    };
  }
  
  #Truncate LL/UL vectors to only include intervals for 0<=x<=x.max
  LL=LL[1:(x.max+1)]
  UL=UL[1:(x.max+1)]
  
  acc.matrix=data.frame(n=0:n,a=a.vec,b=b.vec,cov=cov);acc.matrix
  CMC=data.frame(x=0:x.max, LL, UL);CMC
  
  return(list(acc=acc.matrix, CMC=CMC))
  
}

###############################################################
#-------------------------------------------------------------#
# Examples                                                    #
#-------------------------------------------------------------#
###############################################################

#CI=CMC.binomn(p=0.90,x.max=200,conf.level=.90)
#CI$acc
#CI$CMC

###################################
#---------------------------------#
# CG, OC, BK; NB r unknown        #
# X=# of failures                 #
#---------------------------------#
###################################

CMC.negbinomr <- function(p, x.max, conf.level){

LL=NA
UL=NA

## Determine m(a), for all a<=x.max+1. Additionally, 
# determine rbreak, the smallest value of r
# such that P(a<=X<=b|r)>=conf.level
m.vec=NA
rbreak.vec=NA

a=0; b=0

stop=0
while(a<=x.max | stop<10 ){
  # For each a increment b by 1 until AC(a,b) first exceeds the 
  # confidence level; then set m(a)=b.  For each a we can 
  # start the loop at b=m(a-1)+1 b/c m(a)>=m(a-1)+1
  while(negbinomr(x=a:b,r=max.negbinomr(a,b,p=p),p=p)<conf.level){b=b+1}
  m.vec[a+1]=b

  r.break=max.negbinomr(a,b,p=p)
  #r=0 is give prob 1 in negbinomr regardless of x so no error here when when a=0
  while(negbinomr(x=a:b, r.break-1, p)>=conf.level & r.break-1>=1){
    r.break=r.break-1
  }
  rbreak.vec[a+1]=r.break
  
  a=a+1
  b=max(b,a)  #b is currently at m(a) but when we have card 1 m(a) can be less a+1
  #so in that case we need to start b at a
  #below is old
  #b=b-2 #I think we can start the search for m(a+1) at m(a) b/c m(a+1)>=m(a) (
  #BUT I HAVE NOT PROVEN THIS YET, so will start at m(a)-2 for now).
  
  # I added the stopping feature below in case RB(x.max+1) is skipped. Normally the upper endpoint for
  # x.max is the value of r when we first start using RB(x.max+1) (i.e. x.max not long used); however,
  # it might be possible that we never use RB(x.max+1) if RB(x.max+k) rises above confidence level at
  # same r or earlier r for any k. If same r then UL for x.max wouldn't change anyways, but if earlier 
  # the UL would be smaller.  Therefore we need to find r.break for subsequent rainbows.  We hope 
  # that continuing until 10 values of r.break are > than the rbreak for x.max+1 will be sufficient
  # Not too happy about arbitrarily choosing 10 but should be sufficient in most cases
  # may try to improve this part in future!
  if(a>x.max+1){
    if(rbreak.vec[a]>rbreak.vec[x.max+2]){
      stop=stop+1
    }
  }
}
rbreak.vec
m.vec



## Determine confidence interval endpoints for 0<=x<=x.max.  
#store {a}_r and {b}_r sequences and cov_r at each r
a.vec=NA
b.vec=NA
cov=NA

r=1
a=max(which(rbreak.vec==1))-1
b=m.vec[max(which(rbreak.vec==1))]
a.vec[r]=a; a.vec
b.vec[r]=b; b.vec
cov[r]=negbinomr(x=a:b, r, p)
LL[(a:b)+1]=1


## Determine lower endpoints for 0<=x<=x.max.
# Run until l(x.max) determined
while(is.na(UL[x.max+1])){
  

    ## (1) increment r while the current curve P(a<=X<=b|r) remains
    # in play (i.e. P(a<=X<=b|r) remains above conf.level) and 
    # current rainbow still in play (haven't got to rbreak of subseqent
    # rainbow)
    
    # nextRB is the location where we begin to use subsequent rainbow
    # sometimes rainbows are skipped if the core of any subsequent rainbow 
    # rises above conf.level sooner or at the same location as current RB 
    nextRB=min(rbreak.vec[(a+2):length(rbreak.vec)])
    nexta=max(which(rbreak.vec==nextRB))-1  
      
  
    while(negbinomr(x=a:b, r+1, p)>=conf.level & r+1<nextRB){
      r=r+1;r
      a.vec[r]=a
      b.vec[r]=b
      cov[r]=negbinomr(x=a:b, r, p)
    }
   
    
    ## (2) If the current curve is no longer in play (and thus we have 
    # exited the loop above) and it is not yet time to transition 
    # to the next rainbow then we transition between curves from 
    # the same rainbow, AC(a,b) to AC(a,b+j), j>=1 when A(a,b) drops below 
    # the confidence.  This corresponds the lower endpoint for b+1,...,b+j
    if(r+1<nextRB){
      b.temp=b
      while(negbinomr(x=a:(b+1), r+1, p)<conf.level){
        b=b+1
      }
      
     
      b=b+1;b
      r=r+1;r
      LL[(b.temp+2):(b+1)]=r
      
      a.vec[r]=a
      b.vec[r]=b
      cov[r]=negbinomr(x=a:b, r, p)
    }
    
    ## (3) If its time to transition to "next" rainbows rainbow 
    # RB(a+j).  This transition occurs when the core of RB(a+j) first rises 
    # above the conf.level and this determines the lower endpoints for 
    # x=b+1,...,m(a+j)and the upper end point for a+1,...a+j
    # '+j' instead of '+1' b/c rainbows can be skipped if  a future 
    # rainbow (usually the following one) rises above the conf.level
    # sooner or at the same location as current RB 
    if(r+1==nextRB){
     
      r=r+1
      LL[(b+2):(m.vec[nexta+1]+1)]=nextRB
      UL[(a+1):nexta]=nextRB-1
      
      a=nexta
      b=m.vec[a+1]
      
      a.vec[r]=a
      b.vec[r]=b
      cov[r]=negbinomr(x=a:b, r, p)
      #print(c(a,b,cov[r]))
    }
}

#Truncate LL/UL vectors to only include intervals for 0<=x<=x.max
LL=LL[1:(x.max+1)]
UL=UL[1:(x.max+1)]

acc.matrix=data.frame(r=1:r,a=a.vec,b=b.vec,cov=cov);acc.matrix
CMC=data.frame(x=0:x.max, LL, UL);CMC

return(list(acc=acc.matrix, CMC=CMC))

}


###############################################################
#-------------------------------------------------------------#
# Examples                                                    #
#-------------------------------------------------------------#
###############################################################

# CI=CMC.negbinomr(p=0.01,x.max=10, conf.level=.90);CI
# CI$CMC
# CI$acc

###################################
#---------------------------------#
# CG, OC, BK; NB r unknown        #
# N=Y=# of trials                 #
#---------------------------------#
###################################

CMC.negbinomrN <- function(p, y.max, conf.level){
  
  LL=NA
  UL=NA
  
  ## Determine m(a), for all a<=x.max+1. Additionally, 
  # determine rbreak, the smallest value of r
  # such that P(a<=X<=b|r)>=conf.level
  m.vec=NA
  rbreak.vec=NA
  
  a=1; b=1
  stop=0
  while(a<=y.max+10 | stop<10){
    # For each a increment b by 1 until AC(a,b) first exceeds the 
    # confidence level; then set m(a)=b.  For each a we can 
    # start the loop at b=m(a-1)+1 b/c m(a)>=m(a-1)+1
    while(negbinomrN(a:b,r=max.negbinomrN(a,b,p=p),p=p)<conf.level){b=b+1}
    m.vec[a]=b
    
    r.break=max.negbinomrN(a,b,p=p)
    #note r=0 is give prob 1 in negbinomrN regardless of y 
    while(negbinomrN(a:b, r.break-1, p)>=conf.level & r.break-1>=1){
      r.break=r.break-1
    }
    rbreak.vec[a]=r.break
    
    a=a+1
    b=max(b,a)  #b is currently at m(a) but when we have card 1 m(a) can be less a+1
    #so in that case we need to start b at a
    
    # I added the stopping feature below in case RB(x.max+1) is skipped. Normally the upper endpoint for
    # x.max is the value of r when we first start using RB(x.max+1) (i.e. x.max not long used); however,
    # it might be possible that we never use RB(x.max+1) if RB(x.max+k) rises above confidence level at
    # same r or earlier r for any k. If same r then UL for x.max wouldn't change anyways, but if earlier 
    # the UL would be smaller.  Therefore we need to find r.break for subsequent rainbows.  We hope 
    # that continuing until 10 values of r.break are > than the rbreak for x.max+1 will be sufficient
    # Not too happy about arbitrarily choosing 10 but should be sufficient in most cases
    # may try to improve this part in future!
    if(a>y.max+1){
      if(rbreak.vec[a-1]>rbreak.vec[y.max+1]){
        stop=stop+1
      }
    }
  }
  
  
  ## Determine confidence interval endpoints for 1<=y<=y.max.  
  #store {a}_r and {b}_r sequences and cov_r at each r
  a.vec=NA
  b.vec=NA
  cov=NA
  
  r=1
  a=max(which(rbreak.vec==1))
  b=m.vec[max(which(rbreak.vec==1))]
  a.vec[r]=a; a.vec
  b.vec[r]=b; b.vec
  cov[r]=negbinomrN(a:b, r, p)
  LL[(a:b)]=1
  
  
  ## Determine lower endpoints for 0<=x<=x.max.
  # Run until l(x.max) determined
  while(is.na(UL[y.max])){
    
    
    ## (1) increment r while the current curve P(a<=X<=b|r) remains
    # in play (i.e. P(a<=X<=b|r) remains above conf.level) and 
    # current rainbow still in play (haven't got to rbreak of subseqent
    # rainbow)
    
    # nextRB is the location where we begin to use subsequent rainbow
    # sometimes rainbows are skipped if the core of any subsequent rainbow 
    # rises above conf.level sooner or at the same location as current RB 
    nextRB=min(rbreak.vec[(a+1):length(rbreak.vec)])
    nexta=max(which(rbreak.vec==nextRB))  
    
    
    while(negbinomrN(a:b, r+1, p)>=conf.level & r+1<nextRB){
      r=r+1;r
      a.vec[r]=a
      b.vec[r]=b
      cov[r]=negbinomrN(a:b, r, p)
    }
    
    
    ## (2) If the current curve is no longer in play (and thus we have 
    # exited the loop above) and it is not yet time to transition 
    # to the next rainbow then we transition between curves from 
    # the same rainbow, AC(a,b) to AC(a,b+j), j>=1 when A(a,b) drops below 
    # the confidence.  This corresponds the lower endpoint for b+1,...,b+j
    if(r+1<nextRB){
      b.temp=b
      while(negbinomrN(a:(b+1), r+1, p)<conf.level){
        b=b+1
      }
      
      
      b=b+1;b
      r=r+1;r
      LL[(b.temp+1):b]=r
      
      a.vec[r]=a
      b.vec[r]=b
      cov[r]=negbinomrN(a:b, r, p)
    }
    
    ## (3) If its time to transition to "next" rainbows rainbow 
    # RB(a+j).  This transition occurs when the core of RB(a+j) first rises 
    # above the conf.level and this determines the lower endpoints for 
    # x=b+1,...,m(a+j)and the upper end point for a+1,...a+j
    # '+j' instead of '+1' b/c rainbows can be skipped if  a future 
    # rainbow (usually the following one) rises above the conf.level
    # sooner or at the same location as current RB 
    if(r+1==nextRB){
      
      r=r+1
      LL[(b+1):(m.vec[nexta+1])]=nextRB
      UL[a:(nexta-1)]=nextRB-1
      
      a=nexta
      b=m.vec[a]   
      
      a.vec[r]=a
      b.vec[r]=b
      cov[r]=negbinomrN(a:b, r, p)
    }
  }
  
  #Truncate LL/UL vectors to only include intervals for 1<=y<=y.max
  LL=LL[1:(y.max)]
  UL=UL[1:(y.max)]
  
  acc.matrix=data.frame(r=1:r,a=a.vec,b=b.vec,cov=cov);acc.matrix
  CMC=data.frame(N=1:y.max, LL, UL);CMC
  
  return(list(acc=acc.matrix, CMC=CMC))
  
}


###############################################################
#-------------------------------------------------------------#
# Examples                                                    #
#-------------------------------------------------------------#
###############################################################

#CI=CMC.negbinomrN(p=0.80,y.max=80, conf.level=.95);CI$acc
# CI$CMC
# CI$acc

CI1=CMC.negbinomrN(p=0.02,y.max=40, conf.level=.90);CI1$acc
CI2=MC.negbinomrN(p=0.02,y.max=40, conf.level=.90);CI2$acc.CG














