###############################################################
#-------------------------------------------------------------#
# Minimal Cardinality Procedures (CG, OC, BK)                 #   
# binom n unknown AND NB r unknown                            #
#-------------------------------------------------------------#
###############################################################
#CG:  Crow and Gardner (1959)
#OC: Schilling and Holladay (2017)
#BK: Byrne and Kabaila (2005)

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
#p=0.90;conf.level=.90; x.max=200;acc.list[1:5]


MC.binomn <-function(p,x.max,conf.level=.90){
  #x.max is the largest value of x of interest

  #To store candidate acceptance sets of minimum cardinality at each r
  acc.list=list(NA)
  
  #Only acceptance set when n=0 is AS(n=0)={0}
  acc.list[[1]]=data.frame(n=0,a=0,b=0,card=1,prob=1)
  
  #When n=1 AS(n=1) is either {0}, {0,1} or {1} depending on 
  # p<1-conf.level, 1-conf.level<=p<conf.level, or p>=conf.level
  if(p>=conf.level){
    a.min=0;b.min=0
    cardn.start=0
    cardn.end=0
  }else if(p<1-conf.level){
    a.min=0;b.min=0
    cardn.start=0
    cardn.end=0
  }else{
    a.min=0;b.min=1
    cardn.start=1
    cardn.end=1
  }

  a.max=a.min
  b.max=b.min
  
  #Run until upper endpoint for x.max can be determined for all min. card. procedures
  while(a.min<=x.max){
    
    #find largest a  such that AC(a-b) matching the current cardinality exceeds the confidence level
    while(binomn(x=(a.max+1):(b.max+1), n=max.binomn(a.max+1,b.max+1,p=p), p=p)>=conf.level){
      a.max=a.max+1
      b.max=b.max+1
    }
    
    #determine where the curve found above falls below the confidence level
    # it is at this point where cardinality needs to be increased possibily 
    # by more than one step.  increment cardr.end until curve found above falls
    # below level.
    cardn.end=max.binomn(a.max,b.max,p=p)
    while(binomn(x=a.max:b.max, n=cardn.end+1, p=p)>=conf.level){
      cardn.end=cardn.end+1
    }
    
    # For each n such that the current cardinality is in play (i.e. is still the min. card.) 
    # create matrix with all candidate acceptance sets of minimal cardinality 
    # Each row of the matrix records n, values a,b of the acceptance curves P(a<=X<=b), 
    # cardinality=b-a+1, and coverage P(a<=X<=b|n)
    for(n in cardn.start:cardn.end){
      a=a.min
      b=b.min
      acc.matrix=data.frame(n=NA,a=NA,b=NA,card=NA,prob=NA)
      j=1
      
      # starting from the smallest values of a,b  such that AC(a-b) matching the current cardinality 
      # exceeds the confidence level, increment each a,b by 1 (a=a+1, b=b+1) until the acceptance 
      # curves no longer exceed conf.level.  Store these candidate acceptance sets in acc.list
      while(binomn(x=a:b, n=n, p=p)>=conf.level){
        acc.matrix[j,]=c(n,a,b,b-a+1,binomn(a:b,n=n,p=p))
        a=a+1
        b=b+1
        j=j+1
      }
      acc.list[[n+1]]=acc.matrix
      
      # if we still have iterations left in the 'for' loop (i.e. not yet time to change cardinality), 
      # determine the smallest a  such that AC(a-b) matching the current cardinality exceeds the confidence level
      if(n+1<=cardn.end){
        while(binomn(x=a.min:b.min, n=n+1, p=p)<conf.level){
          a.min=a.min+1
          b.min=b.min+1
        }
      }
    }
    
    #determine a.min for next cardinality in use (may be more than one step above current card.)
    
    #the start location for the next cardinality is value of r after where the curves of the current
    # cardinality dropped below conf.level
    cardn.start=cardn.end+1
    
    #see if any curves of the next cardinality are above the confidence level at the new cardr.start
    # if not increase cardinality until curves above the confidence level our found
    # some cardinalities might be skipped
    stop1=0
    a=a.min
    b=b.min+1
    while(stop1==0){
      while(binomn(x=a:b, n=max.binomn(a,b,p=p), p=p)>=conf.level){
        if(binomn(x=a:b, n=cardn.start, p=p)>=conf.level){
          stop1=1
          break
        }
        a=a+1
        b=b+1
      }
      if(stop1==0){
        a=a.min
        b.min=b.min+1
        b=b.min
      }
    }
    #set the new a.min, b.min, a.max, b.max for next cardinality and repeat all of the above
    a.min=a
    b.min=b
    a.max=a
    b.max=b
  }; acc.list[1:5]
  
  
  # Now that we have list of candidate acceptance set of minimal cardinality at each n
  # we can determine the CI's for CG, OC, and BK
  
  #To store acceptance sets used at each r by CG, OC, and BK
  #These will be selected from the candidate acceptance sets determined above.
  n.vec=0:n
  acc.matrix.CG=data.frame(n=rep(NA,n+1),a=rep(NA,n+1),b=rep(NA,n+1),cov=rep(NA,n+1))
  acc.matrix.OC=data.frame(n=rep(NA,n+1),a=rep(NA,n+1),b=rep(NA,n+1),cov=rep(NA,n+1))
  acc.matrix.BK=data.frame(n=rep(NA,n+1),a=rep(NA,n+1),b=rep(NA,n+1),cov=rep(NA,n+1))
  
  #Find acceptance curve for CG, OC, BK
  for(n in n.vec){
    acc.matrix=acc.list[[n+1]]; acc.matrix
    
    #CG
    #a value of last acceptance curve used
    # to ensure we keep sequences {a} and {b} nondecreasing
    if(n>=1){a.prev=acc.matrix.CG$a[n]
    }else{a.prev=0}
    a.prev
    
    #which min card curves a:b have maximal value of a
    ind=which(acc.matrix$a==max(acc.matrix$a[acc.matrix$a>=a.prev]));ind
    acc.matrix.CG$n[n+1]=n
    acc.matrix.CG$a[n+1]=acc.matrix$a[ind]
    acc.matrix.CG$b[n+1]=acc.matrix$b[ind]
    acc.matrix.CG$cov[n+1]=acc.matrix$prob[ind]
    
    #OC
    # a value of last acceptance curve used
    # to ensure we keep sequences {a} and {b} nondecreasing
    if(n>=1){a.prev=acc.matrix.OC$a[n]
    }else{a.prev=0}
    a.prev
    
    #which min card curve give highest coverage
    ind=which(acc.matrix$prob==max(acc.matrix$prob[acc.matrix$a>=a.prev]));ind
    
    #whenlength(ind)>1 his means there is a tie for largest prob.  
    #In this case choose the acceptance set a-b with largest a
    ind=ind[length(ind)]
    
    acc.matrix.OC$n[n+1]=n
    acc.matrix.OC$a[n+1]=acc.matrix$a[ind]
    acc.matrix.OC$b[n+1]=acc.matrix$b[ind]
    acc.matrix.OC$cov[n+1]=acc.matrix$prob[ind]
    
    #BK
    #a value of last acceptance curve used
    # to ensure we keep sequences {a} and {b} nondecreasing
    if(n>=1){a.prev=acc.matrix.BK$a[n]
    }else{a.prev=0}
    a.prev
    
    #which min card curves a:b have minimal value of a
    ind=which(acc.matrix$a==min(acc.matrix$a[acc.matrix$a>=a.prev]));ind
    acc.matrix.BK$n[n+1]=n
    acc.matrix.BK$a[n+1]=acc.matrix$a[ind]
    acc.matrix.BK$b[n+1]=acc.matrix$b[ind]
    acc.matrix.BK$cov[n+1]=acc.matrix$prob[ind]
  }
  
  x.vec=0:x.max
  CI.CG=list(NA)
  CI.OC=list(NA)
  CI.BK=list(NA)
  for(x in x.vec){
    CI.CG[[x+1]]=acc.matrix.CG$n[x>=acc.matrix.CG$a & x<=acc.matrix.CG$b]
    CI.OC[[x+1]]=acc.matrix.OC$n[x>=acc.matrix.OC$a & x<=acc.matrix.OC$b]
    CI.BK[[x+1]]=acc.matrix.BK$n[x>=acc.matrix.BK$a & x<=acc.matrix.BK$b]
  }
  
  
  CI.CG.LL=rep(NA,length(x.vec))
  CI.CG.UL=rep(NA,length(x.vec))
  CI.OC.LL=rep(NA,length(x.vec))
  CI.OC.UL=rep(NA,length(x.vec))
  CI.BK.LL=rep(NA,length(x.vec))
  CI.BK.UL=rep(NA,length(x.vec))
  for(x in x.vec){
    CI.CG.LL[x+1]=(CI.CG[[x+1]])[1]
    CI.CG.UL[x+1]=(CI.CG[[x+1]])[length(CI.CG[[x+1]])]
    
    CI.OC.LL[x+1]=(CI.OC[[x+1]])[1]
    CI.OC.UL[x+1]=(CI.OC[[x+1]])[length(CI.OC[[x+1]])]
    
    CI.BK.LL[x+1]=(CI.BK[[x+1]])[1]
    CI.BK.UL[x+1]=(CI.BK[[x+1]])[length(CI.BK[[x+1]])]
  }
  
  return(list(#CG.points=CI.CG,OC.points=CI.OC,BK.points=CI.BK,
    CG=data.frame(method='CG',conf.level=conf.level, p=p,x=x.vec,LL=CI.CG.LL,UL=CI.CG.UL),
    OC=data.frame(method='OC',conf.level=conf.level, p=p,x=x.vec,LL=CI.OC.LL,UL=CI.OC.UL),
    BK=data.frame(method='BK',conf.level=conf.level, p=p,x=x.vec,LL=CI.BK.LL,UL=CI.BK.UL),
    acc.CG=acc.matrix.CG,acc.OC=acc.matrix.OC,acc.BK=acc.matrix.BK,
    acc.all=acc.list))
  
}
############
#Example

# CI=MC.binomn(p=0.90,x.max=100, conf.level=.90) 
# CI$OC
# CI$BK
# CI$acc.CG
# CI$acc.OC
# CI$acc.BK
# CI$acc.all

###################################
#---------------------------------#
# CG, OC, BK; NB r unknown        #
# X=# of failures                 #
#---------------------------------#
###################################

MC.negbinomr <- function(p, x.max, conf.level=conf.level){
  #x.max is the largest value of x of interest
  
  #To store candidate acceptance sets of minimum cardinality at each r
  acc.list=list(NA)
  
  #r=1 is just the geometric distribution, which has decreasing probabilities. 
  # So it follows immediately that the min. card. acceptance sets are all of 
  # the form {0,...,b} for all 0<p<1.  So we need to find smallest b such that
  # P(a<=X<=b)>=1-alpha
  a.min=0
  b.min=0
  while(negbinomr(x=a.min:b.min, r=1, p=p)<conf.level){b.min=b.min+1}
  a.max=0
  b.max=b.min
  
  cardr.start=1
  cardr.end=1
  
  #Run until upper endpoint for x.max can be determined for all min. card. procedures
  while(a.min<=x.max){
  
    #find largest a  such that AC(a-b) matching the current cardinality exceeds the confidence level
    while(negbinomr(x=(a.max+1):(b.max+1), r=max.negbinomr(a.max+1,b.max+1,p=p), p=p)>=conf.level){
      a.max=a.max+1
      b.max=b.max+1
    }
    
    #determine where the curve found above falls below the confidence level
    # it is at this point where cardinality needs to be increased possibily 
    # by more than one step.  increment cardr.end until curve found above falls
    # below level.
    cardr.end=max.negbinomr(a.max,b.max,p=p)
    while(negbinomr(x=a.max:b.max, r=cardr.end+1, p=p)>=conf.level){
      cardr.end=cardr.end+1
    }

    # For each r such that the current cardinality is in play (i.e. is still the min. card.) 
    # create matrix with all candidate acceptance sets of minimal cardinality 
    # Each row of the matrix records r, values a,b of the acceptance curves P(a<=X<=b), 
    # cardinality=b-a+1, and coverage P(a<=X<=b|r)
    for(r in cardr.start:cardr.end){
      a=a.min
      b=b.min
      acc.matrix=data.frame(r=NA,a=NA,b=NA,card=NA,prob=NA)
      j=1
      
      # starting from the smallest values of a,b  such that AC(a-b) matching the current cardinality 
      # exceeds the confidence level, increment each a,b by 1 (a=a+1, b=b+1) until the acceptance 
      # curves no longer exceed conf.level.  Store these candidate acceptance sets in acc.list
      while(negbinomr(x=a:b, r=r, p=p)>=conf.level){
        acc.matrix[j,]=c(r,a,b,b-a+1,negbinomr(a:b,r=r,p=p))
        a=a+1
        b=b+1
        j=j+1
      }
      acc.list[[r]]=acc.matrix
      
      # if we still have iterations left in the 'for' loop (i.e. not yet time to change cardinality), 
      # determine the smallest a  such that AC(a-b) matching the current cardinality exceeds the confidence level
      if(r+1<=cardr.end){
        while(negbinomr(x=a.min:b.min, r=r+1, p=p)<conf.level){
          a.min=a.min+1
          b.min=b.min+1
        }
      }
    }
      
      #determine a.min for next cardinality in use (may be more than one step above current card.)
    
      #the start location for the next cardinality is value of r after where the curves of the current
      # cardinality dropped below conf.level
      cardr.start=cardr.end+1
      
      #see if any curves of the next cardinality are above the confidence level at the new cardr.start
      # if not increase cardinality until curves above the confidence level our found
      # some cardinalities might be skipped
      stop1=0
      a=a.min
      b=b.min+1
      while(stop1==0){
        while(negbinomr(x=a:b, r=max.negbinomr(a,b,p=p), p=p)>=conf.level){
          if(negbinomr(x=a:b, r=cardr.start, p=p)>=conf.level){
            stop1=1
            break
          }
          a=a+1
          b=b+1
        }
        if(stop1==0){
          a=a.min
          b.min=b.min+1
          b=b.min
        }
      }
      #set the new a.min, b.min, a.max, b.max for next cardinality and repeat all of the above
      a.min=a
      b.min=b
      a.max=a
      b.max=b
  }
  
  
  #To store acceptance sets used at each r by CG, OC, and BK
  #These will be selected from the candidate acceptance sets determined above.
  r.vec=1:r
  acc.matrix.CG=data.frame(r=rep(NA,r),a=rep(NA,r),b=rep(NA,r),cov=rep(NA,r))
  acc.matrix.OC=data.frame(r=rep(NA,r),a=rep(NA,r),b=rep(NA,r),cov=rep(NA,r))
  acc.matrix.BK=data.frame(r=rep(NA,r),a=rep(NA,r),b=rep(NA,r),cov=rep(NA,r))
  
  #Find acceptance curve for CG, OC, BK
  for(r in r.vec){
    acc.matrix=acc.list[[r]];acc.matrix
    
    #CG
    #a value of last acceptance curve used
    # to ensure we keep sequences {a} and {b} nondecreasing
    if(r>=2){a.prev=acc.matrix.CG$a[r-1]
    }else{a.prev=0}
    a.prev
    
    #which min card curves a:b have maximal value of a
    ind=which(acc.matrix$a==max(acc.matrix$a[acc.matrix$a>=a.prev]));ind
    acc.matrix.CG$r[r]=r
    acc.matrix.CG$a[r]=acc.matrix$a[ind]
    acc.matrix.CG$b[r]=acc.matrix$b[ind]
    acc.matrix.CG$cov[r]=acc.matrix$prob[ind]
    
    #OC
    #a value of last acceptance curve used
    # to ensure we keep sequences {a} and {b} nondecreasing
    if(r>=2){a.prev=acc.matrix.OC$a[r-1]
    }else{a.prev=0}
    a.prev
    
    #which min card curve give highest coverage
    ind=which(acc.matrix$prob==max(acc.matrix$prob[acc.matrix$a>=a.prev]));ind
    
    #whenlength(ind)>1 his means there is a tie for largest prob.  
    #In this case choose the acceptance set a-b with largest a
    ind=ind[length(ind)]
    
    acc.matrix.OC$r[r]=r
    acc.matrix.OC$a[r]=acc.matrix$a[ind]
    acc.matrix.OC$b[r]=acc.matrix$b[ind]
    acc.matrix.OC$cov[r]=acc.matrix$prob[ind]
    
    #BK
    #a value of last acceptance curve used
    # to ensure we keep sequences {a} and {b} nondecreasing
    if(r>=2){a.prev=acc.matrix.BK$a[r-1]
    }else{a.prev=0}
    a.prev
    
    #which min card curves a:b have minimal value of a
    ind=which(acc.matrix$a==min(acc.matrix$a[acc.matrix$a>=a.prev]));ind
    acc.matrix.BK$r[r]=r
    acc.matrix.BK$a[r]=acc.matrix$a[ind]
    acc.matrix.BK$b[r]=acc.matrix$b[ind]
    acc.matrix.BK$cov[r]=acc.matrix$prob[ind]
  }
  
  x.vec=0:x.max
  CI.CG=list(NA)
  CI.OC=list(NA)
  CI.BK=list(NA)
  for(x in x.vec){
    CI.CG[[x+1]]=acc.matrix.CG$r[x>=acc.matrix.CG$a & x<=acc.matrix.CG$b]
    CI.OC[[x+1]]=acc.matrix.OC$r[x>=acc.matrix.OC$a & x<=acc.matrix.OC$b]
    CI.BK[[x+1]]=acc.matrix.BK$r[x>=acc.matrix.BK$a & x<=acc.matrix.BK$b]
  }
  
  CI.CG.LL=rep(NA,length(x.vec))
  CI.CG.UL=rep(NA,length(x.vec))
  CI.OC.LL=rep(NA,length(x.vec))
  CI.OC.UL=rep(NA,length(x.vec))
  CI.BK.LL=rep(NA,length(x.vec))
  CI.BK.UL=rep(NA,length(x.vec))
  for(x in 0:x.max){
    CI.CG.LL[x+1]=(CI.CG[[x+1]])[1]
    CI.CG.UL[x+1]=(CI.CG[[x+1]])[length(CI.CG[[x+1]])]
    
    CI.OC.LL[x+1]=(CI.OC[[x+1]])[1]
    CI.OC.UL[x+1]=(CI.OC[[x+1]])[length(CI.OC[[x+1]])]
    
    CI.BK.LL[x+1]=(CI.BK[[x+1]])[1]
    CI.BK.UL[x+1]=(CI.BK[[x+1]])[length(CI.BK[[x+1]])]
  }
  
  return(list(#CG.points=CI.CG,OC.points=CI.OC,BK.points=CI.BK,
    CG=data.frame(method='CG',conf.level=conf.level, p=p,x=x.vec,LL=CI.CG.LL,UL=CI.CG.UL),
    OC=data.frame(method='OC',conf.level=conf.level, p=p,x=x.vec,LL=CI.OC.LL,UL=CI.OC.UL),
    BK=data.frame(method='BK',conf.level=conf.level, p=p,x=x.vec,LL=CI.BK.LL,UL=CI.BK.UL),
    acc.CG=acc.matrix.CG,acc.OC=acc.matrix.OC,acc.BK=acc.matrix.BK,
    acc.all=acc.list))
  
}

############
#Example

#CI=MC.negbinomr(p=0.01,x.max=10, conf.level=.90)
# CI$CG
# CI$OC
# CI$BK
# CI$acc.CG
# CI$acc.OC
# CI$acc.BK
# CI$acc.all


###################################
#---------------------------------#
# CG, OC, BK; NB r unknown        #
# N=Y=# of trails                 #
#---------------------------------#
###################################

MC.negbinomrN <- function(p, y.max, conf.level=conf.level){
  #y.max is the largest value of y of interest
  
  #To store candidate acceptance sets of minimum cardinality at each r
  acc.list=list(NA)
  
  #r=1 is just the geometric distribution, which has decreasing probabilities. 
  # So it follows immediately that the min. card. acceptance sets are all of 
  # the form {r,...,b} for all 0<p<1.  So we need to find smallest b such that
  # P(a<=X<=b)>=1-alpha
  a.min=1
  b.min=1
  while(negbinomrN(a.min:b.min, r=1, p=p)<conf.level){b.min=b.min+1}
  a.max=1
  b.max=b.min
  
  cardr.start=1
  cardr.end=1
  
  #Run until upper endpoint for x.max can be determined for all min. card. procedures
  while(a.min<=y.max){
    
    #find largest a  such that AC(a-b) matching the current cardinality exceeds the confidence level
    while(negbinomrN((a.max+1):(b.max+1), r=max.negbinomrN(a.max+1,b.max+1,p=p), p=p)>=conf.level){
      a.max=a.max+1
      b.max=b.max+1
    }
    
    #determine where the curve found above falls below the confidence level
    # it is at this point where cardinality needs to be increased possibily 
    # by more than one step.  increment cardr.end until curve found above falls
    # below level.
    cardr.end=max.negbinomrN(a.max,b.max,p=p)
    while(negbinomrN(a.max:b.max, r=cardr.end+1, p=p)>=conf.level){
      cardr.end=cardr.end+1
    }
    
    # For each r such that the current cardinality is in play (i.e. is still the min. card.) 
    # create matrix with all candidate acceptance sets of minimal cardinality 
    # Each row of the matrix records r, values a,b of the acceptance curves P(a<=X<=b), 
    # cardinality=b-a+1, and coverage P(a<=X<=b|r)
    for(r in cardr.start:cardr.end){
      a=a.min
      b=b.min
      acc.matrix=data.frame(r=NA,a=NA,b=NA,card=NA,prob=NA)
      j=1
      
      # starting from the smallest values of a,b  such that AC(a-b) matching the current cardinality 
      # exceeds the confidence level, increment each a,b by 1 (a=a+1, b=b+1) until the acceptance 
      # curves no longer exceed conf.level.  Store these candidate acceptance sets in acc.list
      while(negbinomrN(a:b, r=r, p=p)>=conf.level){
        acc.matrix[j,]=c(r,a,b,b-a+1,negbinomrN(a:b,r=r,p=p))
        a=a+1
        b=b+1
        j=j+1
      }
      acc.list[[r]]=acc.matrix
      
      # if we still have iterations left in the 'for' loop (i.e. not yet time to change cardinality), 
      # determine the smallest a  such that AC(a-b) matching the current cardinality exceeds the confidence level
      if(r+1<=cardr.end){
        while(negbinomrN(a.min:b.min, r=r+1, p=p)<conf.level){
          a.min=a.min+1
          b.min=b.min+1
        }
      }
    }
    
    #determine a.min for next cardinality in use (may be more than one step above current card.)
    
    #the start location for the next cardinality is value of r after where the curves of the current
    # cardinality dropped below conf.level
    cardr.start=cardr.end+1
    
    #see if any curves of the next cardinality are above the confidence level at the new cardr.start
    # if not increase cardinality until curves above the confidence level our found
    # some cardinalities might be skipped
    stop1=0
    a=a.min
    b=b.min+1
    while(stop1==0){
      while(negbinomrN(a:b, r=max.negbinomrN(a,b,p=p), p=p)>=conf.level){
        if(negbinomrN(a:b, r=cardr.start, p=p)>=conf.level){
          stop1=1
          break
        }
        a=a+1
        b=b+1
      }
      if(stop1==0){
        a=a.min
        b.min=b.min+1
        b=b.min
      }
    }
    #set the new a.min, b.min, a.max, b.max for next cardinality and repeat all of the above
    a.min=a
    b.min=b
    a.max=a
    b.max=b
  }
  
  
  #To store acceptance sets used at each r by CG, OC, and BK
  #These will be selected from the candidate acceptance sets determined above.
  r.vec=1:r
  acc.matrix.CG=data.frame(r=rep(NA,r),a=rep(NA,r),b=rep(NA,r),cov=rep(NA,r))
  acc.matrix.OC=data.frame(r=rep(NA,r),a=rep(NA,r),b=rep(NA,r),cov=rep(NA,r))
  acc.matrix.BK=data.frame(r=rep(NA,r),a=rep(NA,r),b=rep(NA,r),cov=rep(NA,r))
  
  #Find acceptance curve for CG, OC, BK
  for(r in r.vec){
    acc.matrix=acc.list[[r]];acc.matrix
    
    #CG
    #a value of last acceptance curve used
    # to ensure we keep sequences {a} and {b} nondecreasing
    if(r>=2){a.prev=acc.matrix.CG$a[r-1]
    }else{a.prev=1}
    a.prev
    
    #which min card curves a:b have maximal value of a
    ind=which(acc.matrix$a==max(acc.matrix$a[acc.matrix$a>=a.prev]));ind
    acc.matrix.CG$r[r]=r
    acc.matrix.CG$a[r]=acc.matrix$a[ind]
    acc.matrix.CG$b[r]=acc.matrix$b[ind]
    acc.matrix.CG$cov[r]=acc.matrix$prob[ind]
    
    #OC
    #a value of last acceptance curve used
    # to ensure we keep sequences {a} and {b} nondecreasing
    if(r>=2){a.prev=acc.matrix.OC$a[r-1]
    }else{a.prev=1}
    a.prev
    
    #which min card curve give highest coverage
    ind=which(acc.matrix$prob==max(acc.matrix$prob[acc.matrix$a>=a.prev]));ind
    
    #whenlength(ind)>1 his means there is a tie for largest prob.  
    #In this case choose the acceptance set a-b with largest a
    ind=ind[length(ind)]
    
    acc.matrix.OC$r[r]=r
    acc.matrix.OC$a[r]=acc.matrix$a[ind]
    acc.matrix.OC$b[r]=acc.matrix$b[ind]
    acc.matrix.OC$cov[r]=acc.matrix$prob[ind]
    
    #BK
    #a value of last acceptance curve used
    # to ensure we keep sequences {a} and {b} nondecreasing
    if(r>=2){a.prev=acc.matrix.BK$a[r-1]
    }else{a.prev=1}
    a.prev
    
    #which min card curves a:b have minimal value of a
    ind=which(acc.matrix$a==min(acc.matrix$a[acc.matrix$a>=a.prev]));ind
    acc.matrix.BK$r[r]=r
    acc.matrix.BK$a[r]=acc.matrix$a[ind]
    acc.matrix.BK$b[r]=acc.matrix$b[ind]
    acc.matrix.BK$cov[r]=acc.matrix$prob[ind]
    
  }
  
  y.vec=1:y.max
  CI.CG=list(NA)
  CI.OC=list(NA)
  CI.BK=list(NA)
  for(y in y.vec){
    CI.CG[[y]]=acc.matrix.CG$r[y>=acc.matrix.CG$a & y<=acc.matrix.CG$b]
    CI.OC[[y]]=acc.matrix.OC$r[y>=acc.matrix.OC$a & y<=acc.matrix.OC$b]
    CI.BK[[y]]=acc.matrix.BK$r[y>=acc.matrix.BK$a & y<=acc.matrix.BK$b]
  }
  
  CI.CG.LL=rep(NA,length(y.vec))
  CI.CG.UL=rep(NA,length(y.vec))
  CI.OC.LL=rep(NA,length(y.vec))
  CI.OC.UL=rep(NA,length(y.vec))
  CI.BK.LL=rep(NA,length(y.vec))
  CI.BK.UL=rep(NA,length(y.vec))
  for(y in 1:y.max){
    CI.CG.LL[y]=(CI.CG[[y]])[1]
    CI.CG.UL[y]=(CI.CG[[y]])[length(CI.CG[[y]])]
    
    CI.OC.LL[y]=(CI.OC[[y]])[1]
    CI.OC.UL[y]=(CI.OC[[y]])[length(CI.OC[[y]])]
    
    CI.BK.LL[y]=(CI.BK[[y]])[1]
    CI.BK.UL[y]=(CI.BK[[y]])[length(CI.BK[[y]])]
  }
  
  return(list(#CG.points=CI.CG,OC.points=CI.OC,BK.points=CI.BK,
    CG=data.frame(method='CG',conf.level=conf.level, p=p,N=y.vec,LL=CI.CG.LL,UL=CI.CG.UL),
    OC=data.frame(method='OC',conf.level=conf.level, p=p,N=y.vec,LL=CI.OC.LL,UL=CI.OC.UL),
    BK=data.frame(method='BK',conf.level=conf.level, p=p,N=y.vec,LL=CI.BK.LL,UL=CI.BK.UL),
    acc.CG=acc.matrix.CG,acc.OC=acc.matrix.OC,acc.BK=acc.matrix.BK,
    acc.all=acc.list))
  
}

############
#Example

#CI=MC.negbinomrN(p=0.80,y.max=25, conf.level=.90)
# CI$CG
# CI$OC
# CI$BK
# CI$acc.CG
# CI$acc.OC
# CI$acc.BK
# CI$acc.all












