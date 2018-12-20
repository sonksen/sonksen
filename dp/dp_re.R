library(msm)

covernt3words <- function(c,phi){
  k <- length(c)
  c.old <- c
  phi.t <- phi[c]
  c <- rep(2,k)
  c[phi.t==min(phi.t)] <-length(unique(c.old)) 
  c[phi.t==max(phi.t)]<-1
  word <- as.integer(paste(c[1],c[2],c[3],sep=''))
  return(word)
}

rphi<-function(){
	rnorm(1)
}

dphi<-function(phi){
	dnorm(phi,log=TRUE)
}
dc<-function(i,ci,c,alpha){
	k<-length(c)
	log.l<-0
	if(i==1){return(0)}
	nic<-sum(c[1:(i-1)]==ci)
	if(nic==0){
		log.l<-log.l+log(alpha/(i-1+alpha))
	}else{
		log.l<-log.l+log(nic/(i-1+alpha))
	}
	return(log.l)
}

rc<-function(i,c,alpha){
	if(i==1){return(1)}
	l<-sort(unique(c[1:(i-1)]))
	k<-length(l)
	nic<-rep(0,k)
	for(j in 1:k){
	      nic[j]<-sum(c[1:j]==l[j])
	}
	ci<-c(1:k,k+1)[rmultinom(1,1,c(nic,alpha)/(i-1+alpha))==1]
	return(ci)
}

sync<-function(c,phi){
  k<-max(c)
  c.new<-c
  for(i in 1:k){
    c.new[ c==((1:k)[phi[1:k]==sort(phi[1:k],decreasing = TRUE)[i]]) ]<- i
  }
  phi[1:k]<-sort(phi[1:k],decreasing=TRUE)
  list(phi=phi,c=c.new)
}

psi <- function(x,n){
  sqrt((n+1)+1/(1+n*x)^2)
}

Cn <- function(n){
  (1-sqrt(n-1))/(sqrt(n)+sqrt(n-1))^3
}





dalpha<-function(alpha,a=1,b=1){
  dgamma(alpha,a,scale=b,log=TRUE)
}



update.alpha<-function(alpha,phi,c,a=1,b=1,n=3){
  k<-length(unique(c))
  eta<-rbeta(1,alpha+1,n)
  temp<-(a+k-1)/(n*(b-log(eta)))
  p<-exp(temp)/(1+exp(temp))
  if(rbinom(1,1,p)==1){
    alpha <- rgamma(1,a+k,  rate=b-log(eta))
  }else{
    alpha <- rgamma(1,a+k-1,rate=b-log(eta))
  }
  return(alpha)
}



update.c<-function(theta,c,phi,y,alpha,re){
       n.y <- length(y)
	k<-length(c)
        mx<-length(y)
	for(i in 1:k){
	      if(sum(c[i]==c)!=1){
		c_star<-c
		c_star[i]<-max(c)+1
		phi_star<-c(phi,rphi())
		a<-min(1, (alpha/(k-1))*exp(pi.theta(theta,n.y,
                       sync(c_star,phi_star))-
		        pi.theta(theta,n.y,sync(c,phi))))
		if(rbinom(1,1,a)==1){ 
			c<-c_star
			phi<-phi_star
		}
	      }else{
	        c_star<-c
		c_star[i]<-sample(c[1:k!=i],1)
		a<-min(1, ((k-1)/(alpha))*exp(pi.theta(theta,n.y,                                                   sync(c_star,phi))-
			        pi.theta(theta,n.y,sync(c,phi))  ))
		if(rbinom(1,1,a)==1){ 
			phi<-phi[1:length(phi)!=c[i]]
			c<-c_star
		}
		
	      }
	      c.list<-sort(unique(c))
	      c.temp<-c
	      for(j in 1:length(c.list)){
	      	    c.temp[c.list[j]==c]<-j
	      }
	      c<-c.temp
	}
	for(i in 1:k){
	      if(sum(c[i]==c)!=1){
		probs<-rep(0,k)
		for(j in 1:k){
		      c_alt<-c
		      c_alt[i]<-c[j]
		      probs[j]<-(1/(k-1))*exp(pi.theta(theta,n.y,sync(c_alt,phi)))
	      	}
		probs<-probs/sum(probs)
		c[i]<-sample(c,1,prob=probs)
	      }
	}
        sy<-sync(c,phi)
        c<-sy$c
        phi <- sy$phi
	list(c=c,phi=phi)
}



update.theta <- function(theta,re,c,y,sd.t,phi){
  syn <- NULL
  syn$c <- c
  syn$phi <- phi
  n <- length(y)
  prop <- theta
  prop[1] <- rnorm(1,theta[1],sd.t[1])
  log.r <- log.lik(y,prop,re)-log.lik(y,theta,re)+
    pi.theta(prop,n,syn)-pi.theta(theta,n,syn)
  if(rbinom(1,1,min(exp(log.r),1))==1) theta <- prop

  prop <- theta
  prop[2] <- rtnorm(1,theta[2],sd.t[2],0,Inf)
  log.r <- log.lik(y,prop,re)-log.lik(y,theta,re)+
    pi.theta(prop,n,syn)-pi.theta(theta,n,syn)
  if(rbinom(1,1,min(exp(log.r),1))==1) theta <- prop

  prop <- theta
  prop[3] <- rtnorm(1,theta[3],sd.t[3],0,Inf)
  log.r <- log.lik(y,prop,re)-log.lik(y,theta,re)+
    pi.theta(prop,n,syn)-pi.theta(theta,n,syn)+
      pi.re(re,prop)-pi.re(re,theta)
  if(rbinom(1,1,min(exp(log.r),1))==1) theta <- prop

  return(theta)

}





update.phi<-function(theta,c,phi,y){
	k<-length(phi)
        n.y <- length(y)
	for(i in 1:k){
	      prop<-phi
	      prop[i]<-rphi()
	log.l<-pi.theta(theta,n.y,sync(c,prop))+dphi(prop[i])-
	       pi.theta(theta,n.y,sync(c, phi))-dphi(phi[i])
	if(rbinom(1,1,min(1,exp(log.l)))==1) phi<-prop	
	      
	}
        sy<-sync(c,phi)
        c<-sy$c
        phi <- sy$phi
	list(c=c,phi=phi)
}












pi.theta <- function(theta,n,syn){
  sigma2 <- theta[2]
  tau2 <- theta[3]
  mu <- theta[1]
  phi <- syn$phi
  c <- syn$c
  k <- length(unique(c))
  c.old <- c
  phi.t <- phi[c]
  c <- rep(2,k)
  c[phi.t==min(phi.t)] <-length(unique(c.old)) 
  c[phi.t==max(phi.t)]<-1
  word <- as.integer(paste(c[1],c[2],c[3],sep=''))
  if(word==111){
    log.l <- -log(sigma2)-(3/2)*log(n*tau2 + sigma2)
  }
  if(word==112){
    log.l <- -2.5*log(sqrt(sigma2)) -log(n*tau2 + sigma2)
  }
  if(word==121){
    log.l <- -2*(Cn(n)/2)*log(sqrt(tau2)) -log(sigma2) +
      log(psi(tau2/sigma2,n))
  }
  if(word==212){
    log.l <- -log(sqrt(sigma2)) -(3/2)*log(n*tau2 + sigma2)
  }
  if(word==221){
    log.l <- -log(sqrt(tau2))-log(sigma2)-.5*log(n*tau2+sigma2)+
      log(psi(tau2/sigma2,n))
  }
  if(word==122 | word==211| word==123| word==213| word==312){
    log.l <- -log(sigma2) -log(n*tau2 + sigma2)
  }
  if(word==132 | word==231| word==321){
    log.l <- -Cn(n)*log(sqrt(tau2)) -log(sigma2) - log(psi(tau2/sigma2,n))
  }
  return(log.l)
}

log.lik <- function(y,theta,alpha){
  k <- dim(y)[2]
  n <- dim(y)[1]
  means <- matrix(0,n,k)
  for(i in 1:k)  means[,i] <- rep(theta[1]+alpha[i],n)
  sum(dnorm(y,means,sqrt(theta[2]),log=TRUE))
}
  
pi.re <- function(alpha,theta){
  sum(dnorm(alpha,0,sqrt(theta[3]),log=TRUE))
}
  

update.re <- function(theta,alpha,y,sd.a){
 mu <- theta[1]
 sigma2 <- theta[2]
 tau2 <- theta[3]
 vars <- (1/sigma2+1/tau2)^(-1)
 means <- ((y-mu)/sigma2)*vars
 alpha <- rnorm(length(alpha),means,sqrt(vars))
  return(alpha)
}




mcmc_alg7<-function(iter,thin,y,theta,c,phi,alpha,re,sd.t,sd.a){
  n.y <- length(y)
  theta.sims<-matrix(0,iter/thin,3)
  c.sims<-matrix(0,iter/thin,3)
  phi.sims<-matrix(0,iter/thin,10)
  re.sims <- matrix(,iter/thin,n.y)
  for(i in 1:iter){
#      print(i/iter)
      theta<-update.theta(theta,re,c,y,sd.t,phi)
      trans<-update.c(theta,c,phi,y,alpha,re)
      c<-trans$c
      phi<-trans$phi
      trans<-update.phi(theta,c,phi,y)
      c<-trans$c
      phi<-trans$phi
      re <- update.re(theta,re,y,sd.a)
      if(i%%thin==0){
      print(i/iter)
	theta.sims[i/thin,]<-theta
      	c.sims[i/thin,]<-c
        #print(phi)
      	phi.sims[i/thin,1:length(phi)]<-phi
        re.sims[i/thin,] <- re
      } 
    }
  list(theta.sims=theta.sims,c.sims=c.sims,phi.sims=phi.sims)
}


set.seed(10)

n <- 10
k <- 3

mu.true <- 2
tau2.true <- 2^2
sigma2.true <- 4^2
alpha.true <- rnorm(k,0,sqrt(tau2.true))
theta.true <- c(mu.true,sigma2.true,tau2.true)


y <- matrix(0,n,k)
for(i in 1:k){
  y[,i] <- rnorm(n,mu.true+alpha.true[i],sqrt(sigma2.true))
}

iter <- 10000000
thin <- 10

re <- alpha.true
theta <- theta.true
phi <- c(1,2,3)
alpha <- 100

sd.t <- c(2,9,10)
sd.a <- rep(1,n)

c <- rep(1,3)

draws <- mcmc_alg7(iter,thin,y,theta,c,phi,alpha,re,sd.t,sd.a)

draws$c


save.image("re2.RData")

############################################
# Post
############################################

words.sims <- rep("",iter/thin)

for(t in 1:(iter/thin)){
  words.sims[t] <- covernt3words(draws$c[t,],draws$phi[t,])
}

plot(words.sims)

unique(words.sims)

pdf("norm_barplot.pdf",height=5,width=5)
barplot(summary(as.factor(words.sims))[order(summary(as.factor(words.sims)),decreasing=TRUE)]/(iter/thin),
        ylab="probability",xlab=expression(phi))
dev.off()

#plot(as.factor(words.sims))
