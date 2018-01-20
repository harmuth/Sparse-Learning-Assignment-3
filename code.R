library(msgl)
data(PrimaryCancers)
ord_class <- order(classes)
ord_mir <- order(colMeans(x), decreasing = TRUE)
x <- x[ord_class, ord_mir]
classes <- classes[ord_class]
image(Matrix(x))

library(ggplot2)
library(dplyr)
library(tidyr)



class.col.means <- data.frame(class = classes) %>% 
  bind_cols(as.data.frame(x)) %>%
  group_by(class) %>%
  summarise_all(mean)


plot.cols <- class.col.means %>%
  gather(col,means,-class) %>%
  group_by(class) %>%
  mutate(col.num = row_number()) %>%
  arrange(class,col.num)


ggplot(plot.cols,aes(x=col.num,y=means)) +
  geom_line() + facet_wrap(~class)

X.mean.matrix <- data.frame(class = classes) %>% 
  left_join(class.col.means) %>% 
  select(-class)
X.residuals <- x-X.mean.matrix
X.residuals.svd <- svd(X.residuals)
pc <- t(X.residuals.svd$d * t(X.residuals.svd$u)) %>%
  as.data.frame() %>% 
  bind_cols(data.frame(class = classes))

ggplot(pc, aes(V1,V2,col=class))+geom_point()

class.pc <- data.frame(class = rep(unique(classes),each=length(pc$class)),
           target =  rep(pc$class,length(unique(classes))),
           pc1 = rep(pc$V1,length(unique(classes))),
           pc2 = rep(pc$V2,length(unique(classes)))) %>%
  mutate(same = class == target)

ggplot(class.pc, aes(pc1,pc2,col=same)) +
  geom_point() +
  facet_wrap(~class)



ORS.PCA <- function(X, t, max.iter, convergence.threshold) {
  over.threshold <- TRUE
  v <- matrix(1/dim(X)[2], ncol = 1, nrow = dim(X)[2])
  m <- 1
  u <- X %*% v / sqrt(sum((X %*% v)^2))
  v.f <- function (lambda, y) {
    sign(t(X) %*% y) * pmax(t(X) %*% y - lambda, 0)
  }
  while(over.threshold | m > max.iter) {
    new.u <- X %*% v / sqrt(sum((X %*% v) ^ 2))
    if(max(v.f(0, new.u)) < t){
      lambda <- 0
    } else {
      opt.func <- function(lambda){ abs(max(v.f(lambda, new.u)) - t)}
      lambda <- optimise(opt.func, c(0, max(v.f(0, new.u))))$minimum
    }
    
    new.v <- v.f(lambda, new.u)
    
    if((max(abs(new.v-v))<convergence.threshold) & (max(abs(new.u-u)) < convergence.threshold)){
      over.threshold = FALSE
    }
    cat("iter =",m,"\n","v convergence:",max(abs(new.v-v)),"\n","u convergence:",max(abs(new.u-u)),"\n", "lambda:", lambda,"\n")
    v <- new.v
    u <- new.u
    m <- m + 1
  }
  return(list(v = v, u = u, convergence = !over.threshold, n_iter = m - 1))
}

res <- ORS.PCA(as.matrix(X.residuals), 10, 100, 0.0001)

MRS.PCA <- function(X, k, t, max.iter, convergence.threshold) {
  
  # initial projection
  proj <- diag(nrow(X))
  
  # ready for output
  u <- matrix(NA, ncol = k, nrow = nrow(X))
  v <- matrix(NA, ncol = k, nrow = ncol(X))
  for(i in 1:k){
    new.X <- proj %*% X
    tmp.res <- ORS.PCA(new.X, t, max.iter, convergence.threshold)
    u[,i] <- tmp.res$u
    v[,i] <- tmp.res$v
    
    # update projection
    proj <- proj-u[,i] %*% t(u[,i])
  }
  
  return(list(v = v, u = u))
}

MRS.PCA(as.matrix(X.residuals), 2, 1000, 100, 0.0001)

