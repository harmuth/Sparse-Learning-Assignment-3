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

dim(X.mean.matrix)
dim(class.col.means)
X.residuals <- x-X.mean.matrix

X.residuals.svd <- svd(X.residuals)

pc <- t(X.residuals.svd$d * t(X.residuals.svd$u)) %>%
  as.data.frame() %>% 
  bind_cols(data.frame(class = classes))
  
ggplot(pc, aes(V1,V2,col=class))+geom_point()


