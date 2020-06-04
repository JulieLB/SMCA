data(cheese)

G1 <- partition_select (cheese [,3])
G2 <- partition_select (cheese [,1])
G3 <- partition_select (cheese [,2])
Gvar <- partition_variables(cheese)

G <- list()
for (i in 1:length(G1)) {
  for (j in 1:length(G2)) {
    for (k in 1:length(G3)){
      name <- paste(names(G1)[i], "_", names(G2)[j], "_", names(G3)[k], sep = "")
      G$new <- G1[[i]][G1[[i]] %in% G2[[j]]]
      names(G)[length(G)] <- name
    }
  }
}


##parameters
#row constraint
i <- nrow(cheese)
c1 <- c(1, 0.5*sqrt(i), sqrt(i))

#columns constraints
j <- length(unlist(Gvar))
c2 <- c(1, 0.5*sqrt(j), sqrt(j))



time_compare <- data.frame(noG = rep(0, 3) ,
                           Gcity = rep(0, 3),
                           Gcity.sex.age = rep(0, 3),
                           Gvar = rep(0,3),
                           Gcity_Gvar = rep(0, 3),
                           Gcity.sex.age_Gvar = rep(0, 3), row.names = c("high", "medium", "none"))

##########comparing the time
## without sparsity
tic("analysis 0")
res1 <- SMCA(cheese, c1 = c1[3], c2 = c2[3], Grow = G, Gcol = Gvar, n = 10, init="svd")
t <- toc()

time_compare[3,6] <- t$toc - t$tic

tic("analysis 0")
res1 <- SMCA(cheese, c1 = c1[3], c2 = c2[3], Grow = G1, Gcol = Gvar, n = 10, init="svd")
t <- toc()

time_compare[3,5] <- t$toc - t$tic
tic("analysis 0")
res1 <- SMCA(cheese, c1 = c1[3], c2 = c2[3], Gcol = Gvar, n = 10, init="svd")
t <- toc()

time_compare[3,4] <- t$toc - t$tic

tic("analysis 0")
res1 <- SMCA(cheese, c1 = c1[3], c2 = c2[3], Grow = G, n = 10, init="svd")
t <- toc()

time_compare[3,3] <- t$toc - t$tic
tic("analysis 0")
res1 <- SMCA(cheese, c1 = c1[3], c2 = c2[3], Grow = G1, n = 10, init="svd")
t <- toc()

time_compare[3,2] <- t$toc - t$tic

tic("analysis 0")
res1 <- SMCA(cheese, c1 = c1[3], c2 = c2[3], n = 10, init="svd")
t <- toc()

time_compare[3,1] <- t$toc - t$tic

## high sparsity
tic("analysis 1")
res1 <- SMCA(cheese, c1 = c1[1], c2 = c2[1], n = 10, init="svd")
t <- toc()
time_compare[1,1] <- t$toc - t$tic

tic("analysis 1")
res1 <- SMCA(cheese, c1 = c1[1], c2 = c2[1], Grow = G1, n = 10, init="svd")
t <- toc()
time_compare[1,2] <- t$toc - t$tic

tic("analysis 1")
res1 <- SMCA(cheese, c1 = c1[1], c2 = c2[1], Grow = G, n = 10, init="svd")
t <- toc()
time_compare[1,3] <- t$toc - t$tic

tic("analysis 1")
res1 <- SMCA(cheese, c1 = c1[1], c2 = c2[1], Gcol = Gvar, n = 10, init="svd")
t <- toc()
time_compare[1,4] <- t$toc - t$tic

tic("analysis 1")
res1 <- SMCA(cheese, c1 = c1[1], c2 = c2[1], Grow = G1, Gcol = Gvar, n = 10, init="svd")
t <- toc()
time_compare[1,5] <- t$toc - t$tic

tic("analysis 1")
res1 <- SMCA(cheese, c1 = c1[1], c2 = c2[1], Grow = G, Gcol = Gvar, n = 10, init="svd")
t <- toc()
time_compare[1,6] <- t$toc - t$tic

## medium sparsity
tic("analysis 2")
res1 <- SMCA(cheese, c1 = c1[2], c2 = c2[2], n = 10, init="svd")
t <- toc()
time_compare[2,1] <- t$toc - t$tic

tic("analysis 2")
res1 <- SMCA(cheese, c1 = c1[2], c2 = c2[2], Grow = G1, n = 10, init="svd")
t <- toc()
time_compare[2,2] <- t$toc - t$tic

tic("analysis 2")
res1 <- SMCA(cheese, c1 = c1[2], c2 = c2[2], Grow = G, n = 10, init="svd")
t <- toc()
time_compare[2,3] <- t$toc - t$tic

tic("analysis 2")
res1 <- SMCA(cheese, c1 = c1[2], c2 = c2[2], Gcol = Gvar, n = 10, init="svd")
t <- toc()
time_compare[2,4] <- t$toc - t$tic

tic("analysis 2")
res1 <- SMCA(cheese, c1 = c1[2], c2 = c2[2], Grow = G1, Gcol = Gvar, n = 10, init="svd")
t <- toc()
time_compare[2,5] <- t$toc - t$tic

tic("analysis 2")
res1 <- SMCA(cheese, c1 = c1[2], c2 = c2[2], Grow = G, Gcol = Gvar, n = 10, init="svd")
t <- toc()
time_compare[2,6] <- t$toc - t$tic




###
time_compare


