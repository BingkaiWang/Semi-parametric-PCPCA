set.seed(123)
p <-4
k <- 2
n <- 3
lambda1 <- 1:4
lambda2 <- sample(1:4, 4, replace = F)
lambda3 <- sample(1:4, 4, replace = F)
lambda <- rbind(lambda1, lambda2, lambda3)
gamma_true <- svd(matrix(rnorm(p^2), p, p))$u[,1:k]

# CPCs
png(filename="demonstration/cpc1.png")
heatmap(gamma_true[,1] %*% t(gamma_true[,1]), Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()

png(filename="demonstration/cpc2.png")
heatmap(gamma_true[,2] %*% t(gamma_true[,2]), Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()

# subject 1
nCPC1 <- gramSchmidt(cbind(gamma_true, matrix(rnorm(p*(p-k)), p, p-k)))$Q[,(k+1):p]
png(filename="demonstration/ncpc1-1.png")
heatmap(nCPC1[,1] %*% t(nCPC1[,1]), Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()
png(filename="demonstration/ncpc1-2.png")
heatmap(nCPC1[,2] %*% t(nCPC1[,2]), Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()
png(filename="demonstration/ncpc1.png")
heatmap(nCPC1 %*% diag(lambda1[(k+1):p], nrow = p-k) %*% t(nCPC1), Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()
cov1 <- gamma_true %*% diag(lambda1[1:k], nrow = k) %*% t(gamma_true) +
  nCPC1 %*% diag(lambda1[(k+1):p], nrow = p-k) %*% t(nCPC1)
png(filename="demonstration/cov1.png")
heatmap(cov1, Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()

# subject 2
nCPC2 <- gramSchmidt(cbind(gamma_true, matrix(rnorm(p*(p-k)), p, p-k)))$Q[,(k+1):p]
png(filename="demonstration/ncpc2-1.png")
heatmap(nCPC2[,1] %*% t(nCPC2[,1]), Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()
png(filename="demonstration/ncpc2-2.png")
heatmap(nCPC2[,2] %*% t(nCPC2[,2]), Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()
png(filename="demonstration/nCPC2.png")
heatmap(nCPC2 %*% diag(lambda2[(k+1):p], nrow = p-k) %*% t(nCPC2), Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()
cov2 <- gamma_true %*% diag(lambda2[1:k], nrow = k) %*% t(gamma_true) +
  nCPC2 %*% diag(lambda2[(k+1):p], nrow = p-k) %*% t(nCPC2)
png(filename="demonstration/cov2.png")
heatmap(cov2, Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()

# subject 3
nCPC3 <- gramSchmidt(cbind(gamma_true, matrix(rnorm(p*(p-k)), p, p-k)))$Q[,(k+1):p]
png(filename="demonstration/nCPC3-1.png")
heatmap(nCPC3[,1] %*% t(nCPC3[,1]), Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()
png(filename="demonstration/nCPC3-2.png")
heatmap(nCPC3[,2] %*% t(nCPC3[,2]), Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()
png(filename="demonstration/nCPC3.png")
heatmap(nCPC3 %*% diag(lambda3[(k+1):p], nrow = p-k) %*% t(nCPC3), Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()
cov3 <- gamma_true %*% diag(lambda3[1:k], nrow = k) %*% t(gamma_true) +
  nCPC3 %*% diag(lambda3[(k+1):p], nrow = p-k) %*% t(nCPC3)
png(filename="demonstration/cov3.png")
heatmap(cov3, Colv = NA, Rowv = NA, scale = "none", revC = T, cexRow = NA, cexCol = NA)
dev.off()


