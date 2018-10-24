## Bastazini, V.A.G., Ferreira, P.M., Azambuja, B.O., Casas, G., Debastiani, V.J., Guimarães, P.R. & Pillar, V.D. 2017. Untangling the Tangled Bank: A Novel Method for Partitioning the Effects of Phylogenies and Traits on Ecological Networks. 
## Evolutionary Biology  44(3): 312–324. DOI: 10.1007/s11692-017-9409-8
## Available at github.com/bastazini/Networks
matrix.p1<-function (comm, dist.spp) {
  matrix.w <- sweep(comm, 1, rowSums(comm), "/")
  for(i in 1:dim(comm)[1]){
    for(j in 1:dim(comm)[2]){
      if(is.nan(matrix.w[i,j])==TRUE){
        matrix.w[i,j]=as.numeric(0)
      }
    }
  }
  similar.phy <- 1 - (dist.spp/max(dist.spp))
  matrix.phy <- 1/colSums(similar.phy)
  matrix.q <- sweep(similar.phy, 1, matrix.phy, "*")
  matrix.P <- matrix.w %*% matrix.q
  return(list(matrix.w = matrix.w, matrix.q = matrix.q, matrix.P = matrix.P))
}