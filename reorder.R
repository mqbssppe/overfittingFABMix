#first you have to set wd
setwd("q_2_run_1/")
library('label.switching')
library('RColorBrewer')
z <- as.matrix(read.table("zValues.txt"))
logl <- read.table("kValues.txt", header=T)
burn <- 1:1
logl <- logl[-burn, ]
z <- z[-burn,]
print(table(logl[,1]))
K <- as.numeric(names(sort(table(logl[,1]),decreasing=TRUE)[1]))
kSelected <- K
index <- which(logl[,1] == K)
Kindex <- index
logl <- logl[index,2]
z <- z[index,]
mapIndex <- which(logl == max(logl))[1]
zPivot <- as.numeric(z[mapIndex,])

K <- max(z)
m <- length(logl)
ls <- label.switching(method = c("ECR", "ECR-ITERATIVE-1"), zpivot = zPivot, z = z, K = K)
oldLabels <- 1:K
index <- Kindex

allocationsECR <- z
for (i in 1:m){
        myPerm <- order(ls$permutations$"ECR"[i,])
        allocationsECR[i,] <- myPerm[z[i,]]
}

l <- as.matrix(read.table(paste0("LambdaValues",1,".txt")))
J <- dim(l)[2]
mcmc <- array(data = NA, dim = c(m,K,J))
for(k in 1:K){
#       lMean <- apply(l,2,mean)
#       lMean.matrix <- matrix(lMean,nrow = p, ncol = q, byrow=TRUE)
        l <- as.matrix(read.table(paste0("LambdaValues",k,".txt")))
        l <- l[-burn,]
        mcmc[,k,] <- l[index,]
}
lambda.perm.mcmc <- permute.mcmc(mcmc, ls$permutations$ECR)$output
lambda.mean <- array(data = NA, dim = c(K,p,q))
lambda.map <- array(data = NA, dim = c(K,p,q))
for(k in 1:K){
        lMean <- apply(lambda.perm.mcmc[,k,],2,mean)
        lambda.mean[k,,] <- matrix(lMean,nrow = p, ncol = q, byrow=TRUE)
        lambda.map[k,,] <- matrix(lambda.perm.mcmc[mapIndex,k,],nrow = p, ncol = q, byrow=TRUE)
}

mu <- read.table("muValues.txt")[max(burn) + Kindex,] # auto to grafei ws mu_{11},mu_{12},...,mu_{1K}, ...., mu_{p1},mu_{p2},...,mu_{pK} gia kathe grammi
mu.mcmc <- array(data = NA, dim = c(m,K,p))
for(k in 1:K){
        mu.mcmc[,k,] <- as.matrix(mu[,k + K*((1:p)-1)])
}
mu.mcmc <- permute.mcmc(mu.mcmc, ls$permutations$ECR)$output
mu.mean <- array(data = NA, dim = c(K,p))
mu.map <- array(data = NA, dim = c(K,p))
for(k in 1:K){
        for(j in 1:p){
                mu.mean[k,j] <- mean(mu.mcmc[,k,j])
                mu.map[k,j] <- mu.mcmc[mapIndex,k,j]
        }
}



w.mcmc <- array(as.matrix(read.table("wValues.txt")[max(burn) + Kindex, ]),dim = c(length(Kindex),K,1))
w.mcmc_raw <- w.mcmc
w.mcmc <- array(w.mcmc[,1:K,],dim = c(length(Kindex),K,1))
w.mcmc <- permute.mcmc(w.mcmc, ls$permutations$"ECR")$output
w.mean <- numeric(K)
w.map <- numeric(K)
for(k in 1:K){
        w.mean[k] <- mean(w.mcmc[,k,1])
        w.map[k] <- w.mcmc[mapIndex,k,1]
}

matplot(w.mcmc_raw[,,1], type = "p", pch = 16)
matplot(w.mcmc[,as.numeric(names(table(zPivot))),1], type = "p", pch = 16)



nonTrivialClusters <- as.numeric(names(which(table(ls$clusters[1,]) > 0)))
par(mfrow = c(2,2))
for(j in 1:q){
	matplot(lambda.mean[nonTrivialClusters,,j],main = paste0("Factor ",j))
}
par(mfrow = c(2,2))
for(j in 1:q){
	matplot(t(lambda.mean[nonTrivialClusters,,j]),main = paste0("Factor ",j),type = "l",lty = 1,xlim = c(1,p+1))
	rect(p,-10,p+10,10,col = "lightgray",border = NA)
	text(p + 0.7, t(lambda.mean[nonTrivialClusters,p,j]),paste0(nonTrivialClusters, ": ", as.numeric(table(ls$clusters[1,])[as.character(nonTrivialClusters)])),col = 1:length(nonTrivialClusters))
	abline(v = c(10.5,20.5, 30.5),lty = 2)
}


catNames <- 1:p

pdf(file = "regulonExpressionPerCluster.pdf",width = 16, height = 8)
colors <- brewer.pal(kSelected, name="Set1")
myCols <- colorRampPalette(colors)

perm <- order(table(ls$clusters[1,]))
par(mfrow = c(1,2),mar = c(3,2,2,2))
aa <- 0
aoua <- 0
myFinalIndex <- as.numeric(names(which(sort(table(ls$clusters[1,])) > 0)))
nColors <- length(myFinalIndex)
for(k in as.numeric(names(sort(table(ls$clusters[1,]))))){
        aa <-  aa + 1
        aoua <- aoua + 1
	if( length(which(ls$clusters[1,] == k)) > 0){

		#plot(c(0,p+1),c(-1,1),type = "n",xaxt="n", main = paste0("cluster \'", k, "\': ", as.numeric(table(ls$clusters[1,]))[perm][aa]," genes"))
		keepPoints <- c()
		for(j in 1:q){
		        myIndex <- j + (1:p)*q - q
		        tmp <- rep(0,p)
		        m1 <- 0
		        for(iter in 1:m){
		                myVec <- perm.mcmc[iter,k,myIndex]
		                if(is.nan(sum(myVec)) == FALSE){
		                        m1 <- m1 + 1;tmp <- tmp + myVec
					#if(iter%%20 == 0 ){
					#	points(1:p,myVec,col = j,type = "l")
					#}
		                }
		        }
		        tmp <- tmp/m1
			keepPoints <- rbind(keepPoints,tmp)
		        #points(1:p,tmp,type = "b",col = j,lwd = 2)
			#text(p+0.5,tmp[p],labels = as.character(j),col = myCols(kSelected)[perm[aoua]],lwd = 2)
		}
		matplot(t(keepPoints),type = "b",col = myCols(nColors),xaxt="n", main = paste0("cluster \'", k, "\': ", as.numeric(table(ls$clusters[1,]))[perm][aa]," guys"))
		abline(v = c(20.5),lty = 2,col = "gray")
		axis(1, at= 1:p,labels=catNames)

	}
}
dev.off()





################################################################################################
################################################################################################
################################################################################################

write.table(ls$clusters[1,], file = "singleBestClusterings.txt", quote = FALSE, row.names = FALSE, col.names = "clusterLabel")
zMAP <- read.table("singleBestClusterings.txt",header=TRUE)[,1]



mapAllocationsPosteriorProbs <- numeric(n)
for(i in 1:n){
        mapAllocationsPosteriorProbs[i] <- length(which(allocationsECR[,i] == zMAP[i]))
}
mapAllocationsPosteriorProbs <- mapAllocationsPosteriorProbs/dim(allocationsECR)[1]
subIndex <- which(mapAllocationsPosteriorProbs < 0.9)

pdf(file = "allocationPosteriorProbs.pdf",width = 16, height = 8)
par(mfrow = c(2,2),mar = c(4,4,3,2))
for(k in myFinalIndex){
	index <-which(zMAP == k)
	if( abs(diff(range(mapAllocationsPosteriorProbs[index]))) > 0 ){
		hist(mapAllocationsPosteriorProbs[index],main = paste0("cluster `",k, "` (",length(index)," guys)"), col = 'gray', xlab = 'MAP classification probability')
	}else{
		plot(c(0,1), c(0,length(index)), type = "n", main = paste0("cluster `",k, "` (",length(index)," guys)"), xlab = 'MAP classification probability', ylab = "Frequency")
		rect(mapAllocationsPosteriorProbs[index][1] - 0.025, 0, mapAllocationsPosteriorProbs[index][1], length(index), col = 'gray')
	}
}
dev.off()

# if you enable this, then you filter out observations with lower MAP probs
#zMAP[subIndex] = rep(0,length(subIndex))
pdf(file = "clusterEmpiricalProfiles.pdf",width = 16, height = 8)
par(mfrow = c(2,2),mar = c(3,3,2,2))
aa <- 0
aoua <- 0
keepMean <- c()
for(k in myFinalIndex){
        print(paste("cluster ",k,", n-guys = ",length(which(zMAP == k))))
        aa <-  aa + 1
        aoua <- aoua + 1
   #     plot(c(0,p+1),range(originalX),type = "n", main = paste0("cluster \'", k, "\': ", as.numeric(table(zMAP))[perm][aa]," genes"))
	index <- intersect(which(zMAP == k),which(mapAllocationsPosteriorProbs > 0))
	if(length(index) == 0){index<-which(zMAP == k);cat(paste0("warning: 0 observations with high prob"),"\n")}
#	index <- which(zMAP == k)
        matplot(t(array(originalX[index,], dim = c( length(index), dim(originalX)[2] ))),type = "l",col = "gray",main = paste0("cluster \'", k, "\': ", length(index)," genes (filtered > 0.9)"),xaxt="n")
        myMean <- colSums(array(originalX[index,], dim = c( length(index), dim(originalX)[2] )))/length(index)
	keepMean <- rbind(keepMean,myMean)
        points(1:p,myMean,type = "l",col = myCols(nColors)[aa],lwd = 4)
        points(1:p,myMean,type = "l",col = "white",lwd = 1)
	points(1:p,myMean,pch = 16,col = myCols(nColors)[aa])
	abline(v = c(10.5,20.5,30.5),lty = 2)
	axis(1, at= 1:p,labels=catNames)
}
dev.off()

matplot(t(keepMean),type = "l",pch=16,lty=1,col = myCols(kSelected)[perm],main = "means",xaxt = "n")
matpoints(t(keepMean),pch=16,lty=1,col = myCols(kSelected)[perm])
abline(v = c(3.5,6.5,9.5,15.5,21.5),lty = 2)
axis(1, at= c(2,5,8,12.5,18.5,23),labels=catNames)




colors <- brewer.pal(10, name="Set1")
myCols <- colorRampPalette(colors)


pdf(file = "factor-scores.pdf",width = 16, height = 8)
par(mfrow = c(2,2))
for(j in 1:q){
	matplot(t(lambda.mean[myFinalIndex,,j]),main = paste0("Factor ",j),type = "n",lty = 1,xlim = c(1,p+1),xaxt="n")
	#rect(0,-10,p+10,10,col = "gray90",border = NA)
	matplot(t(lambda.mean[myFinalIndex,,j]),main = paste0("Factor ",j),type = "l",lty = 1,xlim = c(1,p+1),col=myCols(myFinalIndex),add = TRUE)
#	matplot(t(lambda.mean[myFinalIndex,,j]),main = paste0("Factor ",j),type = "p",lty = 1,xlim = c(1,p+1),col=lane,add = TRUE)
	text(p + 0.7, t(lambda.mean[myFinalIndex,p,j]),paste0(myFinalIndex, ": ", as.numeric(table(zMAP)[as.character(myFinalIndex)])),col = myCols(myFinalIndex))
#	abline(v = c(3.5,6.5,9.5,12.5,15.5),lty = 2)
	axis(1, at= 1:p,labels=catNames)
}
dev.off()

j <- j+1;matplot(lambda.mean[myFinalIndex,,j],main = paste0("Factor ",j),type = "b",lty = 1,col = mycol,pch = mypch)

pdf(file = "barplot.pdf",width = 16,height = 8)
par(mfrow = c(2,2))
for(j in 1:q){
	barplot(lambda.mean[myFinalIndex,,j],beside=T,main = paste0("Factor ",j))
	abline(v = 9*c(3,6,9,15,21)+0.5,lty = 2)
}
dev.off()

mypch = c(rep(16, 9), rep(2, 3), rep(3, 6))
mycol = c("darkorange","darkorange","darkorange",2,2,2,3,3,3,4,4,4,5,5,5,6,6,6) 

mypch = 16
mycol = c(rep('blue',10), rep('red',10), rep('darkorange',10), rep('darkgreen',10))



par(mfrow = c(2,5))
for(k in 1:10){
	plot(lambda.mean[nonTrivialClusters[k],,1],lane,col = mycol , pch = mypch)
}

##### get Y
yTable <- read.table("yValues.txt")
yTable <- yTable[-burn,]
yColSums <- colSums(yTable)
m <- dim(yTable)[1]
y.estimate <- array(data = 0, dim = c(n,q))
for(i in 1:n){
        y.estimate[i,] <- as.numeric(yColSums[seq(i,q*n,by = n)])
}
y.estimate <- y.estimate/m

i<- i +1;matplot(yTable[,(1:q - 1)*n + i]);abline(h = y.estimate[i,],col = 1:q);mapAllocationsPosteriorProbs[i]
# c(i,n+i,2*n+i,3*n+i)


my_line <- function(x,y){
    points(x,y, col = mycol , pch = mypch,cex = 1.5,)
    abline(h = 0,v = 0, lty = 2, col = 'gray')
}

for(k in myFinalIndex){
#	pdf(file = paste0("cluster-",k,".pdf"),width = 16,height = 8)
	pairs(lambda.mean[k,,],lower.panel = my_line, upper.panel = my_line, labels = paste0('factor ',1:q, ', k = ', k, ' (', length(which(zMAP == k)),')'))
	
	cat ("Press [enter] to continue")
	line <- readline()
#	dev.off()
}


par(mfrow = c(6,q-1),mar = c(4,4,2,2))
for(k in as.numeric(names(sort(table(zMAP))))){
	if(length(which(zMAP==k))>26){
	for(i in 2:q){
		plot(lambda.mean[k,,c(i,1)],col = mycol , pch = mypch,cex = 1.5,main = paste0("cluster ",k),xlab = "tf 1",ylab = paste0("tf ",i))
	}
	}
}

geneNames <- read.table("../log.expression.sybaris8.txt",header=TRUE,row.names=NULL)[,1]
x_data <- read.table("../log.expression.sybaris8.txt",header=TRUE,row.names=NULL)
n <- dim(x_data)[1]
p <- dim(x_data)[2]
x_data <- x_data[,-1]
rowMeans <- apply(x_data,1,mean)
rS <- which(rowMeans < 1)
x_data <- as.matrix(x_data[-rS,])
n <- dim(x_data)[1]
p <- dim(x_data)[2]
geneNames <- geneNames[-rS]


write.table(cbind(geneNames,zMAP,mapAllocationsPosteriorProbs),file = paste0("clustersOfGenes.txt"),quote = FALSE,row.names=FALSE,col.names = c("geneID","clusterID","prob"))

par(mfrow = c(3,3))
for(k in as.numeric(names(sort(table(zMAP))))){
	if(k > 0){
#	pdf(file = paste0("cluster-",k,".pdf"),width = 9,height = 6)
		matplot(t(y.estimate[which(zMAP == k),]),type = "l",col = "gray",main = paste0("cluster \'", k, "\': ", length(which(zMAP == k))," genes"))
		myMean <- colSums(y.estimate[which(zMAP == k),])/length(which(zMAP == k))
		keepMean <- rbind(keepMean,myMean)
		points(1:q,myMean,type = "l",col = myCols(kSelected)[perm[aoua]],lwd = 4)
		points(1:q,myMean,type = "l",col = "white",lwd = 1)
		points(1:q,myMean,pch = aa,col = myCols(kSelected)[perm[aoua]])
#	dev.off()
	}
}












