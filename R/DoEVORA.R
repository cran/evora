DoEVORA <-
function(train.m,copa.m,pheno.v,copath.v=1:10,ntop.v=seq(100,1000,50),nup=5000){

library(Hmisc);
library(qvalue);

### construct 10Fold CV sets
print("Constructing 10 fold sets");
nfolds <- 10;
train10Fart.li <- list();
test10Fart.li <- list();

nspg.v <- summary(factor(pheno.v));
nstest.v <- floor(nspg.v/nfolds);

rem.idx <- 1:length(pheno.v);
for(f in 1:(nfolds-1)){
  rem1.idx <- rem.idx[which(pheno.v[rem.idx]==1)];
  rem2.idx <- rem.idx[which(pheno.v[rem.idx]==2)];  
  test10Fart.li[[f]] <- c(sample(rem1.idx,nstest.v[1],replace=FALSE),sample(rem2.idx,nstest.v[2],replace=FALSE));
  train10Fart.li[[f]] <- setdiff(1:length(pheno.v),test10Fart.li[[f]]);
  rem.idx <- setdiff(rem.idx,test10Fart.li[[f]]);
}
f <- nfolds;
test10Fart.li[[f]] <- setdiff(1:length(pheno.v),unlist(test10Fart.li[1:9]));
train10Fart.li[[f]] <- setdiff(1:length(pheno.v),test10Fart.li[[f]]);

print("Ranking all features");
c.idx <- which(pheno.v==2);
h.idx <- which(pheno.v==1);
tmp.m <- matrix(nrow=nrow(train.m),ncol=2);
for(r in 1:nrow(train.m)){
 vt.o <- bartlett.test(train.m[r,] ~ pheno.v);
 lratio <- log(var(train.m[r,c.idx])/var(train.m[r,h.idx]));
 tmp.m[r,] <- c(lratio,vt.o$p.value);
}
tmp.s <- sort(tmp.m[,1],decreasing=TRUE,index.return=TRUE);
qv.v <- qvalue(tmp.m[,2])$qvalue;
sqv.v <- sort(qv.v,decreasing=FALSE);
nfdr <- length(which(sqv.v<0.05));
print(paste("Number of features differentially variable at FDR<0.05 is ",nfdr,sep=""));
if ( nfdr < 10 ){
  print("This number is far too small for EVORA to be practical!");
}

### do ranking for each fold
print("Doing folds");
selCpG.idx <- tmp.s$ix[1:nup];
rankedCpG.idx <- tmp.s$ix;

dvc.li <- list();
for(f in 1:nfolds){
 c.idx <- intersect(train10Fart.li[[f]],which(pheno.v==2));
 h.idx <- intersect(train10Fart.li[[f]],which(pheno.v==1)); 
 tmp.m <- matrix(nrow=length(selCpG.idx),ncol=2);
 for(r in 1:length(selCpG.idx)){
  vt.o <- bartlett.test(train.m[selCpG.idx[r],train10Fart.li[[f]]] ~ pheno.v[train10Fart.li[[f]]]);
 lratio <- log(var(train.m[selCpG.idx[r],c.idx])/var(train.m[selCpG.idx[r],h.idx]));
  tmp.m[r,] <- c(lratio,vt.o$p.value);
 }
 tmp.s <- sort(tmp.m[,1],decreasing=TRUE,index.return=TRUE);
 dvc.li[[f]] <- selCpG.idx[tmp.s$ix];
 print(f);
}

tmpART.m <- matrix(nrow=length(ntop.v),ncol=length(copath.v));
colnames(tmpART.m) <- copath.v;
rownames(tmpART.m) <- ntop.v;

for(ni in 1:length(ntop.v)){
 
 for(copath in copath.v){

  allscores.v <- rep(NA,length(pheno.v));
  for(f in 1:nfolds){

    sel.idx <- dvc.li[[f]][1:ntop.v[ni]]; 
    tmp.m <- copa.m[sel.idx,test10Fart.li[[f]]];
    tmpmean.v <- apply(tmp.m,2,mean);
    noutpSpth.v <- vector();
    for(s in 1:ncol(tmp.m)){
     noutpSpth.v[s] <- length(which(tmp.m[,s]>copath));###counts outlier CpGs
    }

    nsprg.v <- summary(factor(noutpSpth.v));
    rg.v <- as.numeric(names(nsprg.v));
    scores.v <- vector();
    for(rg in 1:length(rg.v)){
     tmpS.idx <- which(noutpSpth.v==rg.v[rg]); ## select samples with same risk score
     if (length(tmpS.idx)>=2){ ### then need to discriminate
      tmp.v <- tmpmean.v[tmpS.idx];
      maxv <- max(tmp.v); minv <- min(tmp.v);
      tmpN.v <- (tmp.v-minv)/(maxv-minv); ## scale to 0 and 1
      epsilon <- 0.000001;
      tmpN.v[which.max(tmpN.v)] <- 1-epsilon;
      tmpN.v[which.min(tmpN.v)] <- epsilon;
      scores.v[tmpS.idx] <- rg.v[rg]+tmpN.v;
     }
     else {
      scores.v[tmpS.idx] <- rg.v[rg];
     }
    }
    noutpSpth.v <- scores.v/nrow(tmp.m);
    allscores.v[test10Fart.li[[f]]] <- noutpSpth.v;
  }
  if( length(which(is.na(allscores.v)))>0){ print("PROBLEM WITH ASSIGNED SCORES");}

  rc.o <- rcorr.cens(allscores.v,pheno.v);
  tmpART.m[ni,copath] <- rc.o[[1]];
 }
 print(ntop.v[ni]);
}

maxI <- max(which(tmpART.m==max(tmpART.m)));
colI <- ceiling(maxI/nrow(tmpART.m));
rowI <- as.integer(10*( maxI/nrow(tmpART.m)  - floor(maxI/nrow(tmpART.m)) ));
if(rowI==0){
  rowI <- nrow(tmpART.m);
}

optN <- ntop.v[rowI];
optCOPA <- copath.v[colI];

return(list(data=train.m,riskCpG=selCpG.idx[1:optN],th=optCOPA,ranked=rankedCpG.idx,auc=tmpART.m));

} ## END of FUNCTION

