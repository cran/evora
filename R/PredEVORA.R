PredEVORA <-
function(evora.o,test.m,testCOPA.m,phenoTEST.v){

library(Hmisc);
  
selCpG.idx <- evora.o$riskCpG;
match(rownames(evora.o$data)[selCpG.idx],rownames(test.m)) -> map.idx;
map.idx <- map.idx[which(is.na(map.idx)==FALSE)];
print(paste("The number of risk features not found in test data: nNA=",length(which(is.na(map.idx))),sep=""));
tmp.m <- test.m[map.idx,];
tmpmean.v <- apply(tmp.m,2,mean);
copa.m <- testCOPA.m[map.idx,];
copath <- evora.o$th;

noutpSpth.v <- vector();
for(s in 1:ncol(testCOPA.m)){
 noutpSpth.v[s] <- length(which(copa.m[,s]>copath));###counts outlier CpGs
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

evora.v <- vector();
rc.o <- rcorr.cens(noutpSpth.v,phenoTEST.v);
evora.v[1] <- 0.5*(rc.o[[2]]-1.96*rc.o[[3]]+1);  
evora.v[2] <- rc.o[[1]];
evora.v[3] <- 0.5*(rc.o[[2]]+1.96*rc.o[[3]]+1);

return(list(auc=evora.v,riskS=noutpSpth.v));
 
} ### END 

