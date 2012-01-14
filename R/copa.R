copa <-
function(data.m,lowercent=0.8,overexp.log=TRUE){


  ## auxilliary function for sorting with NA's
  sortNA <- function(v,decreasing=TRUE){
  which(is.na(v)==FALSE) -> true.idx;
  length(v) -> ns;
  if( decreasing ){
   tmpD.s <- sort(v,decreasing=TRUE,index.return=TRUE);
   vals.v <- tmpD.s$x;
   index.v <- true.idx[tmpD.s$ix];
  }
  if( !decreasing ){
   tmpI.s <- sort(v,decreasing=FALSE,index.return=TRUE);
   vals.v <- tmpI.s$x;
   index.v <- true.idx[tmpI.s$ix];
  }
  vals.v <- c(vals.v,rep(NA,ns-length(true.idx)));
  index.v <- c(index.v,rep(NA,ns-length(true.idx)));
  return(list(x=vals.v,ix=index.v));
  }

  Ns <- ncol(data.m);
  # do copa transformation
  copa.m <- data.m
  copaS.m <- matrix(nrow=nrow(data.m),ncol=ncol(data.m));
  print("Computing COPA transformation");
  for( g in 1:nrow(data.m)){
  tmp.v <- data.m[g,] - median(data.m[g,],na.rm=TRUE);
  tmp.v[c(which(tmp.v > 0),which(tmp.v < 0))] -> tmp2.v;
  mad <- median(abs(tmp2.v),na.rm=TRUE);
  copa.m[g,] <- as.numeric(tmp.v/mad);
  tmp.s <- sortNA(as.vector(as.numeric(copa.m[g,])),decreasing=overexp.log);
  copaS.m[g,] <- tmp.s$x ;
  } 
  rownames(copaS.m) <- rownames(copa.m);
  print("Chcckpt1");
  # now pick percentiles and rank genes according to their copa scores
  centiles.v <- 1- c(1:Ns)/Ns ;
  NtopS <- min(which(centiles.v <= lowercent));
  copaOUT.lm <- list();
  for( ns in 1:NtopS ){
   tmp.v <- as.vector(as.numeric(copaS.m[,ns]));
   tmp.s <- sortNA(tmp.v,decreasing=overexp.log);
   copaOUT.lm[[ns]] <- cbind(rownames(data.m)[tmp.s$ix],tmp.s$x);
   colnames(copaOUT.lm[[ns]]) <- c("Gene","CopaScore");
   print(paste("Doing for centile ",centiles.v[ns],sep=""));
  }
  names(copaOUT.lm) <- 1:NtopS;

  return(list(copa=copa.m,out=copaOUT.lm,centiles=centiles.v));

}

