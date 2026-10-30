# Summarize simulation results or compare them with the saved summary.
options(stringsAsFactors=FALSE)
invisible(try(suppressWarnings(Sys.setlocale('LC_ALL','Chinese (Simplified)_China.utf8')),silent=TRUE))
script<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/')
data_dir<-file.path(root,'results')
dest<-data_dir
mode<-if('--collect'%in%commandArgs(TRUE))'collect'else'verify'
groups<-read.csv(file.path(root,'code/settings.csv'))
saved_path<-file.path(data_dir,'summary.csv')
methods<-c('lr','ols','dof','loo','cloo','ploo','oracle');rows<-list()
for(i in seq_len(nrow(groups))) {
  x<-readRDS(file.path(data_dir,'raw',paste0(groups$id[i],'.rds')))
  g<-x$group;n<-g$n;B<-dim(x$raw)[1]
  for(p in dimnames(x$raw)[[2]])for(m in methods) {
    u<-x$raw[,p,paste0('u_',m)];v<-x$raw[,p,paste0('v_',m)]
    ok<-is.finite(u)&is.finite(v)&v>0;fu<-is.finite(u);fv<-is.finite(v)
    z<-sqrt(n)*u[ok]/sqrt(v[ok]);nn<-length(z);rate<-mean(z<qnorm(.025));cv<-qnorm(.975)
    den<-1+cv^2/nn;mid<-(rate+cv^2/(2*nn))/den;half<-cv*sqrt(rate*(1-rate)/nn+cv^2/(4*nn^2))/den
    rows[[length(rows)+1L]]<-data.frame(group=groups$id[i],n=n,design=g$design,scope=g$scope,ratio=g$ratio,experiment=g$experiment,p=as.integer(p),method=m,
      B=B,valid=nn,invalid=B-nn,nonpositive_v=sum(fv&v<=0),nonfinite_v=sum(!fv),nonfinite_u=sum(!fu),rejects=sum(z<qnorm(.025)),
      rate=rate,mcse=sqrt(rate*(1-rate)/nn),wilson_lower=mid-half,wilson_upper=mid+half,
      mean_rootn_u=sqrt(n)*mean(u[fu]),empirical_v=n*var(u[fu]),mean_v=mean(v[fv]),mean_v_ratio=mean(v[fv])/(n*var(u[fu])),
      median_v=median(v[fv]),q99_v=unname(quantile(v[fv],.99)),min_v=min(v[fv]),max_v=max(v[fv]),z_sd=sd(z))
  }
}
actual<-do.call(rbind,rows);stopifnot(nrow(actual)==700L)
if(mode=='collect') {
  write.csv(actual,saved_path,row.names=FALSE)
  cat('COLLECTED',nrow(actual),'method settings.\n')
}else {
  saved<-read.csv(saved_path)
key<-function(x)paste(x$group,x$p,x$method,sep='|')
stopifnot(!anyDuplicated(key(actual)),setequal(key(actual),key(saved)))
actual<-actual[match(key(saved),key(actual)),names(saved)]
checks<-lapply(names(saved),function(nm) {
  a<-actual[[nm]];b<-saved[[nm]]
  if(is.numeric(b)) {
    stopifnot(identical(is.na(a),is.na(b)),identical(is.finite(a),is.finite(b)))
    k<-is.finite(b);err<-if(any(k))max(abs(a[k]-b[k])/(1+abs(b[k])))else 0
    stopifnot(err<1e-12)
  }else {stopifnot(identical(as.character(a),as.character(b)));err<-0}
  data.frame(column=nm,max_relative_error=err,passed=TRUE)
})
  write.csv(do.call(rbind,checks),file.path(dest,'summary_check.csv'),row.names=FALSE)
  cat('PASS: 24 raw files, 700 settings, every summary column matches.\n')
}
