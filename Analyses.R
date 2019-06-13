## Resetting workspace
if(F){
  rm(list=ls(all=T))
  gc()
  cat("\014")
  ls(all=T)
}

## runder av til nøyaktig k desimaler
round2 <- function(x, k){format(round(x, k), nsmall=k)}

## Load .data. Takes about two minutes
## Data on CL, CP, CLP and CL_CLP for smoke, vitamin and drink
date();system.time({load("../data.RData")});date()


## Race
## 1-European, 2-African, 3-Asian, 4-Hispanic, 
## 5-Native Amer, 6-Malaysian, 7-Other

## Counting trios etc.
if(F){
## Race
## 1-European, 2-African, 3-Asian, 4-Hispanic, 
## 5-Native Amer, 6-Malaysian, 7-Other
all.data.cl=.data[["CL_CLP"]][["all"]]
## Table 1
tmp.cl=all.data.cl
tmp.cl=all.data.cl[phdata(all.data.cl)$proband_race==1]
tmp.cl=all.data.cl[phdata(all.data.cl)$proband_race==3]
tmp.cl=all.data.cl[phdata(all.data.cl)$proband_race%in%c(2,4:7)]
## All individuals
nids(tmp.cl)
##Missing=-9, not affected=1, affected (parents and children)=2
cc=table(phdata(tmp.cl)$cc);cc
## Number of individuals in a family
family.size=table(table(substr(idnames(tmp.cl),1,5)));family.size
sum(family.size)
## No of inds is OK
sum(as.numeric(names(family.size))*family.size)
sum(family.size["2"]*2)
sum(family.size["3"]*3)
sum((4:6)*family.size[c("4","5","6")],na.rm=T)

## Table 2
## De med ID som slutter på _1, _4, _5 eller _6 er barn
tmp.cl=all.data.cl
## Maternal characteristics, based on "first" child (_1, not _4, etc.)
child=grep("_1",idnames(tmp.cl));length(child)
# child=substr(idnames(tmp.cl),6,7)%in%paste("_",c(1,4:6),sep="");sum(child)
tmp.cl=tmp.cl[child]
## All
smoke=phdata(tmp.cl)$proband_smoke
table(smoke,useNA="ifany")
vit=phdata(tmp.cl)$proband_vitamin
table(vit,useNA="ifany")
drink=phdata(tmp.cl)$proband_drink
table(drink,useNA="ifany")
## Race
race=phdata(tmp.cl)$proband_race;table(race,useNA="ifany")
## Euro
tmp.race=tmp.cl[race==1]
smoke=phdata(tmp.race)$proband_smoke
table(smoke,useNA="ifany")
vit=phdata(tmp.race)$proband_vitamin
table(vit,useNA="ifany")
drink=phdata(tmp.race)$proband_drink
table(drink,useNA="ifany")
## Asia
tmp.race=tmp.cl[race==3]
smoke=phdata(tmp.race)$proband_smoke
table(smoke,useNA="ifany")
vit=phdata(tmp.race)$proband_vitamin
table(vit,useNA="ifany")
drink=phdata(tmp.race)$proband_drink
table(drink,useNA="ifany")

}

## With parent of origin.
## strata: 13=race, child; 16=smoke; 19=vitamin; 22=drink
system.time({

## Asia runs last and vitamin first, because code stops
## for Asia + smoke or drink
## cleft: CL, CP, CLP, CL_CLP
## environ: smoke, vitamin, drink
cleft="CL_CLP"
environ="drink";race="all";cleft;environ;race
for(race in c("all","eur","asia")){
  race2 <- 0
  if(race=="eur") race2=1
  if(race=="asia") race2=3
  for(environ in c("vitamin","smoke","drink")){
    if(race=="asia" & environ!="vitamin") break
    write(cleft,"cleft.txt")
    write(environ,"environ.txt")
    .tmp=.data[[cleft]][[environ]]
    if(race2>0) .tmp <- .tmp[phdata(.tmp)$proband_race==race2,];gc()
    ## Check that cleft and environ are set correctly
    names(.data)[names(.data)==cleft];names(.data[[cleft]])[names(.data[[cleft]])==environ]
    if(F){rm(.data);gc()}
    result <- c()
    strata <- 16 ## Smoke
    if(environ=="vitamin") strata <- 19
    if(environ=="drink") strata <- 22
    strata
    write(strata,"strata.txt")
    pedInd <- paste(cleft,"samlet_modified.pedIndex",sep="_");pedInd
    write(pedInd,"pedIndex.txt")
    print(Sys.time())
    ## Ran on Haplin v. 6.2. 
    result <- haplinSlide(
      data=.tmp, pedIndex=scan("pedIndex.txt",what="character"),
      markers="ALL",
      winlength=1, 
      strata=as.integer(scan("strata.txt",what="integer")), poo=T, 
      response="mult", cpus=detectCores()-1,
      slaveOutfile = paste(cleft,"_",race,"_",environ,"_slave.txt",sep=""),
      reference="ref.cat", use.missing=T, verbose=F)
    print(Sys.time())
  save(result,file=paste("result_CL_CLP_",race,"_",environ,".RData",sep=""))
  } #environ
} #race

})

## Number of runs with errors
sum(is.na(result))

# head(result,10)

## Må kommentere ut resten når vi skal kjøre på Fimm
top.snps <- vector("list",3)
names(top.snps) <- c("smoke","drink","vitamin")
n.row <- 20
top.snps <- lapply(top.snps,function(x){
  x <- data.frame(matrix(NA,nrow=n.row,ncol=4))
  colnames(x) <- c("SNP","p","q","RRR")
  x}
)
top.snps.q <- top.snps.p <- top.snps;top.snps
## Setter miljø- og utfallsvariabel
race <- "all";cleft <- "CL_CLP";environ <- "smoke";cleft;environ
for(race in c("all","eur","asia")){
  tiff.name <- paste("CL_CLP_",race,".tiff",sep="");tiff.name
  tiff(file=tiff.name,units="cm",width=18,height=6,res=300)
  par(mar=c(3.1,2.9,.5,.1),cex=2,mfrow=c(1,3))
  for(environ in c("smoke","drink","vitamin")){
    load.data <- paste("result_CL_CLP_",race,"_",environ,".RData",sep="");load.data
    load(load.data)
    not.ok <- which(sapply(result,function(x){!is.data.frame(x)}))
    result2 <- result
    if(length(not.ok)>0) result2 <- result[-not.ok]
    pval <- sapply(result2,function(x){x[1,"gxe_poo_pval"]})
    names(pval) <- names(result2);head(pval);length(pval)
    if(race!="asia"|environ=="vitamin"){
      pQQ(pval,nlabs=0,mark=FALSE,cex.axis=1,cex.lab=1,lim=c(0,6.2))
      header <- "SMOKING";if(environ=="drink"){header <- "ALCOHOL USE"};if(environ=="vitamin"){header <- "VITAMIN USE"};header
      x.pos <- .5;y.pos=6;text(x.pos,y.pos,header,cex=.9,pos=4)
      if(environ=="drink"|race=="asia") mtext("Expected p-value (-log10 scale)",line=1.8,side=1,cex=.8)
      if(environ=="smoke"|race=="asia") mtext("Observed p-value (-log10 scale)",line=1.8,side=2,srt=90,cex=.8)
    ## Henter qvalue
    # source("http://bioconductor.org/biocLite.R")
    # biocLite("qvalue")
    require(qvalue)
    ## Henter q-verdier
    adj.p <- qvalue(pval)
    names(adj.p)
    summary(adj.p)
    head(sort(adj.p$qvalues),20)
    tail(sort(adj.p$qvalues))
    ## Skriv relevante SNP-er inn i hapmap.org
    ## Deretter google "entrez gene" og genecards.org for å finne genets betydning
    ## haplo.freq
    top.snps[[environ]][,"SNP"] <- names(head(sort(adj.p$pvalues),n.row))
    top.snps[[environ]][,"p"] <- head(sort(adj.p$pvalues),n.row)
    top.snps[[environ]][,"q"] <- head(sort(adj.p$qvalues),n.row);top.snps[[environ]]
    ## Calculate RRRs with 95% CIs
    ## Tabell
    tab <- data.frame(SNP=top.snps[[environ]][,"SNP"], 
                      p=signif(top.snps[[environ]][,"p"],2), 
                      q=signif(top.snps[[environ]][,"q"],2),
                      RRR=top.snps[[environ]][,"RRR"]);tab
    tab$SNP <- as.character(tab$SNP)
    tab.out <- paste("SNP","p","q","RRR",sep="\t");tab.out
    for(i in 1:n.row){
      tab.out <- paste(tab.out,paste(tab[i,],collapse="\t"),sep="\n")
    };tab.out
    tab.name <- paste("CL_CLP_",environ,"_",race,".doc",sep="");tab.name
    write(tab.out,file=tab.name)
    shell(paste("start",tab.name))
    }
  }
  dev.off()
  shell(paste("start",tiff.name))
}

## Tabell smoke
tab.tmp <- data.frame(matrix(NA,nrow=n.row,ncol=12));tab.tmp
tab.tmp[,seq(1,12,3)] <- top.snps[["smoke"]]
tab.tmp[,seq(2,12,3)] <- round(10^5*top.snps.p[["smoke"]],2);tab.tmp
tab.tmp[,seq(3,12,3)] <- round(top.snps.q[["smoke"]],2);tab.tmp
tab.smoke <- paste("CL\t\t\tCP\t\t\tCLP\t\t\tCL_CLP\t\t",sep="");tab.smoke
tab.smoke <- paste(tab.smoke,"SNP\tp(x10-5)\tq\tSNP\tp(x10-5)\tq\tSNP\tp(x10-5)\tq\tSNP\tp(x10-5)\tq",sep="\n");tab.smoke
i <- 1
for(i in 1:n.row){
  tab.smoke <- paste(tab.smoke,paste(tab.tmp[i,],collapse="\t"),sep="\n");tab.smoke
}
tab.smoke
write(tab.smoke,file="tab_smoke.doc")
shell("start tab_smoke.doc")

## Fra list til data.frame
result2 <- do.call(rbind.data.frame,result[1:30])
result2 <- result2[seq(2,nrow(result2),2),]
## Skriv til csv
write.csv2(result2,file="test.csv")

## RRcm.est er RR assosiert med allelet barnet arvet fra mor
## RRcf.est er RR assosiert med allelet barnet arvet fra far
## RRcm_RRcf.est er RRcm/RRcf
## RRdd er "dobbel-dose", dvs. RRcm*RRcf (under response="mult")


## Finne RR for top-SNP-ene
library(GenABEL)

marker <- vector("list",3);names(marker)=c("all","eur","asia");marker
marker <- lapply(marker,function(x){v=vector("list",3);names(v) <- c("smoke","drink","vitamin");v})
result <- cis <- marker;result;cis
## Get names of top SNPs from top.snps
marker[["all"]][["smoke"]] <- toupper(c("rs10097386","rs2383162","rs10738571","rs7419201","rs7541537","rs7042192","rs4977848","rs7920088","rs12740826","rs13173741","rs10757168","rs8181543","rs168283","rs17408603","rs11624380","rs2177971","rs7943401","rs3793861","rs4394682","rs7087489"))
marker[["all"]][["drink"]] <- toupper(c("rs7964474","rs999783","rs4982619","rs7945550","rs880813","rs2280025","rs11584506","rs10897066","rs2032442","rs163684","rs8025763","rs13008096","rs4699228","rs2723057","rs7201659","rs2151225","rs7197476","rs2367283","rs2914354","rs7209652"))
marker[["all"]][["vitamin"]] <- toupper(c("rs2302304","rs2689128","rs9572250","rs4875398","rs3909551","rs9371494","rs8101981","rs7939975","rs10495767","rs11673884","rs6489630","rs3815311","rs358017","rs7082286","rs921743","rs10764037","rs8112256","rs4569521","rs6830509","rs9503155"))
marker[["eur"]][["smoke"]] <- toupper(c("rs10763707","rs7541537","rs7419201","rs3793861","rs7087489","rs814518","rs4693142","rs4454616","rs2904096","rs2290682","rs6532013","rs1868368","rs2177971","rs6807522","rs17604550","rs12883776","rs7234787","rs8181543","rs4310561","rs3800036"))
marker[["eur"]][["drink"]] <- toupper(c("rs10496410","rs7579926","rs2294035","rs6975650","rs4876274","rs2245225","rs927318","rs10735337","rs6427247","rs12669493","rs13255561","rs12242535","rs943881","rs10491327","rs7945550","rs7232492","rs11242213","rs34352212","rs1990185","rs521419"))
marker[["eur"]][["vitamin"]] <- toupper(c("rs2689128","rs2237360","rs7793050","rs7766106","rs2809964","rs3859121","rs1092733","rs7559678","rs2366837","rs10084852","rs6446389","rs2242909","rs595536","rs6726527","rs12733019","rs8101981","rs17793145","rs4973310","rs3815311","rs8072885"))
marker[["asia"]][["smoke"]] <- NULL
marker[["asia"]][["drink"]] <- NULL
marker[["asia"]][["vitamin"]] <- toupper(c("rs1889976","rs259395","rs10798004","rs12431484","rs10518981","rs1940698","rs171477","rs9862866","rs865585","rs17591732","rs12630106","rs7316350","rs7336296","rs1499916","rs7153574","rs6439772","rs1348564","rs2360838","rs12204808","rs1407555"))
marker

## In common
## Smoke
marker$all$smoke[which(marker$all$smoke%in%marker$eur$smoke)];which(marker$all$smoke%in%marker$eur$smoke)
marker$eur$smoke[which(marker$eur$smoke%in%marker$all$smoke)];which(marker$eur$smoke%in%marker$all$smoke)
marker$all$smoke[which(marker$all$smoke%in%marker$asia$smoke)];which(marker$all$smoke%in%marker$asia$smoke)
marker$eur$smoke[which(marker$eur$smoke%in%marker$asia$smoke)];which(marker$eur$smoke%in%marker$asia$smoke)
## Drink
marker$all$drink[which(marker$all$drink%in%marker$eur$drink)];which(marker$all$drink%in%marker$eur$drink)
marker$eur$drink[which(marker$eur$drink%in%marker$all$drink)];which(marker$eur$drink%in%marker$all$drink)
marker$all$drink[which(marker$all$drink%in%marker$asia$drink)];which(marker$all$drink%in%marker$asia$drink)
marker$eur$drink[which(marker$eur$drink%in%marker$asia$drink)];which(marker$eur$drink%in%marker$asia$drink)
## Vitamin
marker$all$vitamin[which(marker$all$vitamin%in%marker$eur$vitamin)];which(marker$all$vitamin%in%marker$eur$vitamin)
marker$eur$vitamin[which(marker$eur$vitamin%in%marker$all$vitamin)];which(marker$eur$vitamin%in%marker$all$vitamin)
marker$all$vitamin[which(marker$all$vitamin%in%marker$asia$vitamin)];which(marker$all$vitamin%in%marker$asia$vitamin)
marker$eur$vitamin[which(marker$eur$vitamin%in%marker$asia$vitamin)];which(marker$eur$vitamin%in%marker$asia$vitamin)

## Functions from Miriam
.fconf.int <- function(.beta,.var,.alpha){
  c(exp(.beta-qnorm(1-.alpha/2)*sqrt(.var)),exp(.beta+qnorm(1-.alpha/2)*sqrt(.var)))
}
.fbeta <- function(.coef,x){
  .coef[["1"]][x,]-.coef[["0"]][x,]
} 
.fvar <- function(.cov,x){
  .cov[["1"]][x,x]+.cov[["0"]][x,x]
}
## Fra Miriam
rr.ci <- function(res,r=2){
  .cov <- attr(res,"cov")
  .coef <- attr(res,"coef")
  
  .beta <- .fbeta(.coef,r) #"c1" instead of 2?
  .var <- .fvar(.cov,r) #"c1" instead of 2?
  .chisq <- .beta^2/.var
  .RRR_c1 <- exp(.beta)
  .conf.int_c1 <- .fconf.int(.beta,.var,0.05)
  res <- c(.RRR_c1,.conf.int_c1,2*(1-pnorm(sqrt(.chisq))))
  names(res) <- c("RRR","lower","upper","p")
  return(res)
}

## RRRs for top SNPs
cleft <- "CL_CLP";environ <- "vitamin";race <- "eur";cleft;environ;race
for(environ in c("vitamin","smoke","drink")){
  for(race in c("all","eur","asia")){
    if(race=="asia"&environ%in%c("smoke","drink")) break
    if(race=="eur") race2 <- 1
    if(race=="asia") race2 <- 3
    strata=16 ## Smoke
    if(environ=="vitamin") strata=19
    if(environ=="drink") strata=22
    
    .tmp=.data[["CL_CLP"]][[environ]]
    if(race=="eur"|race=="asia") .tmp<-.tmp[phdata(.tmp)$proband_race==race2,]
    
    snps <- marker[[race]][[environ]]
    for(snp in snps){
      test <- haplinStrat(
        data=.tmp,pedIndex="CL_CLP_samlet_modified.pedIndex",
        markers=which(snpnames(.tmp)==toupper(snp)),
        strata=strata,poo=T,response="mult",cpus=1,
        reference="ref.cat",table.output = T,
        use.missing=T,verbose=F)
      
      ## Only works on Miriam's enhanced version of Haplin
      res <- gxe(test)
      # attributes(res)
      
      ## Fra Miriam
      .cov <- attr(res,"cov")
      .coef <- attr(res,"coef")
      
      .beta <- .fbeta(.coef,4)
      .var <- .fvar(.cov,4)
      .chisq <- .beta^2/.var
      .RRR_c1 <- exp(.beta)
      .conf.int_c1 <- .fconf.int(.beta,.var,0.05)
      
      ## Output
      .RRR_c1; .conf.int_c1 
      cis[[race]][[environ]] <- paste(cis[[race]][[environ]],paste(snp,"\t",round2(.RRR_c1,2)," (",round2(.conf.int_c1[1],2),"-",round2(.conf.int_c1[2],2),")",sep=""),sep="\n")
    }
  }
}

## Lager tabeller
for(race in c("all","eur","asia")){
  fil <- "SNP\tRR (95%CI)"
  for(environ in c("vitamin","smoke","drink")){
    fil <- paste(fil,"\n",environ,"\t",cis[[race]][[environ]],sep="")
  }
  file <- paste("RR_",race,".doc",sep="");file
  write(fil,file=file)
  shell(paste("start ",file,sep=""))
}

## Haplotyper
## SNPs in ANK3 (smoke, Europe)
## Order on chromosome is correct
ank <- c("rs3793861","rs7087489","rs4310561");ank
## Strata
## Smoke: 16; vitamin: 19; drink: 22
## Race
## 1-European, 2-African, 3-Asian, 4-Hispanic, 
## 5-Native Amer, 6-Malaysian, 7-Other
## Single SNP
.dataEur <- .data[["CL_CLP"]][["smoke"]][phdata(.data[["CL_CLP"]][["smoke"]])$proband_race==1]
str(.dataEur)
## Child effects, PoO effects and GxSmoke effects
## GxE
testGxE <- lapply(1:3,function(i){
  haplinStrat(
    data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
    markers=which(snpnames(.dataEur)%in%toupper(ank[i])),
    strata=16,poo=F,response="mult",cpus=1,
    reference="ref.cat",table.output = T,
    use.missing=T,verbose=F)
})
names(testGxE) <- ank
tabGxE <- lapply(testGxE,haptable)
## rs3793861 child: 1.06 (0.91,1.24) - MAF=0.30 - p=0.45
## rs7087489 child: 1.06 (0.91,1.24) - MAF=0.30 - p=0.44
## rs4310561 child: 1.12 (0.97,1.24) - MAF=0.34 - p=0.13
lapply(tabGxE,function(x){x[1:2,c("RR.est.","RR.lower","RR.upper","RR.p.value","haplofreq","haplos")]})
resGxE <- lapply(1:3,function(i){gxe(testGxE[[i]])})
names(resGxE) <- ank
## rs3793861 GxE: 1.20 (0.83,1.60) - p=0.41
## rs7087489 GxE: 1.10 (0.82,1.60) - p=0.42
## rs4310561 GxE: 1.20 (0.87,1.60) - p=0.27
## r=2 because there are only two alleles (major vs. minor)
lapply(1:3,function(i){signif(rr.ci(res=resGxE[[i]]),r=2)})

## PoO
testPoO <- lapply(1:3,function(i){
  haplinStrat(
    data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
    markers=which(snpnames(.dataEur)%in%toupper(ank[i])),
    strata=NULL,poo=T,response="mult",cpus=1,
    reference="ref.cat",table.output = T,
    use.missing=T,verbose=F)
})
names(testPoO) <- ank
## rs3793861 PoO: 1.39 (1.09,1.76) - p=0.007
## rs7087489 PoO: 1.40 (1.10,1.77) - p=0.006
## rs4310561 PoO: 1.31 (1.03,1.64) - p=0.02
tabPoO <- lapply(testPoO,haptable)
lapply(tabPoO,function(x){
  signif(x[,c("RRcm_RRcf.est.","RRcm_RRcf.lower","RRcm_RRcf.upper","RRcm_RRcf.p.value")],3)
})
## PoOxSmoke
testPoOxE <- lapply(1:3,function(i){
  haplinStrat(
    data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
    markers=which(snpnames(.dataEur)%in%toupper(ank[i])),
    strata=16,poo=T,response="mult",cpus=1,
    reference="ref.cat",table.output = T,
    use.missing=T,verbose=F)
})
names(testPoOxE) <- ank
tabPoOxE <- lapply(testPoOxE,haptable)
lapply(tabPoOxE,function(x){x[1:2,c("haplos","haplofreq")]})
resPoOxE <- lapply(1:3,function(i){gxe(testPoOxE[[i]])})
names(resPoOxE) <- ank
## rs3793861 PoOxE: 3.67 (2.13,6.32) - p=2.6e-6
## rs7087489 PoOxE: 3.63 (2.11,6.25) - p=3.1e-6
## rs4310561 PoOxE: 2.90 (1.75,4.83) - p=4.0e-5
## r=4 because that is the position of cm_cf1
lapply(1:3,function(i){signif(rr.ci(res=resPoOxE[[i]],r=4),3)})

## rs3793861-rs7087489
testGxE2a <- haplinStrat(
  data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
  markers=which(snpnames(.dataEur)%in%toupper(ank[1:2])),
  strata=16,poo=F,response="mult",cpus=1,
  reference="ref.cat",table.output = T,
  use.missing=T,verbose=F)
testGxE2a
tabGxE2a <- lapply(testGxE2a,haptable)
## rs7087489-rs4310561
testGxE2b <- haplinStrat(
  data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
  markers=which(snpnames(.dataEur)%in%toupper(ank[2:3])),
  strata=16,poo=F,response="mult",cpus=1,
  reference="ref.cat",table.output = T,
  use.missing=T,verbose=F)
testGxE2b
tabGxE2b <- lapply(testGxE2b,haptable)
## rs3793861-rs7087489-rs4310561
testGxE3 <- haplinStrat(
  data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
  markers=which(snpnames(.dataEur)%in%toupper(ank)),
  strata=16,poo=F,response="mult",cpus=1,
  reference="ref.cat",table.output = T,
  use.missing=T,verbose=F)
testGxE3
tabGxE3 <- lapply(testGxE3,haptable)
## rs3793861-rs7087489 child: 1.07 (0.92,1.24) - MAF=0.30 - p=0.39 - c-t/G-A
## rs7087489-rs4310561 child: 1.32 (0.95,1.83) - MAF=0.04 - p=0.10 - A-a/A-T
## rs7087489-rs4310561 child: 1.09 (0.94,1.28) - MAF=0.30 - p=0.26 - t-a/A-T
## rs3793861-rs7087489-rs4310561 child: 1.32 (0.95,1.83) - MAF=0.04 - p=0.10 - G-A-a/G-A-T
## rs3793861-rs7087489-rs4310561 child: 1.09 (0.94,1.28) - MAF=0.30 - p=0.26 - c-t-a/G-A-T
tabGxE2a$all[,c("marker","RR.est.","RR.lower","RR.upper","RR.p.value","haplofreq","haplos")]
tabGxE2b$all[,c("marker","RR.est.","RR.lower","RR.upper","RR.p.value","haplofreq","haplos")]
tabGxE3$all[,c("marker","RR.est.","RR.lower","RR.upper","RR.p.value","haplofreq","haplos")]
resGxE2a <- gxe(testGxE2a)
resGxE2b <- gxe(testGxE2b)
resGxE3 <- gxe(testGxE3)
## rs3793861-rs7087489 GxE: 1.10 (0.81-1.60) - p=0.45
## rs7087489-rs4310561 (A-a) GxE: 1.50 (0.71-3.10) - p=0.29
## rs7087489-rs4310561 (t-a) GxE: 1.20 (0.83-1.60) - p=0.37
## rs3793861-rs7087489-rs4310561 GxE (G-A-a): 1.50 (0.71-3.10) - p=0.29
## rs3793861-rs7087489-rs4310561 GxE (c-t-a): 1.20 (0.83-1.60) - p=0.37
## r=2 because there are only two haplotypes
signif(rr.ci(res=resGxE2a, r=2),2)
## r=3 and r=4 because there are three haplotypes
signif(rr.ci(res=resGxE2b, r=3),2)
signif(rr.ci(res=resGxE2b, r=4),2)
## r=3 and r=4 because there are three haplotypes
signif(rr.ci(res=resGxE3, r=3),2)
signif(rr.ci(res=resGxE3, r=4),2)

testPoO2a <- haplinStrat(
  data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
  markers=which(snpnames(.dataEur)%in%toupper(ank[1:2])),
  strata=NULL,poo=T,response="mult",cpus=1,
  reference="ref.cat",table.output = T,
  use.missing=T,verbose=F)
testPoO2b <- haplinStrat(
  data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
  markers=which(snpnames(.dataEur)%in%toupper(ank[2:3])),
  strata=NULL,poo=T,response="mult",cpus=1,
  reference="ref.cat",table.output = T,
  use.missing=T,verbose=F)
testPoO3 <- haplinStrat(
  data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
  markers=which(snpnames(.dataEur)%in%toupper(ank)),
  strata=NULL,poo=T,response="mult",cpus=1,
  reference="ref.cat",table.output = T,
  use.missing=T,verbose=F)
## rs3793861-rs7087489 PoO: 1.38 (1.08-1.75) - p=0.0083
## rs7087489-rs4310561 PoO (A-a): 0.93 (0.60-1.45) - p=0.74
## rs7087489-rs4310561 PoO (t-a): 1.34 (1.06-1.70) - p=0.02
## rs3793861-rs7087489-rs4310561 PoO (A-a): 0.93 (0.60,1.45) - p=0.75
## rs3793861-rs7087489-rs4310561 PoO (t-a): 1.35 (1.06-1.71) - p=0.01
tabPoO2a <- haptable(testPoO2a)
cbind(tabPoO2a[,"haplos"],signif(tabPoO2a[,c("RRcm_RRcf.est.","RRcm_RRcf.lower","RRcm_RRcf.upper","RRcm_RRcf.p.value")],3))
tabPoO2b <- haptable(testPoO2b)
cbind(tabPoO2b[,"haplos"],signif(tabPoO2b[,c("RRcm_RRcf.est.","RRcm_RRcf.lower","RRcm_RRcf.upper","RRcm_RRcf.p.value")],3))
tabPoO3 <- haptable(testPoO3)
cbind(tabPoO3[,"haplos"],signif(tabPoO3[,c("RRcm_RRcf.est.","RRcm_RRcf.lower","RRcm_RRcf.upper","RRcm_RRcf.p.value")],3))

## rs3793861-rs7087489
testPoOxE2a <- haplinStrat(
  data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
  markers=which(snpnames(.dataEur)%in%toupper(ank[1:2])),
  strata=16,poo=T,response="mult",cpus=1,
  reference="ref.cat",table.output = T,
  use.missing=T,verbose=F)
testPoOxE2a
tabPoOxE2a <- lapply(testPoOxE2a,haptable)
## rs7087489-rs4310561
testPoOxE2b <- haplinStrat(
  data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
  markers=which(snpnames(.dataEur)%in%toupper(ank[2:3])),
  strata=16,poo=T,response="mult",cpus=1,
  reference="ref.cat",table.output = T,
  use.missing=T,verbose=F)
testPoOxE2b
tabPoOxE2b <- lapply(testPoOxE2b,haptable)
## rs3793861-rs7087489-rs4310561
testPoOxE3 <- haplinStrat(
  data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
  markers=which(snpnames(.dataEur)%in%toupper(ank)),
  strata=16,poo=T,response="mult",cpus=1,
  reference="ref.cat",table.output = T,
  use.missing=T,verbose=F)
testPoOxE3
tabPoOxE3 <- lapply(testPoOxE3,haptable)
## rs3793861-rs7087489 child: 1.07 (0.92,1.24) - MAF=0.30 - p=0.39 - c-t/G-A
## rs7087489-rs4310561 child: 1.32 (0.95,1.83) - MAF=0.04 - p=0.10 - A-a/A-T
## rs3793861-rs7087489-rs4310561 child: 1.32 (0.95,1.83) - MAF=0.04 - p=0.10 - G-A-a/G-A-T
tabPoOxE2a$all[,c("RRcm_RRcf.est.","RRcm_RRcf.lower","RRcm_RRcf.upper","RRcm_RRcf.p.value","haplofreq","haplos")]
tabPoOxE2b$all[,c("RRcm_RRcf.est.","RRcm_RRcf.lower","RRcm_RRcf.upper","RRcm_RRcf.p.value","haplofreq","haplos")]
tabPoOxE3$all[,c("RRcm_RRcf.est.","RRcm_RRcf.lower","RRcm_RRcf.upper","RRcm_RRcf.p.value","haplofreq","haplos")]
resPoOxE2a <- gxe(testPoOxE2a)
resPoOxE2b <- gxe(testPoOxE2b)
resPoOxE3 <- gxe(testPoOxE3)
## rs3793861-rs7087489 PoOxE: 3.71 (2.16,6.39) - p=2.2e-6
## rs7087489-rs4310561 PoOxE (A-a): 1.57 (0.62-3.96) - p=0.34
## rs7087489-rs4310561 PoOxE (t-a): 3.65 (2.13-6.28) - p=2.7e-6
## rs3793861-rs7087489-rs4310561 PoOxE (G-A-a): 1.56 (0.62-3.94) - p=0.35
## rs3793861-rs7087489-rs4310561 PoOxE (c-t-a): 3.62 (2.10-6.21) - p=3.27e-6
## r=4 because that is the position of cm_cf2
signif(rr.ci(res=resPoOxE2a,r=4),3)
## r=7 and r=8 because those are the positions of cm_cf1 and cm_cf2
signif(rr.ci(res=resPoOxE2b,r=7),3)
signif(rr.ci(res=resPoOxE2b,r=8),3)
## r=7 and r=8 because those are the positions of cm_cf1 and cm_cf2
signif(rr.ci(res=resPoOxE3,r=7),3)
signif(rr.ci(res=resPoOxE3,r=8),3)

## SNPs in ARHGEF10 (drink, Europe)
arhgef <- c("rs2294035","rs4876274")
## Strata
## Smoke: 16; vitamin: 19; drink: 22
## Race
## 1-European, 2-African, 3-Asian, 4-Hispanic, 
## 5-Native Amer, 6-Malaysian, 7-Other
## Single SNP
.dataEur <- .data[["CL_CLP"]][["drink"]][phdata(.data[["CL_CLP"]][["drink"]])$proband_race==1]
str(.dataEur)
## Child effects, PoO effects and GxSmoke effects
testGxE <- lapply(1:2,function(i){
  haplinStrat(
    data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
    markers=which(snpnames(.dataEur)%in%toupper(arhgef[i])),
    strata=22,poo=F,response="mult",cpus=1,
    reference="ref.cat",table.output = T,
    use.missing=T,verbose=F)
})
names(testGxE) <- arhgef
tabGxE <- lapply(testGxE,haptable)
## rs2294035 child: 0.94 (0.82,1.08) - MAF=0.49 - p=0.38 - a/T
## rs4876274 child: 1.04 (0.90,1.20) - MAF=0.47 - p=0.57 - t/A
lapply(tabGxE,function(x){x[1:2,c("RR.est.","RR.lower","RR.upper","RR.p.value","haplofreq","haplos")]})
resGxE <- lapply(1:2,function(i){gxe(testGxE[[i]])})
names(resGxE) <- arhgef
## rs2294035 GxE: 1.20 (0.87,1.50) - p=0.32
## rs4876274 GxE: 0.90 (0.67,1.20) - p=0.47
## r=2 because there are only two SNPs (minor vs. major)
lapply(1:2,function(i){signif(rr.ci(res=resGxE[[i]], r=2),2)})

testPoO <- lapply(1:2,function(i){
  haplinStrat(
    data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
    markers=which(snpnames(.dataEur)%in%toupper(arhgef[i])),
    strata=NULL,poo=T,response="mult",cpus=1,
    reference="ref.cat",table.output = T,
    use.missing=T,verbose=F)
})
names(testPoO) <- arhgef
## rs2294035 PoO: 0.95 (0.75,1.20) - p=0.67
## rs4876274 PoO: 1.02 (0.80,1.29) - p=0.90
tabPoO <- lapply(testPoO,haptable)
lapply(tabPoO,function(x){
  signif(x[,c("RRcm_RRcf.est.","RRcm_RRcf.lower","RRcm_RRcf.upper","RRcm_RRcf.p.value")],3)
})
## PoOxDrink
testPoOxE <- lapply(1:2,function(i){
  haplinStrat(
    data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
    markers=which(snpnames(.dataEur)%in%toupper(arhgef[i])),
    strata=22,poo=T,response="mult",cpus=1,
    reference="ref.cat",table.output = T,
    use.missing=T,verbose=F)
})
names(testPoOxE) <- arhgef
tabPoOxE <- lapply(testPoOxE,haptable)
lapply(tabPoOxE,function(x){x[1:2,c("haplos","haplofreq")]})
resPoOxE <- lapply(1:2,function(i){gxe(testPoOxE[[i]])})
names(resPoOxE) <- arhgef
## rs2294035 PoOxE: 0.32 (0.19,0.51) - p=2.9e-6
## rs4876274 PoOxE: 2.99 (1.83,4.90) - p=1.3e-5
## r=4 because that is the position of cm_cf1 or cm_cf2
lapply(1:2,function(i){signif(rr.ci(res=resPoOxE[[i]],r=4),3)})

## rs2294035-rs4876274
testGxE2a <- haplinStrat(
  data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
  markers=which(snpnames(.dataEur)%in%toupper(arhgef)),
  strata=22,poo=F,response="mult",cpus=1,
  reference="ref.cat",table.output = T,
  use.missing=T,verbose=F)
testGxE2a
tabGxE2a <- lapply(testGxE2a,haptable)
## rs2294035-rs4876274 child: 1.15 (0.80,1.68) - MAF=0.04 - p=0.44 - T-A/a-A
## rs2294035-rs4876274 child: 1.04 (0.90,1.20) - MAF=0.47 - p=0.63 - T-t/a-A
tabGxE2a$all[,c("marker","RR.est.","RR.lower","RR.upper","RR.p.value","haplofreq","haplos")]
resGxE2a <- gxe(testGxE2a)
## rs2294035-rs4876274 GxE: 0.73 (0.33,1.60) - p=0.45
## rs2294035-rs4876274 GxE: 0.90 (0.67,1.20) - p=0.46
signif(rr.ci(res=resGxE2a,r=3),2)
signif(rr.ci(res=resGxE2a,r=4),2)

## PoO
testPoO2a <- haplinStrat(
  data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
  markers=which(snpnames(.dataEur)%in%toupper(arhgef)),
  strata=NULL,poo=T,response="mult",cpus=1,
  reference="ref.cat",table.output = T,
  use.missing=T,verbose=F)
## rs2294035-rs4876274 PoO: 1.41 (0.85,2.37) - p=0.19
## rs2294035-rs4876274 PoO: 1.00 (0.80,1.27) - p=0.98
tabPoO2a <- haptable(testPoO2a)
cbind(tabPoO2a[,"haplos"],signif(tabPoO2a[,c("RRcm_RRcf.est.","RRcm_RRcf.lower","RRcm_RRcf.upper","RRcm_RRcf.p.value")],3))

## rs2294035-rs4876274
testPoOxE2a <- haplinStrat(
  data=.dataEur,pedIndex="CL_CLP_samlet_modified.pedIndex",
  markers=which(snpnames(.dataEur)%in%toupper(arhgef)),
  strata=22,poo=T,response="mult",cpus=1,
  reference="ref.cat",table.output = T,
  use.missing=T,verbose=F)
testPoOxE2a
tabPoOxE2a <- lapply(testPoOxE2a,haptable)
tabPoOxE2a$all[,c("RRcm_RRcf.est.","RRcm_RRcf.lower","RRcm_RRcf.upper","RRcm_RRcf.p.value","haplofreq","haplos")]
resPoOxE2a <- gxe(testPoOxE2a)
## rs2294035-rs4876274 PoOxE: 1.56 (0.49,4.93) - p=0.45
## rs2294035-rs4876274 PoOxE: 3.20 (1.97,5.21) - p=2.8e-6
## r=7 and r=8 because that is where cm_cf2 and cm_cf3 are
signif(rr.ci(res=resPoOxE2a,r=7),3)
signif(rr.ci(res=resPoOxE2a,r=8),3)





