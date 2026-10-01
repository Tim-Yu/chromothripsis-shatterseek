#!/usr/bin/env Rscript
##################################################################################
## chromothripsis_original_runnable.r
##
## Stand-alone, runnable copy of the legacy GEL-sarcoma chromothripsis caller
## (legacy/chromothripsis_original.r, Dan/Atef).  The algorithm is UNCHANGED;
## only the plumbing was edited so that it runs end to end with Rscript:
##
##   * removed .libPaths()/install.packages() and the vcd / EMT libraries
##     (neither package is actually used: testExp() uses MASS::fitdistr +
##     stats::ks.test, testclass() uses stats::chisq.test);
##   * input/output taken from the command line:
##       Rscript chromothripsis_original_runnable.r --sv SVs.txt --cn CNAs.txt --out calls.tsv
##   * the duplicated mclapply(mc.cores=5) pass was removed (the original ran the
##     per-sample loop twice: once with mclapply, once with lapply) - plain lapply;
##   * no save(allCPs) side-effect; the table is written to --out;
##   * deriveTable() returns NULL gracefully when nothing is flagged;
##   * NA in the final rule (e.g. pRandomClass = NA when the region has no
##     DEL/DUP/INV) is reported as "Cluster" instead of NA, which is what Pui's
##     table shows (ifelse(NA, ...) in the original would give NA).
##
## See legacy/README_legacy.md for a step-by-step description of the algorithm.
##################################################################################
# Library paths for the GEL research environment, applied only when the
# directories exist, so nothing changes elsewhere. Equivalent to
#   .libPaths(c("/home/byu/R/x86_64-pc-linux-gnu-library/4.5", .libPaths(),
#               "/tools/aws-workspace-ubuntu-apps/ce/R/4.5.3", "/tools/aws-workspace-apps/ce/R/4.2.1/"))
# i.e. the personal library holding ShatterSeek first, the shared trees with the
# dependencies last. Override / extend with SHATTERSEEK_LIB (prepended) and
# SHATTERSEEK_EXTRA_LIBS="dir1:dir2" (appended).
local({
  personal <- c(Sys.getenv("SHATTERSEEK_LIB", ""), "/home/byu/R/x86_64-pc-linux-gnu-library/4.5")
  extra <- c("/tools/aws-workspace-ubuntu-apps/ce/R/4.5.3", "/tools/aws-workspace-apps/ce/R/4.2.1/",
             strsplit(Sys.getenv("SHATTERSEEK_EXTRA_LIBS", ""), ":")[[1]])
  personal <- personal[nzchar(personal) & dir.exists(personal)]
  extra <- extra[nzchar(extra) & dir.exists(extra)]
  if (length(personal) || length(extra)) .libPaths(c(personal, .libPaths(), extra))
})

suppressPackageStartupMessages({
    library(MASS)
    library(copynumber)
    library(GenomicRanges)
})

##################################################################################
## minimal CLI
##################################################################################
parseArgs <- function()
{
    args <- commandArgs(trailingOnly=TRUE)
    out <- list(sv=NULL, cn=NULL, out="chromothripsis_calls_original.tsv")
    i <- 1
    while(i <= length(args))
    {
        key <- sub("^--","",args[i])
        if(key=="help")
        {
            cat("Usage: Rscript chromothripsis_original_runnable.r --sv <SVs.txt> --cn <CNAs.txt> --out <calls.tsv>\n")
            quit(status=0)
        }
        if(i==length(args)) stop("Missing value for --",key)
        out[[key]] <- args[i+1]
        i <- i+2
    }
    if(is.null(out$sv) | is.null(out$cn)) stop("--sv and --cn are required (use --help)")
    out
}
opts <- parseArgs()

##################################################################################
N <- NULL
##################################################################################

##################################################################################
## wrapper around try
mytry <- function(x,retVal=NA,...)
{
    kk <- try(x,silent=T,...)
    if(inherits(kk,"try-error")) return(retVal)
    kk
}

## test if distribution of segment sizes is exponential (=random)
testExp <- function(chrCN,N=NULL)
{
    sizes <- diff(chrCN[,"end"])
    if(length(sizes)<4) stop("number of breakpoints too low")
    if(any(sizes<0))
    {
        print(chrCN)
        stop("non-ordered chromosomal positions")
    }
    fit <- fitdistr(sizes,"exponential")
    ks.test(sizes,pexp,rate=fit$estimate)$p.value
}

## performs PCF in interbreakpoint distance
## to identify regions of high density of breakpoints along chromosomes
doPCF <- function(nt, minEvents)
{
    diffs <- c(nt[1,"end"]-nt[1,"start"],diff(nt[,"end"]))
    res <- pcf(data.frame(chr=rep(nt$chr[1],length(diffs)),
                          positions=nt[,"end"],sample1=diffs),
               verbose=F,
               kmin=minEvents)
    res
}

## further annotates flagged regions for copy number states
getCNstates <- function(pcf,whichFlagged,chrCN,maxStates=3)
{
    starts <- c(0,cumsum(pcf$n.probes[-c(nrow(pcf))]))+1
    ends <- cumsum(pcf$n.probes)
    nStates <- list()
    for(i in whichFlagged)
    {
        wf <- as.logical(1:nrow(chrCN)*0)
        wf[starts[i]:ends[i]] <- TRUE
        states <- paste(chrCN[wf,"nMaj1"],chrCN[wf,"nMin1"],sep="-")
        sizes <- chrCN[wf,"end"]-chrCN[wf,"start"]
        isSubclonal <- chrCN[wf,"frac1"]
        ssizes <- tapply(sizes,as.factor(states),sum)/1000
        nsegs <- tapply(sizes,as.factor(states),length)
        fsizes <- ssizes/sum(ssizes)
        nS <- length(ssizes)
        if(nS<maxStates) maxStates <- nS
        nStates[[i]] <- list(states=states,
                             sizes=sizes,ssizes=ssizes,
                             isSubclonal=isSubclonal,
                             subclonalFrac=sum(sizes/100*as.numeric(isSubclonal!=1),na.rm=T)/sum(sizes/100,na.rm=T),
                             fsizes=fsizes,
                             nsegs=nsegs,
                             coveragemaxStates.modeSegs=sum(fsizes[order(nsegs,decreasing=T)][1:maxStates]),
                             coveragemaxStates.maxCov=sum(sort(fsizes,decreasing=T)[1:maxStates]))
    }
    nStates
}

## merges all density-flagged PCF segments of an arm into one region (per arm)
treatPCF <- function(pcf,
                     flagDensity=150000000/50,
                     minEvents=5)
{
    npcf <- pcf
    whichFlaggedP <- which(pcf$mean<flagDensity & pcf$n.probes>=minEvents & pcf$arm=="p")
    if(length(whichFlaggedP)>0)
    {
        whichFlaggedP <- min(whichFlaggedP):max(whichFlaggedP)
        if(length(whichFlaggedP)>1)
            npcf <- npcf[-c(whichFlaggedP[2:length(whichFlaggedP)]),]
        i <- whichFlaggedP[1]
        npcf[i,"start.pos"] <- pcf[min(whichFlaggedP),"start.pos"]
        npcf[i,"end.pos"] <- pcf[max(whichFlaggedP),"end.pos"]
        npcf[i,"n.probes"] <- sum(pcf[whichFlaggedP,"n.probes"])
        wfp <- whichFlaggedP
        npcf[i,"mean"] <- sum(pcf[wfp,"mean"]/1000*pcf[wfp,"n.probes"])/sum(pcf[wfp,"n.probes"])*1000
    }
    pcf <- npcf
    whichFlaggedP <- which(pcf$mean<flagDensity & pcf$n.probes>=minEvents & pcf$arm=="q")
    if(length(whichFlaggedP)>0)
    {
        whichFlaggedP <- min(whichFlaggedP):max(whichFlaggedP)
        if(length(whichFlaggedP)>1)
            npcf <- npcf[-c(whichFlaggedP[2:length(whichFlaggedP)]),]
        i <- whichFlaggedP[1]
        npcf[i,"start.pos"] <- pcf[min(whichFlaggedP),"start.pos"]
        npcf[i,"end.pos"] <- pcf[max(whichFlaggedP),"end.pos"]
        npcf[i,"n.probes"] <- sum(pcf[whichFlaggedP,"n.probes"])
        wfp <- whichFlaggedP
        npcf[i,"mean"] <- sum(pcf[wfp,"mean"]/1000*pcf[wfp,"n.probes"])/sum(pcf[wfp,"n.probes"])*1000
    }
    return(npcf)
}

## takes copynumber as input and flags region of high density of breakpoints
CP <- function(chrCN,
               flagDensity=150000000/50,
               minEvents=5)
{
    if(nrow(chrCN)>minEvents)
    {
        pcf <- doPCF(chrCN,minEvents=minEvents)
        pcf <- treatPCF(pcf,flagDensity,minEvents)
        starts <- cumsum(c(1,pcf$n.probes[-c(nrow(pcf))]))
        ends <- cumsum(pcf$n.probes)
        whichFlagged <- which(pcf$mean<flagDensity & pcf$n.probes>=minEvents)
        tE <- mytry(testExp(chrCN,N=NULL))
        tERegion <- lapply(whichFlagged,function(x)
        {
            res <- try(testExp(chrCN[starts[x]:ends[x],],N=N),silent=T)
            if(inherits(res,"try-error")) return(NULL)
            res
        })
        flaggedtE <- tE<0.05
        flaggedDensity <- length(whichFlagged)>0
        flagged <- flaggedtE | flaggedDensity
        copyNstates <- getCNstates(pcf,whichFlagged,chrCN)
        return(list(starts=chrCN[starts[whichFlagged],"start"],
                    ends=chrCN[ends[whichFlagged],"end"],
                    arms=pcf[whichFlagged,"arm"],
                    densities=pcf$mean[whichFlagged],
                    Nevents=pcf$n.probes[whichFlagged],
                    cnStates=copyNstates,
                    testExpChromosome=tE,
                    testExpRegion=tERegion,
                    flaggedtE=flaggedtE,
                    flaggedDensity=flaggedDensity,
                    flagged=flagged))
    }
    else
        return(list(flagged=FALSE,
                    minEventsReached=FALSE))
}

## check length of homozygous deletions
summaryHD <- function(cn)
{
    hd <- (cn$nMaj1==0 & cn$nMin1==0) | (cn$nMaj2==0 & cn$nMin2==0)
    sizes <- sum(cn[hd,"end"]/1000000-cn[hd,"start"]/1000000,na.rm=T)
}

## helper functions to make table from list of calls
tablise <- function(cp.chr,i,nms)
{
    ww <- which(!sapply(cp.chr$cnStates, is.null))[1:2]
    if(all(c("p","q")%in%cp.chr$arms))
        return(tablise2(cp.chr,i,nms,ww))
    c(as.character(nms),
      as.character(cp.chr$start),
      as.character(cp.chr$end),
      as.character(paste0("chr",i,cp.chr$arms)),
      as.character(cp.chr$densities),
      as.character(cp.chr$Nevents),
      as.character(cp.chr$testExpRegion[[1]]),
      as.character(paste(cp.chr$cnStates[[ww[1]]]$states,sep="",collapse=";")),
      as.character(paste(cp.chr$cnStates[[ww[1]]]$sizes,sep="",collapse=";")),
      as.character(cp.chr$cnStates[[ww[1]]]$coveragemaxStates.modeSegs),
      as.character(cp.chr$cnStates[[ww[1]]]$coveragemaxStates.maxCov))
}

tablise2 <- function(cp.chr,i,nms,ww)
{
    rbind(c(as.character(nms),
            as.character(cp.chr$start[1]),
            as.character(cp.chr$end[1]),
            as.character(paste0("chr",i,cp.chr$arms[1])),
            as.character(cp.chr$densities[1]),
            as.character(cp.chr$Nevents[1]),
            as.character(cp.chr$testExpRegion[[1]]),
            as.character(paste(cp.chr$cnStates[[ww[1]]]$states,sep="",collapse=";")),
            as.character(paste(cp.chr$cnStates[[ww[1]]]$sizes,sep="",collapse=";")),
            as.character(cp.chr$cnStates[[ww[1]]]$coveragemaxStates.modeSegs),
            as.character(cp.chr$cnStates[[ww[1]]]$coveragemaxStates.maxCov)),
          c(as.character(nms),
            as.character(cp.chr$start[2]),
            as.character(cp.chr$end[2]),
            as.character(paste0("chr",i,cp.chr$arms[2])),
            as.character(cp.chr$densities[2]),
            as.character(cp.chr$Nevents[2]),
            as.character(cp.chr$testExpRegion[[2]]),
            as.character(paste(cp.chr$cnStates[[ww[2]]]$states,sep="",collapse=";")),
            as.character(paste(cp.chr$cnStates[[ww[2]]]$sizes,sep="",collapse=";")),
            as.character(cp.chr$cnStates[[ww[2]]]$coveragemaxStates.modeSegs),
            as.character(cp.chr$cnStates[[ww[2]]]$coveragemaxStates.maxCov)))
}

## return min. fraction of region with max. 3 CN states
## (linear in the number of CN segments: 1.1 - 0.006*nb, capped at 1;
##  i.e. = 1 for nb <= 16, 0.8 at nb = 50, 0.5 at nb = 100)
returnCov <- function(nb_cna)
{
    b <- (.5-100/50*.8)/(1-100/50)
    a <- (.8-b)/50
    ## linear cap depending on number fragments; max at 1
    cov <- a*nb_cna+b
    cov[cov>1] <- 1
    cov
}

## derive table from list of calls
deriveTable <- function(allCPs)
{
    tt <- NULL
    for(i in 1:length(allCPs))
    {
        for(j in 1:length(allCPs[[i]]))
        {
            if(allCPs[[i]][[j]]$flagged & !is.na(allCPs[[i]][[j]]$flagged))
            {
                if(allCPs[[i]][[j]]$flaggedDensity)
                {
                    vv <- tablise(allCPs[[i]][[j]],names(allCPs[[i]])[j],names(allCPs)[i])
                    tt <- rbind(tt,vv)
                }
            }
        }
    }
    if(is.null(tt)) return(NULL)
    colnames(tt) <- c("samplename","start","end",
                      "chrArm",
                      "density",
                      "NbreakpointsCNA",
                      "p.KS_exp",
                      "CNA.states",
                      "CNA.sizes",
                      "coverage.mode",
                      "coverage.max")
    tt
}

## Chi^2 to test if random SV classes
testclass <- function(svs)
{
    classes <- svs[,"svclass"]
    types <- c("DEL",
               "h2hINV",
               "t2tINV",
               "DUP")
    observed <- as.vector(table(classes)[types])
    observed[is.na(observed)] <- 0
    Chi2 <- try(suppressWarnings(chisq.test(observed,p=rep(1/4,4))$p.value),silent=T)
    if(inherits(Chi2,"try-error")) Chi2 <- NA
    Chi2
}

## test random mate orders (not used)
testorder <- function(svs)
{
    COR <- try(suppressWarnings(cor.test(svs[,"pos1"],
                        svs[,"pos2"],
                        met="sp")$p.value),silent=T)
    if(inherits(COR,"try-error")) COR <- NA
    COR
}

## run two previous tests on flagged regions
tests <- function(ct)
{
    chr <- gsub("p","",gsub("q","",ct[4]))
    chr <- gsub("chr","",chr)
    if(chr=="23") chr <- "X"
    svs <- lSV[[which(samples==ct["samplename"])]]
    grSV1 <- GRanges(as.character(svs[,2]),IRanges(svs[,4],svs[,4]))
    grSV2 <- GRanges(as.character(svs[,5]),IRanges(svs[,7],svs[,7]))
    grCT <- GRanges(chr,
                    IRanges(as.numeric(ct[2]),as.numeric(ct[3])))
    ov1 <- findOverlaps(grSV1,grCT)
    ov2 <- findOverlaps(grSV2,grCT)
    keep <- unique(c(queryHits(ov1),queryHits(ov2)))
    svs <- svs[keep,]
    list(testclass(svs),
         testorder(svs[as.character(svs[,2])==chr & as.character(svs[,5])==chr,]),
         length(keep))
}
##################################################################################

##################################################################################
## Load SVs and CNAs (format prepared by Atef)
svs <- read.table(opts$sv,sep="\t",header=T,stringsAsFactors=FALSE)
cnas <- read.table(opts$cn,sep="\t",header=T,stringsAsFactors=FALSE)
##################################################################################

##################################################################################
## Unique samples
samples <- unique(as.character(cnas[,1]))
##################################################################################
## list of CNA and SV, columns renamed
lCN <- lapply(samples,function(x)
{
    tmp <- cnas[as.character(cnas[,1])==x,]
    colnames(tmp) <- c("samplename","chr","start","end","nMaj1","nMin1","frac1","nMaj2","nMin2","frac2","SD","ploidy")
    tmp[order(tmp[,2],tmp[,3],decreasing=F),]
})
lSV <- lapply(samples,function(x)
{
    tmp <- svs[as.character(svs[,1])==x,]
    colnames(tmp) <- c("samplename","chr1","strand1","pos1","chr2","strand2","pos2","svclass")
    tmp
})
##################################################################################

##################################################################################
## quick sanity check on size of 0+0 regions
hd <- sapply(lCN,summaryHD)
cat("Total size (Mb) of homozygous deletions per sample:\n"); print(summary(hd))
##################################################################################

##################################################################################
## list of flagged high density breakpoint regions
## also annotates p-value for test for distribution of segment sizes
chrs <- c(1:22,"X")
allCPs <- lapply(lCN,function(t)
{
    CPs <- lapply(chrs,function(x) {
        CP(t[as.character(t$chr)==x,])
    })
    names(CPs) <- chrs
    CPs
})
names(allCPs) <- samples
##################################################################################

##################################################################################
## make table from flagged regions
tt <- deriveTable(allCPs)
if(is.null(tt))
{
    cat("No density-flagged region in any sample; writing empty table to",opts$out,"\n")
    writeLines(paste(c("samplename","start","end","chrArm","density","NbreakpointsCNA","p.KS_exp",
                       "CNA.states","CNA.sizes","coverage.mode","coverage.max","pRandomClass",
                       "pRandomOrder","expectedCov","nbSVs","Calls"),collapse="\t"),opts$out)
    quit(status=0)
}
##################################################################################

##################################################################################
## test regions for random SV classes/types; and random mate position order (not used)
alltest <- apply(tt,1,function(x) tests(x))
pClass <- sapply(alltest,function(x) x[[1]])
pOrder <- sapply(alltest,function(x) x[[2]])
nbSVs <- sapply(alltest,function(x) x[[3]])
##################################################################################
## Augment table with p-values; expected min. region size covered by at most 3 copy number states;
## and number of breakpoints;
tt <- cbind(tt,pClass,pOrder)
colnames(tt)[(ncol(tt)-1):ncol(tt)] <- c("pRandomClass","pRandomOrder")
tt <- cbind(tt,returnCov(as.numeric(tt[,"NbreakpointsCNA"])))
colnames(tt)[ncol(tt)] <- c("expectedCov")
tt <- cbind(tt,nbSVs)
colnames(tt)[ncol(tt)] <- c("nbSVs")
##################################################################################

##################################################################################
## Goes through flagged regions and annotations and annotates as "Chromothripsis" or "Cluster"
isChromothripsis <- apply(tt,1,function(x)
{
    mode <- as.numeric(x["coverage.mode"])>as.numeric(x["expectedCov"])
    max <- as.numeric(x["coverage.max"])>as.numeric(x["expectedCov"])
    class <- as.numeric(x["pRandomClass"])>0.01
    expt <- as.numeric(x["p.KS_exp"])<0.05
    nbsv <- as.numeric(x["NbreakpointsCNA"])>15
    nbcna <- as.numeric(x["nbSVs"])>15
    (mode|max)&class&expt&(nbsv|nbcna)
})
tt <- cbind(tt,ifelse(isChromothripsis %in% TRUE,"Chromothripsis","Cluster"))
colnames(tt)[ncol(tt)] <- "Calls"
##################################################################################

##################################################################################
## Make data.frame
tt <- as.data.frame(tt,stringsAsFactors=FALSE)
tt$start <- as.numeric(as.character(tt$start))
tt$end <- as.numeric(as.character(tt$end))
tt$p.KS_exp <- as.numeric(as.character(tt$p.KS_exp))
tt$pRandomClass <- as.numeric(as.character(tt$pRandomClass))
tt$pRandomOrder <- as.numeric(as.character(tt$pRandomOrder))
rownames(tt) <- NULL
##################################################################################
## save table of annotated regions and calls
write.table(tt,
            file=opts$out,
            sep="\t",
            col.names=T,
            row.names=F,
            quote=F)
cat("Wrote",nrow(tt),"flagged regions to",opts$out,"\n")
##################################################################################
