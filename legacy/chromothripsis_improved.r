#!/usr/bin/env Rscript
##################################################################################
## chromothripsis_improved.r
##
## Improved version of the legacy GEL-sarcoma chromothripsis caller
## (legacy/chromothripsis_original.r).  Same inputs, same output columns (plus
## new ones), same overall structure:
##
##   1. per sample / chromosome: PCF on inter-breakpoint distances flags regions
##      of dense breakpoints (as in the original), now ALSO from SV breakpoints
##      and with a flank-insensitive "span" density;
##   2. each flagged region is annotated (legacy statistics + new ones: CN
##      oscillation statistics, SV breakpoint counts by type, SV clustering
##      test, CN<->SV breakpoint concordance, inter-region links);
##   3. a tiered final rule: Chromothripsis_HighConf / Chromothripsis_LowConf /
##      Cluster.  The legacy rule is still evaluated and reported in
##      Calls_legacyRule for comparison.
##
## Usage:
##   Rscript chromothripsis_improved.r --sv SVs.txt --cn CNAs.txt --out calls.tsv [--option value ...]
##   Rscript chromothripsis_improved.r --help      (lists all tunable thresholds)
##
## Input formats (identical to the legacy script):
##   SVs.txt : samplename chr1 strand1 pos1 chr2 strand2 pos2 svclass   (svclass in DEL,DUP,h2hINV,t2tINV,TRA)
##   CNAs.txt: samplename chr start end nMaj1 nMin1 frac1 nMaj2 nMin2 frac2 SD ploidy
##
## See legacy/README_legacy.md for the rationale behind every change.
##################################################################################
# Extra library trees (GEL research environment): appended when they exist, so
# shared dependencies are found without editing the script. Override / extend
# with the environment variable SHATTERSEEK_EXTRA_LIBS="dir1:dir2".
local({
  extra <- c("/tools/aws-workspace-ubuntu-apps/ce/R/4.5.3", "/tools/aws-workspace-apps/ce/R/4.2.1/",
             strsplit(Sys.getenv("SHATTERSEEK_EXTRA_LIBS", ""), ":")[[1]])
  extra <- extra[nzchar(extra) & dir.exists(extra)]
  if (length(extra)) .libPaths(c(.libPaths(), extra))
})

suppressPackageStartupMessages({
    library(MASS)          # fitdistr (KS-exponential test, kept as informative)
    library(copynumber)    # pcf
    library(GenomicRanges) # overlaps
})

##################################################################################
## 0. Command line / tunable thresholds
##################################################################################
defaults <- list(
    sv = NULL, cn = NULL, out = "chromothripsis_calls_improved.tsv",
    ## --- region flagging (legacy: flagDensity = 150e6/50 = 3 Mb, minEvents = 5) ---
    flagDensity        = 3e6,   # a PCF segment is flagged when its mean OR span inter-breakpoint distance (bp) is below this
    minEvents          = 5,     # kmin of PCF and min. number of CN segments ("probes") in a flagged segment
    flag_from_sv       = TRUE,  # also flag dense regions from SV breakpoints (catches SV clusters with few CN changes)
    minEventsSV        = 8,     # min. number of SV breakpoints in an SV-flagged PCF segment
    trim_gap           = 6e6,   # SV-derived region hulls are trimmed of terminal breakpoints separated by more than this (bp)
    assembly           = "hg38",# genome build used by copynumber::pcf to assign chromosome arms (legacy used hg19)
    ## --- annotation ---
    bp_tol             = 20000, # CN breakpoint is "SV-supported" when an SV breakpoint lies within this distance (bp)
    highAmp_above_ploidy = 5,   # maxCN >= ploidy + this  => high-level amplification flag
    ## --- final rule: HighConf ---
    min_cn_bp_high     = 10,    # min. number of CN segments in region (NbreakpointsCNA)
    min_sv_high        = 15,    # min. number of SVs with >=1 breakpoint in region (nbSVs)
    p_cluster_high     = 0.01,  # SV-breakpoint clustering test (binomial, chromosome-wide uniform) p-value
    osc2_high          = 7,     # >= this many adjacent CN segments alternating between 2 total-CN states
    osc3_high          = 9,     # >= this many adjacent CN segments alternating among <=3 total-CN states
    zigzag_high        = 0.7,   # fraction of CN direction changes (up/down reversals) ...
    zigzag_min_bp_high = 15,    # ... needed together with at least this many CN segments
    ## --- final rule: LowConf ---
    min_cn_bp_low      = 7,
    min_sv_low         = 8,
    p_cluster_low      = 0.05,
    osc2_low           = 4,
    osc3_low           = 6,
    zigzag_low         = 0.6,
    zigzag_min_bp_low  = 8,
    ## --- inter-region linking (multi-chromosome events such as DFSP rings) ---
    min_inter_link     = 5,     # inter-chromosomal SVs joining two flagged regions of the same sample => "linked"
    ## --- legacy criteria, optional (default: informative only) ---
    require_ks         = FALSE, # require p.KS_exp < 0.05 (legacy)
    require_random_class = FALSE, # require pRandomClass > 0.01 (legacy)
    cores              = 1
)

parseArgs <- function(defaults)
{
    args <- commandArgs(trailingOnly=TRUE)
    opts <- defaults
    i <- 1
    while(i <= length(args))
    {
        key <- sub("^--","",args[i])
        if(key=="help")
        {
            cat("Usage: Rscript chromothripsis_improved.r --sv <SVs.txt> --cn <CNAs.txt> --out <calls.tsv> [--option value]\n\nOptions and defaults:\n")
            for(k in setdiff(names(defaults),c("sv","cn","out"))) cat(sprintf("  --%-22s %s\n",k,as.character(defaults[[k]])))
            quit(status=0)
        }
        if(!key %in% names(defaults)) stop("Unknown option --",key," (see --help)")
        if(i==length(args)) stop("Missing value for --",key)
        val <- args[i+1]
        if(is.logical(defaults[[key]])) val <- as.logical(val)
        else if(is.numeric(defaults[[key]])) val <- as.numeric(val)
        opts[[key]] <- val
        i <- i+2
    }
    if(is.null(opts$sv) | is.null(opts$cn)) stop("--sv and --cn are required (use --help)")
    opts
}
opts <- parseArgs(defaults)

##################################################################################
## 1. Helper functions carried over from the legacy script (unchanged logic)
##################################################################################
mytry <- function(x,retVal=NA,...)
{
    kk <- try(x,silent=TRUE,...)
    if(inherits(kk,"try-error")) return(retVal)
    kk
}

## KS test of segment sizes against a fitted exponential (legacy "random = exponential")
testExp <- function(chrCN)
{
    sizes <- diff(chrCN[,"end"])
    if(length(sizes)<4) stop("number of breakpoints too low")
    if(any(sizes<0)) stop("non-ordered chromosomal positions")
    fit <- fitdistr(sizes,"exponential")
    suppressWarnings(ks.test(sizes,pexp,rate=fit$estimate)$p.value)
}

## PCF on inter-breakpoint distances (legacy).  `positions` are breakpoint
## coordinates (CN segment ends, or SV breakpoints); the first "distance" is the
## distance from the chromosome start (legacy used end-start of the 1st segment).
doPCF <- function(chr,positions,firstDiff,minEvents,assembly)
{
    diffs <- c(firstDiff,diff(positions))
    pcf(data.frame(chr=rep(chr,length(diffs)),positions=positions,sample1=diffs),
        verbose=FALSE,kmin=minEvents,assembly=assembly)
}

## legacy coverage-by-<=3-states statistics
cnCoverage <- function(segs,maxStates=3)
{
    states <- paste(segs[,"nMaj1"],segs[,"nMin1"],sep="-")
    sizes <- segs[,"end"]-segs[,"start"]
    ssizes <- tapply(sizes,as.factor(states),sum)/1000
    nsegs <- tapply(sizes,as.factor(states),length)
    fsizes <- ssizes/sum(ssizes)
    maxStates <- min(maxStates,length(ssizes))
    list(states=states,sizes=sizes,
         coverage.mode=sum(fsizes[order(nsegs,decreasing=TRUE)][1:maxStates]),
         coverage.max=sum(sort(fsizes,decreasing=TRUE)[1:maxStates]))
}

## legacy expected coverage: 1.1 - 0.006*nb, capped at 1 (=1 for nb<=16)
returnCov <- function(nb_cna)
{
    b <- (.5-100/50*.8)/(1-100/50)
    a <- (.8-b)/50
    cov <- a*nb_cna+b
    cov[cov>1] <- 1
    cov
}

## legacy chi-square on equiprobable SV classes (TRAs excluded)
testclass <- function(svs)
{
    types <- c("DEL","h2hINV","t2tINV","DUP")
    observed <- as.vector(table(svs[,"svclass"])[types])
    observed[is.na(observed)] <- 0
    Chi2 <- try(suppressWarnings(chisq.test(observed,p=rep(1/4,4))$p.value),silent=TRUE)
    if(inherits(Chi2,"try-error")) Chi2 <- NA
    Chi2
}

## legacy Spearman test on mate positions (never used in the rule)
testorder <- function(svs)
{
    COR <- try(suppressWarnings(cor.test(svs[,"pos1"],svs[,"pos2"],method="spearman")$p.value),silent=TRUE)
    if(inherits(COR,"try-error")) COR <- NA
    COR
}

##################################################################################
## 2. New helper functions
##################################################################################
## Flag PCF segments: legacy criterion (mean inter-breakpoint distance) OR the
## flank-insensitive span density (segment span / number of internal gaps).
## The legacy mean includes the distance from the chromosome/arm start to the
## first breakpoint, which alone pushed e.g. 17 breakpoints over 47 Mb (2.6 Mb
## spacing) above the 3 Mb threshold when preceded by a 15 Mb flank.
flagPCF <- function(pcf,flagDensity,minEvents)
{
    span <- (pcf$end.pos-pcf$start.pos)/pmax(pcf$n.probes-1,1)
    (pcf$mean<flagDensity | span<flagDensity) & pcf$n.probes>=minEvents
}

## Merge all flagged PCF segments of one arm into one (legacy treatPCF logic),
## returning one row per arm with probe index range.
mergeFlaggedByArm <- function(pcf,flagged)
{
    starts <- cumsum(c(1,pcf$n.probes[-nrow(pcf)]))
    ends <- cumsum(pcf$n.probes)
    out <- NULL
    for(arm in c("p","q"))
    {
        w <- which(flagged & pcf$arm==arm)
        if(length(w)==0) next
        w <- min(w):max(w)
        n <- sum(pcf$n.probes[w])
        out <- rbind(out,data.frame(arm=arm,
                                    probeStart=starts[min(w)],probeEnd=ends[max(w)],
                                    n.probes=n,
                                    mean=sum(pcf$mean[w]*pcf$n.probes[w])/n,
                                    stringsAsFactors=FALSE))
    }
    out
}

## Longest run of adjacent segments alternating between exactly 2 total-CN
## states (ShatterSeek-like "max number of oscillating CN segments, 2 states").
oscillation2 <- function(cn)
{
    n <- length(cn)
    if(n<3) return(0L)
    best <- 0L; run <- 0L
    for(i in 3:n)
    {
        if(!is.na(cn[i]) && !is.na(cn[i-1]) && !is.na(cn[i-2]) && cn[i]==cn[i-2] && cn[i]!=cn[i-1])
        {
            run <- run+1L
            best <- max(best,run+2L)
        }
        else run <- 0L
    }
    best
}

## Longest run of adjacent segments with adjacent-differing CN using at most
## 3 distinct total-CN states ("3 states" oscillation).
oscillation3 <- function(cn)
{
    n <- length(cn)
    if(n<3) return(0L)
    best <- 0L
    for(i in 1:(n-1))
    {
        j <- i
        while(j<n && !is.na(cn[j+1]) && !is.na(cn[j]) && cn[j+1]!=cn[j] && length(unique(cn[i:(j+1)]))<=3) j <- j+1
        if(j>i) best <- max(best,j-i+1L)
    }
    best
}

## Fraction of direction reversals (up->down / down->up) along the CN profile.
## A random walk gives ~2/3; a BFB-like staircase or a simple step gives ~0;
## chromothripsis/ring amplicons give values close to 1.
zigzagFraction <- function(cn)
{
    d <- sign(diff(cn))
    d <- d[d!=0 & !is.na(d)]
    if(length(d)<2) return(NA_real_)
    sum(d[-1]!=d[-length(d)])/(length(d)-1)
}

## Binomial test: SV breakpoints inside region vs. expected under a uniform
## distribution of the chromosome's (or genome's) SV breakpoints.
clusterTest <- function(k,n,p)
{
    if(is.na(n) || n==0 || is.na(p) || p<=0 || p>=1) return(NA_real_)
    binom.test(k,n,p,alternative="greater")$p.value
}

##################################################################################
## 3. Read data
##################################################################################
svs <- read.table(opts$sv,sep="\t",header=TRUE,stringsAsFactors=FALSE)
cnas <- read.table(opts$cn,sep="\t",header=TRUE,stringsAsFactors=FALSE)
colnames(cnas)[1:12] <- c("samplename","chr","start","end","nMaj1","nMin1","frac1","nMaj2","nMin2","frac2","SD","ploidy")
colnames(svs)[1:8] <- c("samplename","chr1","strand1","pos1","chr2","strand2","pos2","svclass")
cnas$chr <- sub("^chr","",as.character(cnas$chr))
svs$chr1 <- sub("^chr","",as.character(svs$chr1))
svs$chr2 <- sub("^chr","",as.character(svs$chr2))
samples <- unique(as.character(cnas$samplename))
chrs <- c(1:22,"X")

##################################################################################
## 4. Per sample: flag regions, annotate
##################################################################################
processSample <- function(s)
{
    cn <- cnas[cnas$samplename==s,]
    cn <- cn[order(match(cn$chr,chrs),cn$start),]
    sv <- svs[svs$samplename==s,]
    chrLen <- tapply(cn$end,cn$chr,max)
    genomeLen <- sum(chrLen)
    nSVbp_genome <- 2*nrow(sv)
    ploidy <- suppressWarnings(as.numeric(cn$ploidy[1]))

    regions <- NULL
    for(ch in chrs)
    {
        chrCN <- cn[cn$chr==ch,]
        if(nrow(chrCN)==0) next
        cand <- NULL
        ## --- 4a. CN-breakpoint based flagging (legacy) ---
        if(nrow(chrCN)>opts$minEvents)
        {
            pcfCN <- mytry(doPCF(ch,chrCN$end,chrCN$end[1]-chrCN$start[1],opts$minEvents,opts$assembly),NULL)
            if(!is.null(pcfCN))
            {
                fl <- flagPCF(pcfCN,opts$flagDensity,opts$minEvents)
                m <- mergeFlaggedByArm(pcfCN,fl)
                if(!is.null(m)) for(k in seq_len(nrow(m)))
                {
                    idx <- m$probeStart[k]:m$probeEnd[k]
                    ## hull of CN breakpoints (segment boundaries) inside the PCF segment,
                    ## i.e. without the telomeric/centromeric flank of the first/last segment
                    bnd <- unique(c(chrCN$start[idx][idx>1],chrCN$end[idx][idx<nrow(chrCN)]))
                    if(length(bnd)<2) bnd <- c(chrCN$start[min(idx)],chrCN$end[max(idx)])
                    cand <- rbind(cand,data.frame(arm=m$arm[k],start=min(bnd),end=max(bnd),
                                                  source="CN",density=m$mean[k],stringsAsFactors=FALSE))
                }
            }
        }
        ## --- 4b. SV-breakpoint based flagging (new) ---
        if(isTRUE(opts$flag_from_sv))
        {
            P <- sort(unique(c(sv$pos1[sv$chr1==ch],sv$pos2[sv$chr2==ch])))
            if(length(P)>opts$minEventsSV)
            {
                pcfSV <- mytry(doPCF(ch,P,P[1],opts$minEventsSV,opts$assembly),NULL)
                if(!is.null(pcfSV))
                {
                    fl <- flagPCF(pcfSV,opts$flagDensity,opts$minEventsSV)
                    m <- mergeFlaggedByArm(pcfSV,fl)
                    if(!is.null(m)) for(k in seq_len(nrow(m)))
                    {
                        ## hull of the SV breakpoints, trimmed of isolated terminal breakpoints
                        ## (a PCF segment must contain >= kmin probes and can swallow a distant SV)
                        Q <- P[m$probeStart[k]:m$probeEnd[k]]
                        while(length(Q)>2 && Q[2]-Q[1]>opts$trim_gap) Q <- Q[-1]
                        while(length(Q)>2 && Q[length(Q)]-Q[length(Q)-1]>opts$trim_gap) Q <- Q[-length(Q)]
                        cand <- rbind(cand,data.frame(arm=m$arm[k],start=min(Q),end=max(Q),
                                                      source="SV",density=m$mean[k],stringsAsFactors=FALSE))
                    }
                }
            }
        }
        if(is.null(cand)) next
        ## --- 4c. merge overlapping CN- and SV-derived regions on the same arm ---
        cand <- cand[order(cand$start),]
        merged <- NULL
        for(k in seq_len(nrow(cand)))
        {
            if(!is.null(merged) && cand$start[k]<=merged$end[nrow(merged)] && cand$arm[k]==merged$arm[nrow(merged)])
            {
                j <- nrow(merged)
                merged$end[j] <- max(merged$end[j],cand$end[k])
                merged$source[j] <- paste(sort(unique(c(strsplit(merged$source[j],"\\+")[[1]],cand$source[k]))),collapse="+")
                if(cand$source[k]=="CN") merged$density[j] <- cand$density[k]
                else if(is.na(merged$density[j])) merged$density[j] <- cand$density[k]
            }
            else merged <- rbind(merged,cand[k,])
        }
        merged$chr <- ch
        regions <- rbind(regions,merged)
    }
    if(is.null(regions)) return(NULL)

    ## --- 4d. annotate every region ---
    ann <- lapply(seq_len(nrow(regions)),function(r)
    {
        ch <- regions$chr[r]; st <- regions$start[r]; en <- regions$end[r]
        chrCN <- cn[cn$chr==ch,]
        segs <- chrCN[chrCN$end>=st & chrCN$start<=en,]
        cnTot <- segs$nMaj1+segs$nMin1
        cov <- cnCoverage(segs)
        nSeg <- nrow(segs)
        ## SVs overlapping the region
        in1 <- sv$chr1==ch & sv$pos1>=st & sv$pos1<=en
        in2 <- sv$chr2==ch & sv$pos2>=st & sv$pos2<=en
        hit <- in1|in2
        svR <- sv[hit,]
        inter <- svR[svR$chr1!=svR$chr2,]
        partner <- if(nrow(inter)>0) sort(table(ifelse(inter$chr1==ch,inter$chr2,inter$chr1)),decreasing=TRUE) else NULL
        nSVbp_chr <- sum(sv$chr1==ch)+sum(sv$chr2==ch)
        nSVbp_region <- sum(in1)+sum(in2)
        regLen <- en-st+1
        ## CN breakpoints supported by an SV breakpoint within bp_tol
        bnd <- segs$start[-1]
        svbp <- c(sv$pos1[sv$chr1==ch],sv$pos2[sv$chr2==ch])
        fracSupp <- if(length(bnd)>0 && length(svbp)>0) mean(sapply(bnd,function(b) any(abs(svbp-b)<=opts$bp_tol))) else NA
        osc2 <- oscillation2(cnTot); osc3 <- oscillation3(cnTot); zz <- zigzagFraction(cnTot)
        maxCN <- max(cnTot,na.rm=TRUE)
        data.frame(
            samplename=s, start=st, end=en, chrArm=paste0("chr",ch,regions$arm[r]),
            density=regions$density[r],
            NbreakpointsCNA=nSeg,
            p.KS_exp=mytry(testExp(segs)),
            CNA.states=paste(cov$states,collapse=";"),
            CNA.sizes=paste(cov$sizes,collapse=";"),
            coverage.mode=cov$coverage.mode, coverage.max=cov$coverage.max,
            pRandomClass=testclass(svR),
            pRandomOrder=testorder(svR[svR$chr1==ch & svR$chr2==ch,]),
            expectedCov=returnCov(nSeg),
            nbSVs=nrow(svR),
            ## ---- new columns ----
            flagSource=regions$source[r],
            region_size_Mb=round(regLen/1e6,2),
            span_density=round(regLen/max(nSeg-1,1)),
            nCNbp=max(nSeg-1,0),
            CN_total=paste(cnTot,collapse=";"),
            nCNstates=length(unique(cnTot)),
            minCN=min(cnTot,na.rm=TRUE), maxCN=maxCN, ploidy=ploidy,
            highLevelAmp=!is.na(ploidy) & maxCN>=ploidy+opts$highAmp_above_ploidy,
            osc2=osc2, osc3=osc3, fracDirChange=round(zz,3),
            nSVbp_region=nSVbp_region,
            nSV_intra=sum(in1&in2),
            nSV_partialChr=sum(xor(in1,in2) & sv$chr1==sv$chr2),
            nSV_inter=nrow(inter),
            partnerChr=if(is.null(partner)) "" else paste(paste0(names(partner),":",partner)[1:min(3,length(partner))],collapse=";"),
            frac_inter_topPartner=if(is.null(partner)) NA else round(partner[1]/sum(partner),3),
            nSVbp_chr=nSVbp_chr,
            p.SVcluster=clusterTest(nSVbp_region,nSVbp_chr,regLen/chrLen[[ch]]),
            p.SVcluster_genome=clusterTest(nSVbp_region,nSVbp_genome,regLen/genomeLen),
            fracCNbp_withSV=round(fracSupp,3),
            stringsAsFactors=FALSE)
    })
    ann <- do.call(rbind,ann)

    ## --- 4e. inter-region links (SVs joining two flagged regions on different chromosomes) ---
    ann$linkedTo <- ""; ann$nSV_linked <- 0L
    if(nrow(ann)>1)
    {
        inReg <- function(chr,pos,r) chr==ann$chr[r] & pos>=ann$start[r] & pos<=ann$end[r]
        ann$chr <- sub("[pq]$","",sub("^chr","",ann$chrArm))
        for(a in seq_len(nrow(ann)))
        {
            links <- c()
            for(b in seq_len(nrow(ann)))
            {
                if(a==b || ann$chr[a]==ann$chr[b]) next
                nl <- sum((inReg(sv$chr1,sv$pos1,a) & inReg(sv$chr2,sv$pos2,b)) |
                          (inReg(sv$chr2,sv$pos2,a) & inReg(sv$chr1,sv$pos1,b)))
                if(nl>=opts$min_inter_link) links[ann$chrArm[b]] <- nl
            }
            if(length(links)>0)
            {
                ann$linkedTo[a] <- paste(paste0(names(links),":",links),collapse=";")
                ann$nSV_linked[a] <- sum(links)
            }
        }
        ann$chr <- NULL
    }
    ann
}

workers <- max(1,as.integer(opts$cores))
tt <- if(workers>1) parallel::mclapply(samples,processSample,mc.cores=workers) else lapply(samples,processSample)
tt <- do.call(rbind,tt)

legacyCols <- c("samplename","start","end","chrArm","density","NbreakpointsCNA","p.KS_exp","CNA.states","CNA.sizes",
                "coverage.mode","coverage.max","pRandomClass","pRandomOrder","expectedCov","nbSVs","Calls")
if(is.null(tt) || nrow(tt)==0)
{
    cat("No flagged region in any sample; writing empty table to",opts$out,"\n")
    writeLines(paste(legacyCols,collapse="\t"),opts$out)
    quit(status=0)
}

##################################################################################
## 5. Final classification
##################################################################################
## legacy rule, for comparison
legacyRule <- with(tt,(coverage.mode>expectedCov | coverage.max>expectedCov) & pRandomClass>0.01 & p.KS_exp<0.05 &
                       (NbreakpointsCNA>15 | nbSVs>15))
tt$Calls_legacyRule <- ifelse(legacyRule %in% TRUE,"Chromothripsis","Cluster")

## new tiered rule
oscStrong <- with(tt,(osc2>=opts$osc2_high | osc3>=opts$osc3_high |
                      (!is.na(fracDirChange) & fracDirChange>=opts$zigzag_high & NbreakpointsCNA>=opts$zigzag_min_bp_high)))
oscWeak   <- with(tt,(osc2>=opts$osc2_low | osc3>=opts$osc3_low |
                      (!is.na(fracDirChange) & fracDirChange>=opts$zigzag_low & NbreakpointsCNA>=opts$zigzag_min_bp_low)))
pOK <- function(p,thr) is.na(p) | p<=thr
clustStrong <- with(tt,nbSVs>=opts$min_sv_high & pOK(p.SVcluster,opts$p_cluster_high))
clustWeak   <- with(tt,nbSVs>=opts$min_sv_low & pOK(p.SVcluster,opts$p_cluster_low))
legacyGate <- rep(TRUE,nrow(tt))
if(isTRUE(opts$require_ks)) legacyGate <- legacyGate & (tt$p.KS_exp<0.05) %in% TRUE
if(isTRUE(opts$require_random_class)) legacyGate <- legacyGate & (tt$pRandomClass>0.01) %in% TRUE

## a region is "linked" when >= min_inter_link SVs join it to another flagged region
## of the same sample that itself carries (at least weak) oscillation + clustering evidence
candidateKey <- paste(tt$samplename,tt$chrArm)
supportive <- candidateKey[oscWeak & clustWeak]
linked <- sapply(seq_len(nrow(tt)),function(i)
{
    if(tt$linkedTo[i]=="") return(FALSE)
    partners <- sub(":.*","",strsplit(tt$linkedTo[i],";")[[1]])
    any(paste(tt$samplename[i],partners) %in% supportive)
})

high <- tt$NbreakpointsCNA>=opts$min_cn_bp_high & clustStrong & (oscStrong | (oscWeak & linked)) & legacyGate
low  <- !high & tt$NbreakpointsCNA>=opts$min_cn_bp_low & clustWeak & oscWeak & legacyGate
tt$Calls <- ifelse(high %in% TRUE,"Chromothripsis_HighConf",ifelse(low %in% TRUE,"Chromothripsis_LowConf","Cluster"))

## CN pattern label (descriptive)
tt$CN_pattern <- ifelse(tt$osc2>=opts$osc2_high,"oscillating_2states",
                 ifelse(tt$osc3>=opts$osc3_low | tt$osc2>=opts$osc2_low,"oscillating_2-3states",
                 ifelse(!is.na(tt$fracDirChange) & tt$fracDirChange>=opts$zigzag_low & tt$nCNstates>3,
                        ifelse(tt$highLevelAmp,"multistate_highlevel_amplicon","multistate_zigzag"),"other")))

## human-readable evidence summary
tt$evidence <- paste0("CNbp=",tt$NbreakpointsCNA,
                      ";SV=",tt$nbSVs,"(intra=",tt$nSV_intra,",inter=",tt$nSV_inter,")",
                      ";osc2=",tt$osc2,";osc3=",tt$osc3,";zigzag=",tt$fracDirChange,
                      ";pClust=",signif(tt$p.SVcluster,2),
                      ";osc=",ifelse(oscStrong,"strong",ifelse(oscWeak,"weak","none")),
                      ";clust=",ifelse(clustStrong,"strong",ifelse(clustWeak,"weak","none")),
                      ifelse(linked,";linked","" ))

## column order: legacy columns first, then the new ones
newCols <- setdiff(colnames(tt),legacyCols)
tt <- tt[,c(legacyCols,newCols)]
rownames(tt) <- NULL
write.table(tt,file=opts$out,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
cat("Wrote",nrow(tt),"flagged regions to",opts$out,"\n")
print(table(tt$Calls))
