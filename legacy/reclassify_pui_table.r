#!/usr/bin/env Rscript
##################################################################################
## reclassify_pui_table.r
## Applies the new tiered rule of chromothripsis_improved.r to an EXISTING legacy
## output table (e.g. Re__Pui_Boyu/All_chromothripsis_calls.txt).
##
## What CAN be recomputed from the legacy columns:
##   - total CN per segment from CNA.states ("nMaj-nMin;...") -> osc2, osc3,
##     fracDirChange, nCNstates, maxCN
##   - NbreakpointsCNA, nbSVs (SVs with >=1 breakpoint in region)
##   - legacy rule components (coverage, p.KS_exp, pRandomClass)
## What CANNOT (needs raw SV coordinates):
##   - p.SVcluster (SV breakpoint clustering test)  -> treated as passed because
##     the region was density-flagged; column marked "NA(not computable)"
##   - nSV_intra / nSV_inter / partnerChr / linkedTo (inter-region linking that
##     can upgrade a LowConf call to HighConf) -> not applied
##   - fracCNbp_withSV, highLevelAmp (ploidy not in table; maxCN reported instead)
##
## Usage: Rscript legacy/reclassify_pui_table.r <legacy_calls.txt> <out.tsv> [--option value ...]
##################################################################################
a <- commandArgs(trailingOnly=TRUE)
inF <- if(length(a)>=1) a[1] else "Re__Pui_Boyu/All_chromothripsis_calls.txt"
outF <- if(length(a)>=2) a[2] else "legacy/test_output/pui_calls_reclassified.tsv"
thr <- list(min_cn_bp_high=10,min_sv_high=15,osc2_high=7,osc3_high=9,zigzag_high=0.7,zigzag_min_bp_high=15,
            min_cn_bp_low=7,min_sv_low=8,osc2_low=4,osc3_low=6,zigzag_low=0.6,zigzag_min_bp_low=8)
if(length(a)>2) for(i in seq(3,length(a),by=2)) thr[[sub("^--","",a[i])]] <- as.numeric(a[i+1])

## --- same oscillation functions as in chromothripsis_improved.r ---
oscillation2 <- function(cn){ n<-length(cn); if(n<3) return(0L); best<-0L; run<-0L
    for(i in 3:n){ if(cn[i]==cn[i-2] && cn[i]!=cn[i-1]){ run<-run+1L; best<-max(best,run+2L)} else run<-0L }; best }
oscillation3 <- function(cn){ n<-length(cn); if(n<3) return(0L); best<-0L
    for(i in 1:(n-1)){ j<-i; while(j<n && cn[j+1]!=cn[j] && length(unique(cn[i:(j+1)]))<=3) j<-j+1; if(j>i) best<-max(best,j-i+1L)}; best }
zigzagFraction <- function(cn){ d<-sign(diff(cn)); d<-d[d!=0]; if(length(d)<2) return(NA_real_); sum(d[-1]!=d[-length(d)])/(length(d)-1) }

tt <- read.delim(inF,stringsAsFactors=FALSE,check.names=FALSE)
cnTot <- lapply(strsplit(tt$CNA.states,";"),function(x) sapply(strsplit(x,"-"),function(y) sum(as.numeric(y))))
tt$CN_total <- sapply(cnTot,paste,collapse=";")
tt$nCNstates <- sapply(cnTot,function(x) length(unique(x)))
tt$maxCN <- sapply(cnTot,max)
tt$osc2 <- sapply(cnTot,oscillation2)
tt$osc3 <- sapply(cnTot,oscillation3)
tt$fracDirChange <- round(sapply(cnTot,zigzagFraction),3)
tt$p.SVcluster <- "NA(not computable: needs SV coordinates)"
tt$linkedTo <- "NA(not computable: needs SV coordinates)"

nb <- as.numeric(tt$NbreakpointsCNA); nsv <- as.numeric(tt$nbSVs)
oscStrong <- with(tt,osc2>=thr$osc2_high | osc3>=thr$osc3_high | (!is.na(fracDirChange) & fracDirChange>=thr$zigzag_high & nb>=thr$zigzag_min_bp_high))
oscWeak   <- with(tt,osc2>=thr$osc2_low  | osc3>=thr$osc3_low  | (!is.na(fracDirChange) & fracDirChange>=thr$zigzag_low  & nb>=thr$zigzag_min_bp_low))
clustStrong <- nsv>=thr$min_sv_high   # p.SVcluster assumed passed (region is density-flagged)
clustWeak   <- nsv>=thr$min_sv_low
high <- nb>=thr$min_cn_bp_high & clustStrong & oscStrong
low  <- !high & nb>=thr$min_cn_bp_low & clustWeak & oscWeak
tt$Calls_original <- tt$Calls
tt$Calls_new <- ifelse(high,"Chromothripsis_HighConf",ifelse(low,"Chromothripsis_LowConf","Cluster"))
tt$evidence_new <- paste0("CNbp=",nb,";SV=",nsv,";osc2=",tt$osc2,";osc3=",tt$osc3,";zigzag=",tt$fracDirChange,
                          ";osc=",ifelse(oscStrong,"strong",ifelse(oscWeak,"weak","none")),
                          ";clust=",ifelse(clustStrong,"strong",ifelse(clustWeak,"weak","none")),"(pSVcluster not computable)")
## why the legacy rule failed
tt$legacy_fail_reason <- apply(tt,1,function(x){
    r <- c()
    if(!(as.numeric(x["coverage.mode"])>as.numeric(x["expectedCov"]) | as.numeric(x["coverage.max"])>as.numeric(x["expectedCov"]))) r <- c(r,"coverage<=expectedCov")
    if(is.na(x["pRandomClass"]) || !(as.numeric(x["pRandomClass"])>0.01)) r <- c(r,"pRandomClass NA/<=0.01")
    if(!(as.numeric(x["p.KS_exp"])<0.05)) r <- c(r,"p.KS_exp>=0.05")
    if(!(as.numeric(x["NbreakpointsCNA"])>15 | as.numeric(x["nbSVs"])>15)) r <- c(r,"<=15 CN bps and <=15 SVs")
    if(length(r)==0) "passed" else paste(r,collapse="; ")})
tt$Calls <- NULL
write.table(tt,outF,sep="\t",quote=FALSE,row.names=FALSE)
cat("\n| GMS | region | CNbp | SVs | CN_total | osc2 | osc3 | zigzag | legacy call | legacy failed on | new call |\n|---|---|---|---|---|---|---|---|---|---|---|\n")
for(i in seq_len(nrow(tt))) cat(sprintf("| %s | %s:%.1f-%.1fMb | %s | %s | %s | %d | %d | %s | %s | %s | **%s** |\n",tt$Sample[i],tt$chrArm[i],tt$start[i]/1e6,tt$end[i]/1e6,
    tt$NbreakpointsCNA[i],tt$nbSVs[i],tt$CN_total[i],tt$osc2[i],tt$osc3[i],tt$fracDirChange[i],tt$Calls_original[i],tt$legacy_fail_reason[i],tt$Calls_new[i]))
cat("\nWrote",outF,"\n")
