patched.readsangerseq<-function(filename)
  #This is a slightly modified version of readsangerseq() in the sangerseqR package
  #It fixes a problem with reading some .ab1 files, which appear to have an aberrant last 
  #character in the sequence strings, which sangerseq() cannot cope with. It returns a sangerseq object.
{ require(sangerseqR)
  fc <- file(filename, open = "rb")
  rawdata <- readBin(fc, what = "raw", n = 1.2 * file.info(filename)$size)
  close(fc)
  filetype <- suppressWarnings(rawToChar(rawdata[1:4]))
  if (filetype == ".scf") {
    seq <- read.scf(filename)
  }
  else if (filetype == "ABIF") {
    seq <- read.abif(filename)
    l<-nchar(seq@data$PBAS.1)
    if(! substr(seq@data$PBAS.1,l,l) %in% LETTERS) { #if last character is not uppercase text
      seq@data$PBAS.1<-substr(seq@data$PBAS.1,1,l-1)
    }
    l<-nchar(seq@data$PBAS.2)
    if(! substr(seq@data$PBAS.2,l,l) %in% LETTERS) { #if last character is not uppercase text
      seq@data$PBAS.2<-substr(seq@data$PBAS.2,1,l-1)
    }
  }
  else 
    {stop("Invalid File.")}
  return(sangerseq(seq))
}


align <- function(control, sample, guide, seqstart = 100, seqend = 700, maxshift = 10, rg1 = NA, rg2 = NA, er1 = F, er2 = F) {
  ##Arguments: 
  ##control_file (char) & experimental_file (char): Sanger chromatogram files (.ab1 or .scf) 
  ##   control is typically transfected with Cas9 but without guide seq
  ##   experimental has typically been treated Cas9 + guide RNA
  ## guide (char): 20 bp guide sequence that is recognized by the Cas9 
  ## seqstart (numeric): start of sequence read from where data will be included (because beginning of seq reads tends to be poor quality)
  ## seqend (numeric): last bp to be included in analysis (will be automatically adjusted if reads are shorter, see below)
  ## maxshift (numeric): range of basepair shifts (indels) to be analyzed, both positive and negative
  ## rg1, rg2 (numeric): [optional] the first (rg1) and last (rg2) base of the sequence region that is used for decomposition; will be automatically set if NA
  ##   Note: rg1&rg2 should be after the breaksite (if not, this will be corrected)
  
  require("Biostrings")
  require("sangerseqR")
  
  B<-c("A","C","G","T") #four bases, in the order that is always used by sangerseqR package
  
  #extract primary sequences as called by sequencer:
  sequence_ctr <- primarySeq(control)
  sequence_mut <- primarySeq(sample)
  
  #adjust seqend to shortest sequence if necessary:
  seqend<-min(seqend,length(sequence_ctr), length(sequence_mut))
  
  
  #Alignments:
  
  #find position of gRNA (also if on opposite strand) and calculate breaksite: 
  Dguide<-DNAString(guide)
  if(!length(Dguide)==20){
    er1 = T
    stop("guide sequence should be 20 nucleotides")}
  guide.align.f <- pairwiseAlignment(pattern=Dguide, subject=sequence_ctr, type="local")
  guide.hit.f <- identical(guide, as.character(subject(guide.align.f))) #forward hit if full match found
  guide.align.r <- pairwiseAlignment(pattern=reverseComplement(Dguide), subject=sequence_ctr, type="local")
  guide.hit.r <- identical(as.character(reverseComplement(Dguide)), as.character(subject(guide.align.r))) #reverse hit if full match found
   
  
  #there can only be one match with the top strands or with the bottom strand.
  if(guide.hit.f & !guide.hit.r) {breaksite<-start(subject(guide.align.f))+16}
  if(!guide.hit.f & guide.hit.r) {breaksite<-start(subject(guide.align.r))+3}
  if(guide.hit.f & guide.hit.r){ 
    er1 <- T
    stop("at least two gRNA matches")}
  if(!guide.hit.f & !guide.hit.r){
    er1 <- T
    stop(paste("no gRNA matches 
               \n investigate if there is a mismatch in control sequence due to sanger sequencing. If so, change guide into identical IUPAC nucleotides as the control sequence.
               \n guide forward:", Dguide, 
               "\n guide reverse complement:", reverseComplement(Dguide),
               "\n control sequence:", sequence_ctr))}
  
  #align the sample to the control sequence and calculate offset
  if(seqstart>(breaksite-maxshift)){
    er1 <- T
    stop(paste("the breaksite (",breaksite,") is too close to the start of the sequence read -> 
              If possible set start of sequence read lower.
              \n The sequence start of sequence read is maximal n bp breaksite - n bp of the chosen indel size range."))}
  seq_ctr <- substr(sequence_ctr, seqstart, breaksite-maxshift)
  seq_mut <- substr(sequence_mut, seqstart, breaksite-maxshift)
  align_seq <- pairwiseAlignment(pattern = seq_mut, subject = seq_ctr, type = "local")
  if(align_seq@score<20){
    er1 = T
    stop("there is no good alignment found between the control amd test sample -> 
        The alignment window is too small or is of bad quality or the control and test sample do not match")
  }
  offset_mut <- align_seq@pattern@range@start-align_seq@subject@range@start
  
  
  #extract control data:
  peak_ctr_loc <- peakPosMatrix(control)[,1]-1 #for some reason sangerseq() added 1, so we substract it again
  peak_ctr_height <- traceMatrix(control)[peak_ctr_loc,] #matrix with a column for each base 
  peak_ctr_height <- peak_ctr_height[1:seqend,]
  peak_ctr_height[is.na(peak_ctr_height)]<-0 #set NAs to 0
  colnames(peak_ctr_height)<-B
  
  #extract experimental data:
  peak_mut_loc <- peakPosMatrix(sample)[,1]-1 #for some reason sangerseq() added 1, so we substract it again
  peak_mut_height <- traceMatrix(sample)[peak_mut_loc,] #matrix with a column for each base 
  peak_mut_height <- peak_mut_height[1:seqend,]
  peak_mut_height[is.na(peak_mut_height)]<-0 #set NAs to 0
  colnames(peak_mut_height)<-B
  
  #set rg1 and rg2, if not provided by user:
  rg1<-ifelse(is.na(rg1), breaksite+maxshift+5, rg1)
  rg2<-ifelse(is.na(rg2), seqend-maxshift-5, rg2)
  
  #check if rg1 and rg2 are within meaningful range: 
  if(rg1< breaksite+maxshift+5) {
    rg1<- breaksite+maxshift+5
    warning(paste("left boundary of decomposed region corrected to", rg1))
  }
  
  if(rg2 > seqend-maxshift-5) {
    rg2<- seqend-maxshift-5
    warning(paste("right boundary of decomposed region corrected to",rg2))
  }
  
  if(rg2 > seqend-offset_mut) {
    rg2<- seqend-offset_mut-5
    warning(paste("right boundary of decomposed region corrected to",rg2))
  }
  
  if(rg2<rg1+maxshift*2) {
    er2 <- T
    stop(paste("boundaries of decomposed region are not acceptable -> 
        Set boundaries further apart or use smaller indel size if possible, 
        \n Maximum decomposition window spans from 5bp + n bp indel size range downstream of the break to 5bp + n bp indel size from the end of the shortest sequence read"))}
  
  return(list(
    ctr=peak_ctr_height, 
    mut=peak_mut_height, 
    seqstart=seqstart,
    seqend=seqend,
    maxshift=maxshift,
    rg1=rg1,
    rg2=rg2,
    breaksite=breaksite,
    offset_mut=offset_mut,
    B=B,
    er1=er1,
    er2=er2))
}

quality <- function(tide) {
  ## Needed for this function are the packages "Biostrings", "sangerseqR", "colorspace")
  ## All the arguments are generated in the function 'TIDE_import'. 
  ## ctr = peakheigths of the control sample (e.g. transfected CRISPR without guide seq)
  ## tide$mut = peakheigths of the sample that have had a DSB/repair
  ## tide$breaksite = site the crispr-guide is supposed to break according to literature (3 bp before the PAM sequence)
  ## tide$seqstart = start of sequence read from where data will be included (because beginning of seq reads tends to be poor quality)
  ## tide$seqend = last bp to be included in analysis (will be automatically adjusted if reads are shorter, see below)
  ## tide$maxshift = which basepair shifts (indels) you want to know the percentage of.
  # Note: tide$maxshift is the number to one direction, in the calculation it determines the shift to both direction (deletion & insertion) 
  ## tide$rg1/tide$rg2 = the sequence trace that is used for decomposition
  # Note: rg1&rg2 should be always after the breaksite
  # tide$offset_mut = the offset that seuquence trace of the sample has with repect to the control sequence trace.
  
  ## The function will return a plot of the percentages of aberrant sequence trance per location.
  # plot will indicate expected breaksite location
  # plot will indicate the sequence window that is used for decomposition   
  ## The function will return the difference percentages of aberrant sequences compared to the controL (~ amount DSB happend in sample)
  
  #Calculate percentage of each bp per peak and correct for the offset
  procent_ctr <- tide$ctr  
  procent_mut <- tide$mut;
  if (tide$offset_mut>0){
    procent_mut <- rbind(tide$mut[(1+tide$offset_mut):nrow(tide$mut),],matrix(NA,tide$offset_mut,4))
    procent_mut <- (procent_mut/(rowSums(procent_mut)))*100
    procent_ctr <- (tide$ctr/(rowSums(tide$ctr)))*100
  } else if (tide$offset_mut<0){
    procent_mut <- rbind(matrix(NA,-tide$offset_mut,4),tide$mut[1:(nrow(tide$mut)+tide$offset_mut),])
    procent_mut <- (procent_mut/(rowSums(procent_mut)))*100
    procent_ctr <- (tide$ctr/(rowSums(tide$ctr)))*100
  } else if (tide$offset_mut==0){
    procent_ctr <- (tide$ctr/(rowSums(tide$ctr)))*100
    procent_mut <- (tide$mut/(rowSums(tide$mut)))*100
  }
  
  ## calculate total percentage mutations
  percentage_mutation_ctr <- rowSums(procent_ctr * t(apply(procent_ctr,1,function(x){!(x==max(x))})))
  percentage_mutation_sample <- rowSums(procent_mut * t(apply(procent_ctr,1,function(x){!(x==max(x))})))
  
  #calculate average mutation percentage 
  meanper_ctr_prebreak <- mean(percentage_mutation_ctr[tide$seqstart:(tide$breaksite-20)])
  meanper_mut_prebreak <- mean(percentage_mutation_sample[tide$seqstart:(tide$breaksite-20)])
  meanper_ctr_postbreak <- mean(percentage_mutation_ctr[tide$breaksite:(tide$seqend-20)])
  meanper_mut_postbreak <- mean(percentage_mutation_sample[tide$breaksite:(tide$seqend-20)]) 
  
  #print the percentages of each shift
  percentage_mutation <- data.frame(percentage = round(c(meanper_ctr_prebreak, meanper_mut_prebreak, meanper_ctr_postbreak, meanper_mut_postbreak),1))

  rownames(percentage_mutation) <- c("mean % pre-break control sample", "mean % pre-break test sample", "mean % post-break control sample", "mean % post-break test sample");
  
  return(list(
    percentage_mutation_ctr=percentage_mutation_ctr, 
    percentage_mutation_sample=percentage_mutation_sample,
    percentage_mutation=percentage_mutation
  ))
  
}


decomposition <- function(tide,  p.threshold = 0.001) {   
  
  ## Needed for this function are the packages "nnls", "colorspace"
  ## All the arguments are generated in the function 'TIDE_import', except p.threshold
  ## tide$ctr = peakheigths of the control sample (e.g. transfected CRISPR without guide seq)
  ## tide$mut = peakheigths of the sample that have had a DSB/repair
  ## tide$maxshift = size range of indels to be considered in the decomposition.
  ##   note: calculation is always done in both directions, i.e. for both deletions and insertions of sizes 0:maxshift 
  ## tide$rg1 = first base in the sequence sequence traces that is used for decomposition
  ## tide$rg2 = last base in the sequence sequence traces that is used for decomposition
  # Note: rg1&rg2 should be always after the breaksite
  ## tide$offset_mut = the offset that seuquence trace of the sample has with repect to the control sequence trace.
  ## p.threshold = p-value signicance threshold
  
  ## The function will generate a barplot with the prediction of the most prominent indels in the population of cells  
  ## The function will return the percentages of each indel in the sample with associated p-value
  
  require("colorspace")
  require("nnls")
  
  shiftrange<-c(-tide$maxshift: tide$maxshift) 
  
  #decomposite tide$mut sequence data into indel combinations, 
  #separately for each base in c("A","C","G","T"). Stack up the data for the four bases in one aggragation matrix
  I_matrix <- c()
  I_vec <- c()
  
  for(b in tide$B) #loop through four bases
  {#simulate sequencing peak data for all hypothesized indels from control peaks:
    sim <- matrix(NA, nrow=tide$rg2-tide$rg1+1, ncol=tide$maxshift*2+1)
    colnames(sim)<-shiftrange
    for(i in shiftrange) {sim[,as.character(i)] <- tide$ctr[(tide$rg1:tide$rg2)-i,b]}
    I_matrix <- rbind(I_matrix, sim)
    I_vec <- c(I_vec,tide$mut[(tide$rg1:tide$rg2)+tide$offset_mut,b])
  }
  
  #non-negative linear fit:
  NNFIT <- nnls(I_matrix,I_vec)
  
  ## pvalue calculation (source: https://www.princeton.edu/~slynch/soc504/mult_reg2.pdf)
  # standard error:
  se <- sqrt(diag((sum((NNFIT$fitted-I_vec)^2)/(nrow(I_matrix)-(tide$maxshift*2+1)))*solve(t(I_matrix)%*%I_matrix)));
  
  # p value:
  pv <- 2*pnorm(-abs(NNFIT$x/se))
  
  #R^2
  Rsq <- cor(NNFIT$fit,I_vec)^2
  
  #components in percentages:
  comper<-(Rsq*100*(NNFIT$x/sum(NNFIT$x)))
   
  COL <- ifelse(pv<p.threshold,"red","black")
  COL[tide$maxshift+1] <- ifelse(pv[tide$maxshift+1]<p.threshold,"#FF000080","black")
      
  decomp.summary <- data.frame(percentage = round(comper,1), pvalue = as.character(signif(pv,2)))
  rownames(decomp.summary) <- shiftrange
  rownames(decomp.summary)[which(shiftrange>0)]=paste('+',rownames(decomp.summary)[which(shiftrange>0)],sep='')

  eff <- round((Rsq*100) - comper[tide$maxshift+1],1)
  
  invisible(list(
    pv=pv, 
    p.threshold=p.threshold, 
    NNFIT=NNFIT,
    comper=comper,
    Rsq=Rsq,
    COL=COL,
    shiftrange=shiftrange,
    summary=decomp.summary,
    efficiency=eff))  
}


insertion <- function(tide, indels, er3 = F) {   
  ## Needed for this function are the packages "Biostrings", "sangerseqR", "colorspace")
  ## All the arguments are generated in the function 'TIDE_import' and the function 'TIDE_decomposition'. 
  ## ctr = peakheigths of the control sample (e.g. transfected CRISPR without guide seq)
  ## tide$mut = peakheigths of the sample that have had a DSB/repair
  ## tide$breaksite = site the crispr-guide is supposed to break according to literature (3 bp before the PAM sequence)
  ## tide$seqstart = start of sequence read from where data will be included (because beginning of seq reads tends to be poor quality)
  ## tide$seqend = last bp to be included in analysis (will be automatically adjusted if reads are shorter, see below)
  ## tide$maxshift = which basepair shifts (indels) you want to know the percentage of.
  # Note: maxshift is the number to one direction, in the calculation it determines the shift to both direction (deletion & insertion) 
  ## tide$rg1/tide$rg2 = the sequence trace that is used for decomposition
  # Note: rg1&rg2 should be always after the breaksite
  ## tide$offset_mut = the offset that seuquence trace of the sample has with repect to the control sequence trace.
  ## indels$p.threshold = threshold of what you determine to be significant
  ## indels$pv = p-value of the various indels in the sample
  ## indels$NNFIT = non-negative linear fit
  
  ## The function will return a plot with the ratio of the inserted basepair when there is a +1 insertion due to the CRISPR-Cas9 system
  
  # get insertion bases
  # only continue with the significant shifts
  #which indels (e.g. -2,-1, 0, +1,... etc) are significant
  SigShift <- which(indels$pv<indels$p.threshold)-tide$maxshift-1; 
  #what are the corresponding percentages of the significant peaks
  PerShift <- indels$NNFIT$x[which(indels$pv<indels$p.threshold)];
  names(PerShift) <- SigShift;
  
  ## generate shift table
  #for each +1 insertion, one bp is added on top of the ctr sequence and 1 is removed from the bottom, for each +2 insertion 2 bps are added on top of the ctr sequence and 2 are removed from the bottom of the sequence, etc. 
  #This shifting makes a similution in the ctr sequence of the 'real' happened shifts
  
  SigShiftList <- rep(list(NULL),length(SigShift));
  names(SigShiftList) <- SigShift;
  
  for (i in 1:length(SigShift)){
    agg <- matrix(NA, abs(SigShift[i]),4); 
    
    if (sign(SigShift[i])==1){
      SigShiftList[[i]] <- rbind(agg,tide$ctr[1:(nrow(tide$ctr)-abs(SigShift[i])),]);
    } else if (sign(SigShift[i])== -1){
      SigShiftList[[i]] <- rbind(tide$ctr[(abs(SigShift[i])+1):(nrow(tide$ctr)),],agg);
    } else if (sign(SigShift[i])== 0){
      SigShiftList[[i]] <- tide$ctr;
    }
  } 
  
  #Get insertion
  #From this point, we separate the insertions from the deletions
  #Note: if you want to spot point mutations, you can take the wt along (change ">" to ">=")
  insertion_mtx <- matrix(NA,length(which(SigShift >0)),4);
  rownames(insertion_mtx) <- (SigShift)[which(SigShift >0)];
  colnames(insertion_mtx) <- tide$B;
  
  if(length(insertion_mtx)==0){
    er3 <- T
    stop("no insertions")}
  if(length(insertion_mtx)>0){
    for(i in 1:nrow(insertion_mtx)){
      subS <- SigShiftList[which(SigShift<as.numeric(rownames(insertion_mtx))[i])];
      subP <- PerShift[which(SigShift<as.numeric(rownames(insertion_mtx))[i])]
      
      #correct for the offset_mut in the sample
      lto <- tide$mut;
      if (tide$offset_mut>0){
        lto <- rbind(tide$mut[(1+tide$offset_mut):nrow(tide$mut),],matrix(NA,tide$offset_mut,4))
      } else if (tide$offset_mut<0){
        lto <- rbind(matrix(NA,-tide$offset_mut,4),tide$mut[1:(nrow(tide$mut)+tide$offset_mut),])
      }
      
      #substract the simulated mutated ctr sequence that contain only the deletions (SubS) from the real mutated sequence (lto) --> insertion mutation will remain
      for (j in 1:length(subS)){
        lto <- lto-subS[[names(subS)[j]]]*subP[names(subS)[j]]
      }
      
      insertion_mtx[i,] <- lto[tide$breaksite+as.numeric(rownames(insertion_mtx)[i]),]/PerShift[rownames(insertion_mtx)[i]];
    }
        
    #plot bp composition insertion
    one_insertion <- (insertion_mtx)*(insertion_mtx>0)
    insertion.summary <- round(((one_insertion/rowSums(one_insertion))*100),1)
    
    if(!as.numeric(rownames(insertion.summary))==1){
      er3 <- T
      stop("no +1 insertion")}
    
    return(list(
      insertion.summary=insertion.summary,
      er3=er3))
  }
}

