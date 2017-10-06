
require(matrixcalc)

calculated_percentage_lower<-function(enhanced_dataframe,jadwin_name,lines_name,lines_dataframe){
  #selects specific sh2 distance value from dataframe corresponding to selected sh2 from histogram and barplot
  temp_lines_df<- as.numeric(lines_dataframe%>%filter(name == lines_name)%>%select(euclidean))
  lessthan<-nrow(enhanced_dataframe%>%
                   filter(pair2==jadwin_name)%>%
                   filter(euclidean < temp_lines_df))
  percentage<- (1-lessthan/(nrow(enhanced_dataframe)/2-1))*100 #display percentage below
  paste(jadwin_name, " has a smaller distance than ", percentage, "% of background", sep = "")
}

create_bar_graph<-function(sh2Data,sh2Name){
  ggplot(data=sh2Data, aes(x = reorder(pair1, euclidean), y = euclidean))+ 
    geom_bar(stat = "identity",
             col = "navyblue",
             #changes fill pattern for selected sh2 domain
             aes(fill = pair1 == sh2Name)) +
    scale_fill_manual(values = c('skyblue','white')) +
    theme(axis.text.x =element_blank(),legend.position ="none")+
    geom_text(size = 5,
              aes(label =pair1),
              position = position_stack(vjust = 0.5),
              angle = 90) +
    labs(x = "SH2 Domain") +
    labs(y = "Euclidean Distance") +
    if (sh2Name == "GRB2 ASH"){
      ggtitle("GRB2")
    }else{
      ggtitle(sh2Name)
    }
}

create_histogram<-function(pd, bin, lines_dataframe){
  ggplot(data = pd,aes(euclidean))+
  geom_histogram(binwidth=bin,
                 col="navyblue", 
                 fill="skyblue")+
  geom_vline(data = lines_dataframe,aes(xintercept=euclidean, colour=name),show.legend = FALSE,size = 1)+
  ggtitle(paste(str_split(as.character(pd$pair2[1])," ")[[1]][1],"Pairs")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Euclidean Distance") +
  labs(y = "Frequency")
}


get_distance<-function(A,B){
  subtracted<-A-B
  frobeniusNormSubtracted<-frobenius.norm(subtracted)
  pwmDistance<-frobeniusNormSubtracted/(dim(B)[2]+dim(A)[2])
}

#### The Following function has been taken and modified slightly from https://github.com/omarwagih/rmimp/blob/master/R/pwm-functions.r



# Constants
AA = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
AA_PRIORS_HUMAN  =   c(A=0.070, R=0.056, N=0.036, D=0.048,C=0.023,Q=0.047,E=0.071,G=0.066,H=0.026,I=0.044,
                       L=0.100,K=0.058,M=0.021,F=0.037,P=0.063,S=0.083,T=0.053,W=0.012,Y=0.027,V=0.060)
AA_PRIORS_YEAST  =   c(A=0.055, R=0.045, N=0.061, D=0.058, C=0.013, Q=0.039, E=0.064, G=0.05, H=0.022, I=0.066,
                       L=0.096, K=0.073, M=0.021, F=0.045, P=0.044, S=0.091, T=0.059, W=0.01, Y=0.034, V=0.056)

#' Construct position weight matrix
#' 
#' Makes a position weight matrix given aligned sequences.
#'
#' @param seqs Aligned sequences all of the same length
#' @param pseudocount Pseudocount factor. Final pseudocount is background probability * this factor
#' @param relative.freq Set to TRUE if each column should be divided by the sum
#' @param is.kinase.pwm Set to TRUE if matrix is being built for a kinase
#' @param priors Named character vector containing priors of amino acids.
#' @param do.pseudocounts TRUE if we are to add pseudocounts
#' @keywords internal pwm construct
#' @examples
#' # No examples
PWM <- function(seqs, pseudocount=0.01, relative.freq=T, is.kinase.pwm=T, priors=AA_PRIORS_HUMAN, do.pseudocounts=F){
  
  
  # Ensure same length characters 
  seq.len = sapply(seqs, nchar)
  num.pos = seq.len[1]
  if(! all(seq.len == num.pos)) stop('Unequal length of sequences')
  
  # List of valid amino acids, sorted
  namespace = AA
  
  # Match priors to aa 
  bg.prob = priors[match(namespace, names(priors))]
  
  # Make matrix of letters
  split = unlist( sapply(seqs, function(seq){strsplit(seq, '')}) )
  m = t( matrix(split, seq.len, length(split)/num.pos) )
  
  # Construct PWM
  pwm.matrix = apply(m, 2, function(pos.data){
    
    # Get frequencies 
    t = table(pos.data)
    
    # Match to aa
    ind = match(namespace, names(t))
    
    # Create column
    col = t[ind]
    col[is.na(col)] = 0
    names(col) = namespace
    
    # Do pseudocounts if were logging
    if(do.pseudocounts) col = col + (pseudocount * (20*bg.prob))
    
    # Do relative frequencies
    if(relative.freq) col = col / sum(col)
    
    #take log2 of frequency/background
    log2(col/bg.prob)
  })
  
  # Information content for MATCH score
  #ic2 = apply(pwm.matrix, 2, function(col) sum(col* log2(col/bg.prob), na.rm=T))
  
  #attr(pwm.matrix, 'pseudocount') = pseudocount
  #attr(pwm.matrix, 'match.ic') = ic2
  
  # Assign AA names to rows/pos col
  rownames(pwm.matrix) = namespace
  colnames(pwm.matrix) = 1:num.pos
  #attr(pwm.matrix, 'is.kinase.pwm') = is.kinase.pwm
  
  #attr(pwm.matrix, 'nseqs') = length(seqs)
  return(pwm.matrix)
}


