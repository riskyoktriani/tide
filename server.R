library(shiny)
library(knitr)
library(sangerseqR)
source('polypeakfunctions.R')
source('TIDE_functions.R')

examplefile1 <- "example1.ab1"
examplefile2 <- "example2.ab1"
exampleguide <- "CATGCCGAGAGTGATCCCGG"

shinyServer(function(input, output, session) {

  inputdata <- reactive({ 
    
  if(input$example) {
    
    control <- patched.readsangerseq(examplefile1);
    sample <- patched.readsangerseq(examplefile2);
    guide <- cleanstring(exampleguide);
      
    return(list(
      control=control,
      sample=sample,
      guide=guide));
  } 
  
  else if(!(is.null(input$control) | is.null(input$sample))) {
    
    control <- patched.readsangerseq(input$control$datapath);
    sample <- patched.readsangerseq(input$sample$datapath);
    guide <- cleanstring(as.character(input$guide));
    
    return(list(
      control=control,
      sample=sample,
      guide=guide));
  }
  })
  
#run align function
a <- reactive({
  if((!is.null(input$control) & !is.null(input$sample) & !(input$guide =="")) | (input$example)){
    return(align(inputdata()$control,
                 inputdata()$sample,
                 inputdata()$guide,
                 seqstart=input$range_align,
                 maxshift=input$maxshift,
                 rg1=input$range_decom[1],
                 rg2=input$range_decom[2]))
    } 
 })

#conditions
#file uploaded yes/no?
output$fileUploaded <- reactive({
  return((!is.null(input$control) & !is.null(input$sample) & !(input$guide =="")) | (input$example))
  })
outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)

#error messages
output$error1 <- reactive({
  return(a()$er1==F)
})
outputOptions(output, "error1", suspendWhenHidden = FALSE)

output$error2 <- reactive({
  return(a()$er2==F)
})
outputOptions(output, "error2", suspendWhenHidden = FALSE)

output$error3 <- reactive({
  return(i()$er3==F)
})
outputOptions(output, "error3", suspendWhenHidden = FALSE)


#run quality function
q <- reactive({
  if((!is.null(input$control) & !is.null(input$sample) & !(input$guide =="")) | (input$example)){
    quality(a())}
    })

#output error messages
output$text7 <- renderUI({ 
  
    texta <- paste("Default settings:")
    textb <- paste("alignment window =", input$range_align, "-", a()$breaksite-10)
    textc <- paste("sequence range end =", a()$seqend)
    textd <- paste("decomposition window =",  rg1=a()$rg1, "-", a()$rg2)
    texte <- paste("indel size =", input$maxshift)
    textf <- paste("p threshold =", d()$p.threshold)
    HTML(paste(texta, textb, textd, texte, textf, sep = '<br/>'))
    
    })

#output error messages
output$text1 <- renderText({ 
  if(input$range_align>(a()$breaksite-a()$maxshift)) {
    paste("Warning: the breaksite (",a()$breaksite,") is too close to the start of the sequence read -> If possible set start of sequence read lower.
          The sequence start of sequence read is maximal n bp breaksite - n bp of the chosen indel size range.")}  
}) 
output$text2 <- renderText({ 
  if(!is.na(input$range_decom[1]) & input$range_decom[1]<(a()$breaksite+a()$maxshift+5)) {
  paste("Warning: left boundary of decomposition window was adjusted", a()$rg1,".", 
        "It must be at least 5bp plus the maximum indel size downstream of the expected break site")}
  }) 
output$text3 <- renderText({ 
  if(!is.na(input$range_decom[2]) & input$range_decom[2]>(a()$seqend-a()$maxshift-5)) {
  paste("Warning: right boundary of decomposition window was adjusted to", a()$rg2,".", 
        "It cannot be more than the length of the shortest sequence read minus the maximum indel size minus 5bp.")}
  })

#plot aberrant sequence signal
output$aberrant_sequence <- renderPlot({
  plot(q()$percentage_mutation_sample, 
       type="h", col="green3", 
       xlim=c(input$range_align, a()$seqend), 
       ylim=c(0,100),
       xlab="basepair",
       ylab="% of aberrant sequences")
  lines(q()$percentage_mutation_ctr, type="h", col="black")
  
  legend("topleft",legend=c("control sample", "test sample"), bty="n", pch=15, col=c(1,3))
  
  #show decomposition window
  if(input$parameter) {rect(input$range_decom[1], 110 , input$range_decom[2], 110, density = 1, xpd=TRUE, col="grey", lwd=6) 
                        text(input$range_decom[1]+((input$range_decom[2]-input$range_decom[1])/2), (110+6), xpd=TRUE, col="grey", as.character("region for decomposition")) }    
  
  if(!input$parameter) {rect(a()$rg1, 110 , a()$rg2, 110, density = 1, xpd=TRUE, col="grey", lwd=6) 
                       text(a()$rg1+((a()$rg2-a()$rg1)/2), (110+6), xpd=TRUE, col="grey", as.character("region for decomposition")) }    
 
  
  #indicate theoretical breaksite
  if (a()$breaksite>0){
    abline(v=a()$breaksite, lty=5,lwd=3,col='blue')
    legend("topright",legend=paste('expected cut at ',a()$breaksite,'bp',sep=''),text.col='blue', bty="n")}
  else{
    legend("topright",legend='no cut',text.col='blue', bty="n")
    }
})

#run decomposition function
d <- reactive({
  if((!is.null(input$control) & !is.null(input$sample) & !(input$guide =="")) | (input$example)){
    decomposition(a(), p.threshold =input$p.threshold)}
  })

output$text4 <- renderText({ 
   if(input$range_align>(a()$breaksite-a()$maxshift)) {}
  }) 
  
#plot decomposition graph              
output$decom <- renderPlot({
  #barplot:
  bp <- barplot(d()$comper, 
                col=d()$COL, 
                border = d()$COL, 
                names.arg=d()$shiftrange, 
                ylim=c(0, max(d()$comper+10)), 
                xlab="<--deletion     insertion-->", 
                ylab="% of sequences",  
                xaxt='n')
  
  
  #make x-axis 
  a <- min(ceiling(d()$shiftrange/5)*5)
  p <- pretty (c(a:-a), n=(round((length(d()$shiftrange)-1)/5,0)-1))
  axis(1,at=bp[p+max(p)+1+a+max(d()$shiftrange)],labels=p)
  
  #above each group of bars: show percentage (mean across four bases)
  text(bp[d()$pv<input$p.threshold], (d()$comper+5)[d()$pv<input$p.threshold], as.character(((round(d()$comper,1))[d()$pv<input$p.threshold])))
  
  #display Rsq values as an indication of the accuracy:
  legend("topright",legend=as.expression(c(bquote(p < .(input$p.threshold)), bquote(p >= .(input$p.threshold)))), title= as.expression(bquote(R^2 == .(round(d()$Rsq,2)))), pch=15, col=c("red",'black'), bty="n")
  })

#guide efficiency
output$text5 <- renderText({ 
paste("overall efficiency =", d()$eff, "%")
})

#table decomposition
output$shifts <- renderTable({
  (d()$summary)
  }, digits = c(0,1,4))


#run insertion function
i <- reactive({
  if((!is.null(input$control) & !is.null(input$sample) & !(input$guide =="")) | (input$example)){
  insertion(a(), d()[1:3])}
  })

#plot +1 insertion
output$insertion <- renderPlot({
  if(i()$er3==F & which(as.numeric(rownames(i()$insertion.summary))==1)){
    
  int <- matrix(as.vector(t(i()$insertion.summary)[,which(as.numeric(rownames(i()$insertion.summary))==1)]),4,1)
  barplot(int, col=c(3,4,1,2), names.arg="+1 insertion", width = 0.1, xlim = c(0, 1))
  legend("top",legend=c("T", "G", "C", "A"), bty="n", pch=15, col=c(2,1,4,3))
  }
  })

#table +1 insertion
output$table_insertion <- renderTable({
  if(i()$er3==F & (as.numeric(rownames(i()$insertion.summary))==1)){
    insert <- matrix(i()$insertion.summary[1,],1,4)
    colnames(insert) <- a()$B
    rownames(insert) <- "+1 insertion"
      insert
  }
  }, digits = 1)

observe({
  c(input$range_sequence, input$range_decom[1], input$range_decom[2], input$maxshift, input$p.threshold)
  if((!is.null(input$control) & !is.null(input$sample) & !(input$guide =="")) | (input$example == TRUE)){
    updateTabsetPanel(session, "maintabset", selected = "Decomposition")}
})
  
})

