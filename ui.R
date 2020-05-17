
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# http://www.rstudio.com/shiny/

  
library(shiny)
                                             
shinyUI(fluidPage(                                                          
  titlePanel(
    list(HTML('<a href="http://www.nki.nl/"><img src="NKIlogo.png" class="img.unframed" align="right" alt="NKI" height="150" width="150"></a>'), "TIDE: Tracking of Indels by DEcomposition"),
    windowTitle="TIDE: Tracking of Indels by DEcomposition"),
    tags$meta("robot.txt"),
  sidebarLayout(
   sidebarPanel(
    tags$head( tags$link(rel="stylesheet", type="text/css", href="styleupdates.css")),
    tags$img(src="logo.png", alt="Tide", width="80%", class="unframed"),
    tags$hr(),
    tags$h4('Upload Data:'),
    
    helpText(tags$span(style="color: black", "Guide sequence:")),
    helpText(tags$small("Submit 20nt guide sequence upstream of PAM ('5-'3)")),
    conditionalPanel(condition = "(!input.example)", tags$textarea(id="guide", rows=1, cols=1, "")), 
    conditionalPanel(condition = "(input.example)", tags$textarea(id="guide", rows=1, cols=1, "CATGCCGAGAGTGATCCCGG")), 
    
    conditionalPanel(condition = "!input.example", fileInput('control', 'Control Sample Chromatogram (.ab1 or .scf)')),
    conditionalPanel(condition = "!input.example", fileInput('sample', 'Test Sample Chromatogram (.ab1 or .scf)')),
        
    checkboxInput('example', tags$span(style="color: black", strong('Load Example Data')), FALSE),
    
    tags$hr(),
    
    tags$h4('Parameters'),
    
    helpText(tags$small("All parameters have default settings but can be adjusted by checking the 'advance settings' box. ")),
    checkboxInput('parameter', 'Advanced settings', FALSE),
    
    conditionalPanel(
      condition = "input.parameter",
    
    tags$br(),
    helpText(tags$span(style="color: black", strong("Alignment window (bp)"))),
    helpText(tags$small("The sequence segment used to align the control and test sample")),
    helpText(tags$small(strong("left boundary"))),
    sliderInput("range_align", "", min = 1, max = 700, value = 100, step= 1),
    helpText(tags$small(strong("right boundary"))),
    helpText(tags$small("automatically set at breaksite - 10bp")),
    
    tags$br(),
            
    helpText(tags$span(style="color: black", strong("Decomposition window (bp)"))),
    helpText(tags$small("The sequence segment used for decomposition.", 
    br(),
    "Default is maximum window possible")),
    
    sliderInput("range_decom", "", min = 1, max = 700, value = c(115, 685), step= 1),
    
    tags$br(),
    
    helpText(tags$span(style="color: black", strong("Indel size range"))),
    helpText(tags$small("Maximum size of insertions and deletions modeled in decomposition")),
    
    sliderInput("maxshift", "", min = 2, max = 50, value = 10, step= 1),
    
    tags$br(),
    
    helpText(tags$span(style="color: black", strong("P-value threshold"))),
    helpText(tags$small("Significance cutoff for decomposition")),
    
    numericInput("p.threshold", "", '0.001', min=0, max=1, step=0.001)
    ),
    
    tags$br(),   
    
    submitButton("Update View"),
    
    tags$hr(),
    tags$small('Created by', tags$a(href="http://research.nki.nl/vansteensellab", "Bas van Steensel lab", target="_blank")),
    tags$br(),
    tags$small('Hosted by', tags$a(href="http://www.nki.nl", target="_blank", "NKI")),
    tags$br(),
    tags$small('TIDE version 0.1'),
    tags$br(),
    tags$small('TIDE supports Firefox6, Chrome4, Safari6, IE10 or higher')
    
  ),
  
  mainPanel(
    tabsetPanel(id="maintabset",      
      tabPanel("Introduction", includeHTML("instructions.html")),
      tabPanel("Quality", 
             conditionalPanel(
                 condition = "!output.fileUploaded",
                 tags$p("A plot will be shown here when the valid sequencing files and guide string have been uploaded.")
             ), 
             
             conditionalPanel(
               condition = "output.fileUploaded",
               tags$h4("Remarks")
             ),
             
             conditionalPanel(
               condition = "output.fileUploaded & !input.parameter",
                 tags$div(htmlOutput("text7"))
               
             ),
                 
               conditionalPanel(
                 condition = "output.fileUploaded & input.parameter",
                 tags$div(textOutput("text1"), style = "color:red")
               ),
                 
              conditionalPanel(
                condition = "output.fileUploaded & output.error1 & input.parameter",
                 tags$div(textOutput("text2"), style = "color:red"),
                tags$br(),
                 tags$div(textOutput("text3"), style = "color:red")
               ),
             
               conditionalPanel(
                 condition = "output.fileUploaded & output.error1",
                 tags$h4("Aberrant sequence signal"), 
                 plotOutput('aberrant_sequence'),
                 tableOutput("meanper"),
                 tags$h4("Checklist quality control"), 
                 helpText(
                   p("Check the following criteria to determine the quality of your data (use the plot)",
                   tags$ul(
                     tags$li("Is there a considerable divergent signal between control and test sample after the breaksite?"),
                      "Sequences of good quality show:",
                     tags$ul(
                       tags$li("in the control sample (black) a low and equally distributed aberrant sequence signal"),
                       tags$li("in the test sample (green) a low signal before the breaksite and a higher signal downstream of the breaksite")
                      ),
                     br(),
                     tags$li("Is the breaksite at expected location?",
                     br(),
                     "The aberrant sequence signal should increase around the expected cut site (blue dotted line)"),
                     br(),
                     tags$li("Does the decomposition window covers a representative sequence?",
                     br(),
                     "For optimal decomposition, the window is set maximal (default), adjust boundaries when the sequence trace is locally of poor quality"
                     )
                   )
                   ) ) 
               )
               
      ),
      tabPanel("Decomposition", 
               conditionalPanel(
                 condition = "!output.fileUploaded",
                 tags$p("The indel spectrum and quantification of indel frequencies will be shown here when the valid sequencing files and guide string have been uploaded.")
               ),
               conditionalPanel(
                 condition = "output.fileUploaded",
                 tags$div(textOutput("text4"), style = "color:red")
               ),
               conditionalPanel(
                 condition = "output.fileUploaded & (output.error1 | output.error2)", 
                 tags$h4("Indel Spectrum"), 
                 plotOutput("decom"),
                 tags$h4("Quantification Indel Frequencies"), 
                 tags$div(textOutput("text5")),
                 tags$br(),
                 tableOutput("shifts")
            )
      ),
      tabPanel("+1 insertion", 
               conditionalPanel(
                 condition = "!output.fileUploaded",
                 tags$p("The prediction of the +1 inserted nucleotide will be shown here when the valid sequencing files and guide string have been uploaded.")
               ), 
               conditionalPanel(
                 condition = "output.fileUploaded & (output.error1 | output.error3)",                
                 tags$h4("Inserted nucleotide probability (%)"), 
                 plotOutput("insertion"),
                 tableOutput("table_insertion")
                 
                 
               )
      )
               
    )
  )
)
)
)