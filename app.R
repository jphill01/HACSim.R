#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)

library(devtools)
devtools::install_github("rstudio/shiny-incubator")
library(shinyIncubator)

library(pegas)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Haplotype Accumulation Curve Simulation"),
   
   # Sidebar with numeric inputs for simulation parameters
   sidebarLayout(
      sidebarPanel(
         numericInput("inds",
                     "Number of individuals (N)",
                     value = 2,
                     min = 2
                     ),
         
         numericInput("haps",
                      "Number of haplotypes (H*)",
                      value = 1,
                      min = 1
                      ),
         
         uiOutput("ui"),
         
         numericInput("pops",
                      "Number of (sub)populations (K)",
                      value = 1,
                      min = 1
                      ),
         
         numericInput("perms",
                      "Number of permutations (perms)",
                      value = 10000,
                      min = 10000,
                      max = 10000
                      ),
         
         numericInput("prop",
                      "Proportion of haplotypes to recover (p)",
                      value = 0.95,
                      min = 0.80,
                      max = 1,
                      step = 0.01
                      ), 
         
         fileInput("seqs", "Upload an aligned/trimmed FASTA file (optional)"
                   ),
         
         helpText("Inputted DNA sequences containing missing and/or ambiguous nucleotides may lead to overestimation of the number of observed unique haplotypes.  Consider excluding sequences or alignment sites containing these data. If missing and/or ambiguous bases occur at the ends of sequences, further alignment trimming is an option."
                  ),
         
         actionButton("submit", "Submit"
                      ), 
         width = 5

      ),
      
      # Show plots of the simulated data
      mainPanel(
         plotOutput("distPlot")
      )
   )
)

# Define server logic required to generate plots
server <- function(input, output) {
  
  output$ui <- renderUI(
    matrixInput("probs", "Haplotype frequency distribution (probs)",
                data = data.frame(rep(1/input$haps, input$haps))
                )
  )

}

# Run the application 
shinyApp(ui = ui, server = server)

