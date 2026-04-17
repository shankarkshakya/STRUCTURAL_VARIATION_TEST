library(shiny)
library(ggplot2)

ui <- fluidPage(
  
  titlePanel("SV / CNV Explorer"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      fileInput("file", "Upload *.sv_calls.txt"),
      
      uiOutput("chr_ui"),
      
      helpText("Select chromosome after upload"),
      
      actionButton("plotBtn", "Plot")
    ),
    
    mainPanel(
      plotOutput("plot", height = "600px")
    )
  )
)

server <- function(input, output, session) {
  
  # =========================
  # LOAD DATA
  # =========================
  data <- reactive({
    
    req(input$file)
    
    df <- read.table(input$file$datapath, header = FALSE, stringsAsFactors = FALSE)
    
    colnames(df) <- c(
      "chr","start","end",
      "rpm","median","dosage","call"
    )
    
    df$start <- as.numeric(df$start)
    df$end <- as.numeric(df$end)
    df$dosage <- as.numeric(df$dosage)
    df$mid <- (df$start + df$end) / 2
    
    return(df)
  })
  
  # =========================
  # CHROMOSOME DROPDOWN
  # =========================
  output$chr_ui <- renderUI({
    
    req(data())
    
    selectInput(
      "chr",
      "Chromosome",
      choices = unique(data()$chr),
      selected = unique(data()$chr)[1]
    )
  })
  
  # =========================
  # FILTER DATA BY CHR
  # =========================
  filtered <- eventReactive(input$plotBtn, {
    
    req(data(), input$chr)
    
    subset(data(), chr == input$chr)
  })
  
  # =========================
  # PLOT
  # =========================
  output$plot <- renderPlot({
    
    df <- filtered()
    req(df)
    
    ggplot(df, aes(x = mid, y = dosage)) +
      
      geom_line(color = "black") +
      geom_point(aes(color = call), size = 1) +
      
      geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
      
      labs(
        title = paste("Dosage Profile -", input$chr),
        x = "Genomic position (bp)",
        y = "Dosage",
        color = "SV Call"
      ) +
      
      theme_bw()
    
  })
  
}

shinyApp(ui, server)
