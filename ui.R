source('packages.R')
ui <- shinydashboard::dashboardPage(
  dashboardHeader(
    title = tags$div(
      tags$img(src = "logo2_1.png", height = "45px"),
      "CloEvoScore",
      style = "color: white;background-color: transparent;"
    ),
    tags$li(
      a(href = "https://github.com/bioinformatic-seragnoli/CloEvoScore.git", target = "_blank", "GitHub", class = "dropdown-item"),
      class = "dropdown"
    ),
    tags$li(
      actionLink("help", "Help", icon = icon("question-circle"), class = "dropdown-item"),
      class = "dropdown"
    ),
    titleWidth = 250
  ),
  
  dashboardSidebar(
    disable = TRUE
  ),
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
        .main-header .logo {
          background-color: #1b4f72 !important;
        }
        .main-header .navbar {
          background-color: #1b4f72 !important;
        }
      "))
    ),
    tags$style(HTML("
  .box {
    margin: 5 px;
    padding: 0px;
  }
  .col-sm-6 {
    padding-right: 5px;
    padding-left: 5px;
  }
  .plot-container {
  text-align: auto;
  margin-top: 0px;
  padding-top: 0px;
  }
   .centered-box {
    margin-left: auto;
    margin-right: auto;
    float: none;
  }
  
  #resbox .box-header {
    background-color: #1b4f72 ;
    color: #fdfefe;
  }
  #resbox .box {
  border: 1px solid #1b4f72 !important;
}
")),
tags$style(HTML(".content-wrapper, .main-header, .main-sidebar {background-color: #e3f2fd;}")),

fluidRow(
  column(
    width = 4,
    div(style = "margin-bottom: 30px;",
        fileInput("gene_calls", "Gene Call Table", accept = c(".csv", ".tsv", ".txt")),
        checkboxInput("toggle_ploidy_purity", "Ploidy and Purity Correction", value = FALSE),
        conditionalPanel(
          condition = "input.toggle_ploidy_purity == true",
          fileInput("ploidy_purity", "Ploidy and Purity Table", accept = c(".csv", ".tsv", ".txt"))
        ),
        checkboxInput("toggle_disease_regions", "Analysis with Disease-Related Gene", value = FALSE),
        conditionalPanel(
          condition = "input.toggle_disease_regions == true",
          fileInput("disease_regions", "Disease-Related Gene Regions", accept = c(".csv", ".tsv", ".txt"))
        ),
        p("If no files are uploaded, the application will use the default datasets"),
        actionButton("run_analysis", "Run Analysis", 
                     icon = icon("play"), 
                     class = "btn btn-primary", 
                     style = "color: #ffffff;")
    )
  ),
  column(
    width = 8,
    uiOutput("results_box"),
    uiOutput("download_button")
  )
)

  )
)


