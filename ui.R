library(shiny)
library(DT)
library(shinyjs)
library(plotly)

my_content <- mainPanel(
  width = 12,

  tags$style(HTML("
                  /* Existing plot scrolling and pagination CSS */
                  .line-wrapper {
                  display: flex;
                  flex-direction: column;
                  height: 500px;
                  }
                  .line-scroll-container {
                  width: 100%;
                  overflow-x: auto;
                  white-space: nowrap;
                  flex-grow: 1;
                  min-height: 380px;
                  padding-bottom: 10px;
                  }
                  .line-inner {
                  display: inline-block;
                  width: 2200px;
                  min-width: 2200px;
                  }
                  .pagination-btn {
                  background-color: #158cba !important;
                  color: white;
                  border: none;
                  margin: 0 5px;
                  }
                  .pagination-btn:hover {
                  background-color: #106a8c !important;
                  }

                  /* Adjust wellPanel padding/margin for table */
                  .well {
                  padding-top: 8px !important;
                  padding-bottom: 8px !important;
                  margin-bottom: 8px !important;
                  }
                  table.dataTable {
                  margin-bottom: 0px !important;
                  }
                  table.dataTable tbody tr:hover {
                  background-color: transparent !important;
                  }

                  /* Simple dynamic table panel */
                  #table-panel {
                  min-height: 120px;
                  max-height: 600px;
                  height: auto;
                  overflow-y: auto;
                  }
                  ")),

  # Plot + pagination section
  wellPanel(
    style = "margin-top: 15px;",

    downloadLink('downloadLinePlots', 'Download plots', class = "download-link"),

    div(
      class = "line-wrapper",

      div(
        class = "line-scroll-container",
        div(
          class = "line-inner",
          uiOutput("line_plot_ui")
        )
      ),

      div(
        style = "text-align: center; margin: 10px 0;",
        actionButton("prev_page_line", "< Previous", class = "pagination-btn"),
        span(textOutput("line_pagination_info", inline = TRUE),
             style = "margin: 0 10px; font-weight: 500;"),
        actionButton("next_page_line", "Next >", class = "pagination-btn")
      )
    )
  ),

  br(),

  # Table download and DT table output
  wellPanel(
    id = "table-panel",
    style = "margin-top: 15px;",

    downloadLink("dynamicDownloadData", 'Download Table', class = "download-link"),

    DT::DTOutput("testtable")
  ),

  tags$style(type = "text/css", "tfoot {display: none;}")
  )

# EWAS Catalog Content
# ewas_content <- mainPanel(
#   width = 12,
#
#   tags$style(HTML("
#                   .well {
#                   padding-top: 8px !important;
#                   padding-bottom: 8px !important;
#                   margin-bottom: 8px !important;
#                   }
#                   table.dataTable {
#                   margin-bottom: 0px !important;
#                   }
#                   table.dataTable tbody tr:hover {
#                   background-color: transparent !important;
#                   }
#                   #ewas-table-panel {
#                   min-height: 120px;
#                   max-height: 600px;
#                   height: auto;
#                   overflow-y: auto;
#                   }
#                   ")),
#
#   wellPanel(
#     id = "ewas-table-panel",
#     style = "margin-top: 15px;",
#
#     tags$div(
#       style = "text-align: center; padding: 10px 0; margin-bottom: 15px; border-bottom: 1px solid #ccc;",
#       tags$p(
#         HTML('The following table is generated from the EWAS catalog available at
#              <a href="https://ewascatalog.org/" target="_blank">https://ewascatalog.org/</a>'),
#         style = "margin: 0; font-weight: 500; font-size: 14px; color: #2c3e50;"
#       )
#     ),
#
#     downloadLink('ewasDownloadData', 'Download Table', class = "download-link"),
#     DTOutput("ewasTable")
#   )
#   )

cmr_content <- mainPanel(
  width = 12,

  tags$style(HTML("
                  .cmr-wrapper {
                  display: flex;
                  flex-direction: column;
                  height: 500px;
                  }

                  .cmr-scroll-container {
                  width: 100%;
                  overflow-x: auto;
                  white-space: nowrap;
                  flex-grow: 1;
                  min-height: 380px;
                  padding-bottom: 10px;
                  }

                  .cmr-inner {
                  display: inline-block;
                  width: 2600px;
                  min-width: 2600px;
                  }

                  .pagination-btn {
                  background-color: #158cba !important;
                  color: white;
                  border: none;
                  margin: 0 5px;
                  }

                  .pagination-btn:hover {
                  background-color: #106a8c !important;
                  }

                  /* Adjust wellPanel padding/margin for table */
                  .well {
                  padding-top: 8px !important;
                  padding-bottom: 8px !important;
                  margin-bottom: 8px !important;
                  }

                  table.dataTable {
                  margin-bottom: 0px !important;
                  }

                  table.dataTable tbody tr:hover {
                  background-color: transparent !important;
                  }

                  /* Dynamic height and scrolling for the table panel */
                  #cmr-table-panel {
                  min-height: 120px;
                  max-height: 600px;
                  height: auto;
                  overflow-y: auto;
                  }
                  ")),

  wellPanel(
    style = "margin-top: 15px;",

    downloadLink('downloadCMRPlot', 'Download plots', class = "download-link"),

    div(
      class = "cmr-wrapper",

      div(
        class = "cmr-scroll-container",
        div(
          class = "cmr-inner",
          uiOutput("cmr_plots")
        )
      ),

      div(
        style = "text-align: center; margin: 10px 0;",
        actionButton("cmr_prev", "< Previous", class = "pagination-btn"),
        span(
          textOutput("cmr_pagination", inline = TRUE),
          style = "margin: 0 10px; font-weight: 500;"
        ),
        actionButton("cmr_next", "Next >", class = "pagination-btn")
      )
    )
  ),

  br(),

  wellPanel(
    id = "cmr-table-panel",
    style = "margin-top: 15px;",
    downloadLink('cmrdownloadData', 'Download Table', class = "download-link"),
    DT::DTOutput("cmrtable")
  )
  )

mqtl_content <- mainPanel(
  width = 12,

  tags$style(HTML("
                  .well {
                  padding-top: 8px !important;
                  padding-bottom: 8px !important;
                  margin-bottom: 8px !important;
                  }
                  table.dataTable {
                  margin-bottom: 0px !important;
                  }
                  table.dataTable tbody tr:hover {
                  background-color: transparent !important;
                  }
                  #table-panel {
                  min-height: 40px;
                  max-height: 600px;
                  height: auto;
                  overflow-y: auto;
                  }
                  ")),

  # Summary table (index)
  wellPanel(
    id = "mqtl-index-panel",
    style = "margin-top: 15px;",
    downloadLink('mqtlIndexDownloadData', 'Download Table', class = "download-link"),
    DTOutput("mqtlIndexTable")
  ),

  br(),

  # Detailed table
  wellPanel(
    id = "mqtl-long-panel",
    style = "margin-top: 15px;",
    downloadLink('mqtlLongDownloadData', 'Download Table', class = "download-link"),
    DTOutput("mqtlLongTable")
  )
)
user_content <- mainPanel(
  width = 12,

  # --- Advisory Note ---
  tags$div(
    style = "background-color: #f9f9f9; border-left: 4px solid #158cba; padding: 10px 15px; margin-bottom: 15px; font-size: 17px;",
    strong("Note:"), " Please be advised that COMPUTE does not currently adjust for differences in array platforms, cohort demographics, study design variations, batch effects, or technical variability. Potential confounding factors and cell type heterogeneity are also not explicitly accounted for. As such, outputs should be interpreted as exploratory and qualitative. Validation with independent datasets and complementary analyses is recommended."
  ),

  # Section 1: CpG Beta Value Distribution
  wellPanel(
    style = "margin-top: 15px;",

    tags$div(
      style = "text-align: center; padding: 10px 0; margin-bottom: 15px; border-bottom: 2px solid #158cba;",
      tags$h4("Beta Value Distribution",
              style = "margin: 0; font-weight: 600; color: #158cba; font-size: 20px;")
    ),

    # Plot container
    tags$div(
      style = "width: 100%; overflow-x: auto; white-space: nowrap;",
      div(
        style = "display: inline-block; width: 1200px; min-width: 1200px;",
        plotOutput("user_violin_plot", height = "400px", width = "1200px")
      )
    ),

    br(),

    # Pagination controls
    fluidRow(
      column(
        12, align = "center",
        div(
          style = "margin: 10px 0;",
          actionButton("prev_page", "< Previous", class = "pagination-btn"),
          span(
            textOutput("pagination_info", inline = TRUE),
            style = "margin: 0 10px; font-weight: 500;"
          ),
          actionButton("next_page", "Next >", class = "pagination-btn")
        )
      )
    )
  ),

  br(),

  # Section 2: PCA Analysis
  wellPanel(
    style = "margin-top: 15px;",

    tags$div(
      style = "text-align: center; padding: 10px 0; margin-bottom: 15px; border-bottom: 2px solid #158cba;",
      tags$h4("Principal Component Analysis",
              style = "margin: 0; font-weight: 600; color: #158cba;  font-size: 20px;")
    ),

    # Plot container
    plotOutput("pca_plot", height = "600px")
  ),

  br(),

  # Section 3: Scatter Both Plots
  wellPanel(
    tags$div(
      style = "text-align: center; padding: 10px 0; margin-bottom: 15px; border-bottom: 2px solid #158cba;",
      tags$h4("DNA Methylation Based Clocks",
              style = "margin: 0; font-weight: 600; color: #158cba;  font-size: 20px;")
    ),
    plotOutput("scatter_both", height = "700px")
  ),

  br(),

  # Section 4: Cord and Plac Violin Plots
  wellPanel(
    tags$div(
      style = "text-align: center; padding: 10px 0; margin-bottom: 15px; border-bottom: 2px solid #158cba;",
      tags$h4("DNA Methylation Based Cell Types",
              style = "margin: 0; font-weight: 600; color: #158cba;  font-size: 20px;")
    ),

    fluidRow(
      column(12, plotOutput("cord_violin_plot", height = "400px"))
    ),
    br(),
    fluidRow(
      column(12, plotOutput("plac_violin_plot", height = "400px"))
    )
  )
)

# Cover Page Content
cover_page <- div(
  id = "cover-page",
  class = "cover-container",

  tags$style(HTML("
                  .cover-container {
                  position: fixed;
                  top: 0;
                  left: 0;
                  width: 100%;
                  height: 100vh;
                  background: linear-gradient(135deg, #667eea 0%, #158cba 100%);
                  display: flex;
                  flex-direction: column;
                  justify-content: center;
                  align-items: center;
                  z-index: 9999;
                  color: white;
                  text-align: center;
                  }

                  .cover-title {
                  margin-bottom: 20px;
                  animation: fadeInUp 1s ease-out;
                  }

                  .cover-logo {
                  height: 250px;
                  width: auto;
                  animation: fadeInDown 1s ease-out;
                  }

                  .cover-description {
                  font-size: 1.1rem;
                  max-width: 600px;
                  line-height: 1.6;
                  margin-bottom: 40px;
                  opacity: 0.8;
                  animation: fadeInUp 1s ease-out 0.6s both;
                  }

                  .cover-button {
                  background-color: #e67e22;
                  color: white;
                  border: none;
                  padding: 15px 40px;
                  font-size: 1.2rem;
                  border-radius: 50px;
                  cursor: pointer;
                  transition: all 0.3s ease;
                  box-shadow: 0 4px 15px rgba(0,0,0,0.2);
                  animation: fadeInUp 1s ease-out 0.9s both;
                  }

                  .cover-button:hover {
                  background-color: #d76c1a;
                  transform: translateY(-2px);
                  box-shadow: 0 6px 20px rgba(0,0,0,0.3);
                  }

                  @keyframes fadeInUp {
                  from { opacity: 0; transform: translateY(30px); }
                  to { opacity: 1; transform: translateY(0); }
                  }

                  @keyframes fadeInDown {
                  from { opacity: 0; transform: translateY(-30px); }
                  to { opacity: 1; transform: translateY(0); }
                  }

                  .fade-out {
                  animation: fadeOut 0.5s ease-out forwards;
                  }

                  @keyframes fadeOut {
                  to { opacity: 0; visibility: hidden; }
                  }
                  ")),

  img(
    src = "COMPUTE_logo_NoBG.png",
    style = "height:250px; width:auto; border-radius:12px; margin-top:17px;"  # Changed from 170px to 250px
  )
  ,


  tags$button(
    "Enter Application",
    class = "cover-button",
    onclick = "
      document.getElementById('cover-page').classList.add('fade-out');
      setTimeout(function() {
        document.getElementById('cover-page').style.display = 'none';
        document.getElementById('main-content').style.display = 'block';
      }, 500);
    "
  )
)


ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  theme = shinythemes::shinytheme("lumen"),

  # JavaScript for full end-to-end timing measurement
  tags$script(HTML("
                   var startTime = null;

                   // Receive start time from R
                   Shiny.addCustomMessageHandler('startTiming', function(time) {
                   startTime = new Date(time);
                   console.log('Upload started at: ' + startTime);
                   });

                   // Monitor when the plot is actually visible
                   function checkPlotReady() {
                   var plotElement = document.getElementById('user_violin_plot');
                   if (plotElement && plotElement.innerHTML.trim() !== '' && startTime) {
                   var endTime = new Date();
                   var totalSeconds = (endTime - startTime) / 1000;
                   console.log('FULL END-TO-END TIME: ' + totalSeconds.toFixed(2) + ' seconds');

                   // Send back to R
                   Shiny.setInputValue('total_timing', totalSeconds, {priority: 'event'});

                   startTime = null; // Reset to avoid multiple measurements
                   } else {
                   // Keep checking every 100ms
                   setTimeout(checkPlotReady, 100);
                   }
                   }

                   // Start checking after any plot update
                   $(document).on('shiny:value', function(event) {
                   if (event.target.id === 'user_violin_plot') {
                   setTimeout(checkPlotReady, 100);
                   }
                   });
                   ")),

  # Cover page
  cover_page,

  # Main application content (hidden initially)
  div(
    id = "main-content",
    style = "display: none;",

    # top title panel
    div(
      style = "
  background-color: #158cba;
  padding: 8px 20px;
  display: flex;
  align-items: center;
  justify-content: space-between;
  height: 90px;
  border-bottom: 2px solid #106a8c;
  ",

      # COMPUTE logo (left side)
      img(
        src = "COMPUTE_logo_NoBG.png",
        style = "height:300px; width:auto; border-radius:12px; margin-top:25px;"
      ),

      # Partner logos (right side)
      div(
        style = "display: flex; align-items: center; gap: 15px;",

        img(
          src = "BCCHR_cropped.jpeg",
          style = "height: 80px; width: auto; background-color: white; padding: 2px; border-radius: 8px;"
        ),

        img(
          src = "CMMT.jpeg",
          style = "height: 75px; width: auto; background-color: white; padding: 5px; border-radius: 8px;"
        ),

        img(
          src = "ELCHA_cropped.jpg",
          style = "height: 70px; width: auto; background-color: white; padding: 5px; border-radius: 8px;"
        )
      )
    ),

    tags$head(
      tags$style(
        type = "text/css",
        "
        /* Increase base font size */
         body, html { font-size: 17px !important; }
         p { font-size: 17px !important; }
         label { font-size: 17px !important; }
         .btn, button { font-size: 17px !important; }
         .nav-tabs > li > a { font-size: 17px !important; }

        .well { min-height: 100vh; }
        #line_plot { max-height: 600px; }
        .shiny-input-container { width: 100% !important; }
        .section-header { margin: 20px 0 10px; }
        .download-link { display: inline-block; margin-bottom: 10px; }
        .pagination-btn { background-color: #337ab7; color: white; border: none; }
        .pagination-btn:hover { background-color: #286090; }

        .panel-title {
        font-weight: bold;
        text-align: center;
        font-size: 18px;
        }

        .nav-tabs li a[data-value='tab1'],
        .nav-tabs li a[data-value='tab2'],
        .nav-tabs li a[data-value='tab3'],
        .nav-tabs li a[data-value='tab4'] {
        background-color: #6c757d !important;
        color: white !important;
        font-weight: bold !important;
        }

        .nav-tabs li a[data-value='tab5'] {
        background-color: #e67e22 !important;
        color: white !important;
        font-weight: bold !important;
        }

        .btn-file {
        background-color: #158cba !important;
        color: white !important;
        border-color: #158cba !important
        }
        .btn-file:hover {
        background-color: #106a8c !important;
        border-color: #106a8c !important;
        color: white !important;
        }
        "
      )
    ),

    sidebarLayout(
      sidebarPanel(
        width = 3,

        # Tabs with shared input: tab1, tab2, tab3, tab4
        conditionalPanel(
          condition = "input.tabs == 'tab1' || input.tabs == 'tab2' || input.tabs == 'tab3' || input.tabs == 'tab4'",

          p("Enter CpGs or genes to explore DNA methylation patterns in human cord blood and placenta in the CANDLE cohort."),
          tags$hr(),

          tags$textarea(
            id = 'input_cpgs_genes',
            placeholder = 'Type here',
            rows = 3,
            style = "width: 100%;"
          ),

          actionButton("go", "Go", class = "btn btn-primary btn-block"),
          h4("OR", style = "text-align: center; margin: 10px 0;"),

          p("Upload a CSV using the sample template."),

          downloadLink("downloadSampleTemplate", "Download sample template", class = "download-link"),

          fileInput(
            "upload",
            label = "",
            accept = c(".csv"),
            multiple = FALSE
          ),

          tags$hr(),
          div(
            style = "font-size: 14px; color: #2c3e50;",
            p("Dr. Michael Kobor's Lab"),
            p("The University of British Columbia"),
            p("Edwin S.H. Leong Centre for Healthy Aging")
          )
        ),

        # User Data Overlay (Tab5)
        conditionalPanel(
          condition = "input.tabs == 'tab5'",

          p("Upload a CSV file containing your sample metadata, clinical and gestational epigenetic ages, predicted cell type proportions and CpG beta values (Max 50), using the template provided below."),

          downloadLink("downloadUserTemplate", "Download Template", class = "download-link"),
          br(),

          fileInput(
            "upload_userdata",
            label = "",
            accept = ".csv"
          ),

          # DATA SECURITY NOTICE
          div(
            style = "margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-left: 4px solid #158cba; border-radius: 4px;",
            tags$strong(icon("shield-alt"), " Data Privacy Notice"),
            tags$p(
              style = "margin-top: 5px; margin-bottom: 5px; font-size: 12px; color: #6c757d;",
              "Your data is automatically cleared after 1 hour of inactivity or when you close this window.",
            )
          ),

          tags$hr(),

          div(
            style = "font-size: 14px; color: #2c3e50;",
            p("Dr. Michael Kobor's Lab"),
            p("The University of British Columbia"),
            p("Edwin S.H. Leong Centre for Healthy Aging")
          )
        )

      ),

      mainPanel(
        width = 9,

        tabsetPanel(
          id = "tabs",
          tabPanel("CANDLE Data Overview", value = "tab1", my_content),
          tabPanel("CANDLE Co-Methylated Regions (CMR)", value = "tab2", cmr_content),
          #tabPanel("EWAS Catalog", value = "tab3", ewas_content),
          tabPanel("CANDLE Methylation Quantilation Trait Loci (mQTL)", value = "tab4", mqtl_content),
          tabPanel("Visualize Your Data", value = "tab5", user_content)
        )
      )
    )
  )
)
