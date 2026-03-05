library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(data.table)
library(tictoc)
library(gridExtra)
library(grid)
library(DT)
library(ggtext)
library(cowplot)
library(promises)
library(future)
library(ggpubr)
library(fst)
library(ewascatalog)
library(ggsignif)
library(conflicted)
library(rsconnect)
library(vegan)
conflict_prefer("melt", "data.table")
conflicts_prefer(DT::renderDataTable)
conflicts_prefer(data.table::dcast)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::lag)




# ------------------------- CALLING ALL THE REQUIRED OBJECTS -------------------------
source('functions_to_source.R')

all_objects <- call_anno()
anno <- all_objects$anno

CANDLE_clocks <- readRDS("Data/CANDLE_clocks.rds")
CANDLE_celltypes <- readRDS("Data/CANDLE_celltypes.rds")

server <- function(input, output, session) {

  in_usable <- reactive(input$input_cpgs_genes)

  tmp <- eventReactive(input$go, {
    filter(anno,
           stringr::str_detect(string = tolower(Gene),
                               pattern = stringr::str_split(subset(in_usable(), !grepl("^[0-9]+$", in_usable())), ",") %>%
                                 unlist %>% stringr::str_trim() %>% tolower() %>% paste0(collapse="|")) |
             stringr::str_detect(string = tolower(CpG),
                                 pattern = stringr::str_split(subset(in_usable(), !grepl("^[0-9]+$", in_usable())), ",") %>%
                                   unlist %>% stringr::str_trim() %>% tolower() %>% paste0(collapse="|")) ) %>%
      dplyr::select(CpG) %>% unlist() %>% as.character() %>% unique()
  })

  # Load methylation data only when CpGs are selected
  methylation_data <- reactive({
    req(validated_input())
    cpgs <- validated_input()

    withProgress(message = 'Loading methylation data...', value = 0, {
      meth_data <- get_methylation_data_chunked(cpgs, anno)

      if (!is.null(meth_data$cord)) {
        cat("Cord rownames:", rownames(meth_data$cord), "\n")
        cat("Cord dimensions:", dim(meth_data$cord), "\n")
      }

      if (!is.null(meth_data$plac)) {
        cat("Plac rownames:", rownames(meth_data$plac), "\n")
        cat("Plac dimensions:", dim(meth_data$plac), "\n")
      }

      # Validate that data was actually loaded
      if (is.null(meth_data$cord) && is.null(meth_data$plac)) {
        validate(need(FALSE, "Could not load methylation data for the selected CpGs"))
      }

      return(meth_data)
    })
  })

  # Create convenient accessors for cord/plac
  cord <- reactive({
    methylation_data()$cord
  })

  plac <- reactive({
    methylation_data()$plac
  })

  # ------------------------- IF USER UPLOADS CSV WITH GENE/CPG LIST -------------------------
  user_list <- reactive({
    if(is.null(input$upload))
      return(NULL)
    data <- data.table::fread(input$upload$datapath, header = FALSE)
    # this prints only annotation table
    anno %>% filter(tolower(CpG) %in% tolower(as.character(data$cpg_or_gene_name)) |
                      tolower(Gene) %in% tolower(as.character(data$cpg_or_gene_name)))
  })

  template_for_upload <- data.table::fread("Data/template_for_upload.csv")

  output$downloadSampleTemplate <- downloadHandler(
    filename = function() {
      'template_for_upload.csv'
    },
    content = function(con) {
      write.csv(template_for_upload, con)
    }
  )



 # ------------------------- CREATING CMR PLOTS - TAB 2 -------------------------
  # Reactive values for pagination
  values <- reactiveValues(cmr_sub_df = NULL, cmr_page = 1, cmr_total_pages = 1)

  # Update CMR subset and pagination info on input change
  observeEvent(validated_input(), {
    tmp_df <- validated_input()

    if (is.null(tmp_df) || length(tmp_df) == 0) {
      values$cmr_sub_df <- NULL
      values$cmr_page <- 1
      values$cmr_total_pages <- 1
      return()
    }

    # Filter CMR dataframe
    sub <- anno %>%
      filter(!Is_in_CMR == "CpG not in a CMR") %>%
      filter(tolower(CpG) %in% tolower(tmp_df)) %>%
      mutate(Is_in_CMR = case_when(
        Is_in_CMR %in% c("CpG in a CMR in only cord") ~ "CMR_Cord",
        Is_in_CMR %in% c("CpG in a CMR in only placenta") ~ "CMR_Placenta",
        Is_in_CMR %in% c("CpG in a CMR in both cord and placenta") ~ "CMR_Both"
      ))

    values$cmr_sub_df <- sub
    values$cmr_page <- 1
    values$cmr_total_pages <- if (nrow(sub) > 0) ceiling(nrow(sub) / 4) else 1
  })

  # Pagination navigation buttons
  observeEvent(input$cmr_next, {
    if (!is.null(values$cmr_total_pages))
      values$cmr_page <- min(values$cmr_page + 1, values$cmr_total_pages)
  })

  observeEvent(input$cmr_prev, {
    if (!is.null(values$cmr_total_pages))
      values$cmr_page <- max(values$cmr_page - 1, 1)
  })

  # Render UI output for CMR plots with 4 plots per page
  output$cmr_plots <- renderUI({
    tmp_df <- validated_input()

    # If no validated input, show note
    if (is.null(tmp_df) || length(tmp_df) == 0) {
      return(div("CpGs/genes not available in COMPUTE"))
    }

    sub <- values$cmr_sub_df

    # If filtering resulted in 0 rows, show note
    if (is.null(sub) || nrow(sub) == 0) {
      return(div("No CpGs/genes found in CMRs."))
    }

    n <- nrow(sub)
    per_page <- 4
    page <- values$cmr_page
    idx <- ((page - 1) * per_page + 1):min(page * per_page, n)

    plot_output_list <- lapply(seq_len(per_page), function(j) {
      plotname <- paste0("cmr_plot_", j)
      plotOutput(plotname, height = "400px", width = "600px")
    })

    div(
      style = "width:100%; overflow-x:auto; white-space:nowrap; display:flex; gap:20px;",
      plot_output_list
    )
  })

  # Render plots for each slot on the current page
  observe({
    sub <- values$cmr_sub_df
    if (is.null(sub) || nrow(sub) == 0) return()

    cord_data <- cord()
    plac_data <- plac()

    n <- nrow(sub)
    per_page <- 4
    page <- values$cmr_page
    idx <- ((page - 1) * per_page + 1):min(page * per_page, n)

    for (j in 1:per_page) {
      local({
        jj <- j
        plotname <- paste0("cmr_plot_", jj)
        if (jj <= length(idx)) {
          i <- idx[jj]
          output[[plotname]] <- renderPlot({
            tissue <- sub$Is_in_CMR[i]
            p <- NULL

            if (tissue == "CMR_Placenta") {
              rn <- stringr::str_split(sub[i, "CMR_Placenta"], ",") %>% unlist() %>% as.character()

              cmr_data <- get_methylation_data_chunked(rn, anno)
              plac_data <- cmr_data$plac

              validate(need(!is.null(plac_data), "Placenta CMR data failed to load."))

              forplot <- plac_data[rownames(plac_data) %in% rn, ] %>% as.data.frame() %>% mutate(CpG = rownames(.))
              p <- cmr_prep(anno, forplot, colour = "#87cefa", tissue = "Placenta")

            } else if (tissue == "CMR_Cord") {
              rn <- stringr::str_split(sub[i, "CMR_Cord"], ",") %>% unlist() %>% as.character()

              cmr_data <- get_methylation_data_chunked(rn, anno)
              cord_data <- cmr_data$cord

              validate(need(!is.null(cord_data), "Cord CMR data failed to load."))

              forplot <- cord_data[rownames(cord_data) %in% rn, ] %>% as.data.frame() %>% mutate(CpG = rownames(.))
              p <- cmr_prep(anno, forplot, colour = "#FADA7A", tissue = "Cord Blood")

            } else { # CMR_Both
              rn_c <- stringr::str_split(sub[i, "CMR_Cord"], ",") %>% unlist() %>% as.character()
              rn_p <- stringr::str_split(sub[i, "CMR_Placenta"], ",") %>% unlist() %>% as.character()

              # Get all unique CpGs needed for this CMR
              all_cmr_cpgs <- unique(c(rn_c, rn_p))

              # Load these CpGs from chunks
              cmr_data <- get_methylation_data_chunked(all_cmr_cpgs, anno)

              # Now extract subsets
              forplot_c <- cmr_data$cord[rownames(cmr_data$cord) %in% rn_c, ] %>%
                as.data.frame()
              forplot_c$CpG <- rownames(cmr_data$cord)[rownames(cmr_data$cord) %in% rn_c]
              forplot_c$Tissue <- "CordBlood"

              forplot_p <- cmr_data$plac[rownames(cmr_data$plac) %in% rn_p, ] %>%
                as.data.frame()
              forplot_p$CpG <- rownames(cmr_data$plac)[rownames(cmr_data$plac) %in% rn_p]
              forplot_p$Tissue <- "Placenta"

              forplot <- rbind(forplot_c, forplot_p)
              forplot$Tissue <- as.factor(forplot$Tissue)
              p <- cmr_prep_common(anno, forplot)
            }

            p + theme(
              panel.background = element_rect(fill = "transparent", color = NA),
              plot.background = element_rect(fill = "transparent", color = NA),
              legend.background = element_rect(fill = "transparent", color = NA),
              legend.box.background = element_rect(fill = "transparent", color = NA)
            )
          }, bg = "transparent", execOnResize = TRUE)

        } else {
          output[[plotname]] <- renderPlot({
            plot.new()
          }, bg = "transparent")
        }
      })
    }
  })

  # Pagination info text
  output$cmr_pagination <- renderText({
    total <- values$cmr_total_pages
    page <- values$cmr_page
    paste("Page", page, "of", total)
  })


    # ------------------------- CREATING CMR TABLE - TAB 2 -------------------------
    validated_input <- reactive({
      # ---- Case 1: CSV upload ----
      if (!is.null(input$upload)) {
        df <- tryCatch({
          read.csv(input$upload$datapath, stringsAsFactors = FALSE)
        }, error = function(e) {
          cat("Error reading CSV:", e$message, "\n")
          return(NULL)
        })

        validate(
          need(!is.null(df) && ncol(df) > 0,
               "Uploaded CSV is empty or invalid.")
        )

        # Extract CpG IDs and Genes from the template columns
        cpg_ids <- if ("CpG_ID" %in% colnames(df)) {
          df$CpG_ID[!is.na(df$CpG_ID) & df$CpG_ID != ""]
        } else {
          character(0)
        }

        genes <- if ("Gene" %in% colnames(df)) {
          df$Gene[!is.na(df$Gene) & df$Gene != ""]
        } else {
          character(0)
        }

        # Combine and clean
        tmp_df <- c(cpg_ids, genes)
        tmp_df <- trimws(tmp_df)
        tmp_df <- tmp_df[tmp_df != ""]
        tmp_df <- unique(tmp_df)

        # DEBUG OUTPUT
        cat("=== CSV UPLOAD DEBUG ===\n")
        cat("CSV columns:", paste(colnames(df), collapse = ", "), "\n")
        cat("CpG IDs extracted:", length(cpg_ids), "\n")
        cat("Genes extracted:", length(genes), "\n")
        cat("Total unique entries:", length(tmp_df), "\n")
        cat("Sample entries:", paste(head(tmp_df, 10), collapse = ", "), "\n")

        validate(
          need(length(tmp_df) > 0,
               "No valid CpG IDs or Gene names found in the uploaded CSV. Please ensure the 'CpG_ID' and/or 'Gene' columns contain data.")
        )

        direct_cpg_test <- tolower(anno$CpG) %in% tolower(tmp_df)
        gene_test <- tolower(anno$Gene) %in% tolower(tmp_df)

        cat("Direct CpG matches found:", sum(direct_cpg_test), "\n")
        cat("Gene matches found:", sum(gene_test), "\n")

        available_cpgs <- anno$CpG[direct_cpg_test]
        gene_matched_cpgs <- anno$CpG[gene_test]
        available <- unique(c(available_cpgs, gene_matched_cpgs))

        cat("Final available count:", length(available), "\n")
        cat("========================\n")

        validate(
          need(length(available) > 0,
               "Note: The CpGs/genes in the uploaded CSV are not available in COMPUTE.")
        )

        return(available)
      }

      # ---- Case 2: Text input ----
      if (isTruthy(input$input_cpgs_genes)) {
        tmp_df <- unlist(strsplit(input$input_cpgs_genes, "[,\\s]+"))
        tmp_df <- trimws(tmp_df)
        tmp_df <- unique(tmp_df)

        # DEBUG OUTPUT
        cat("=== TEXT INPUT DEBUG ===\n")
        cat("Raw input:", input$input_cpgs_genes, "\n")
        cat("Processed input:", paste(tmp_df, collapse = ", "), "\n")
        cat("Input length:", length(tmp_df), "\n")

        # Test the exact matching logic
        direct_cpg_test <- tolower(anno$CpG) %in% tolower(tmp_df)
        gene_test <- tolower(anno$Gene) %in% tolower(tmp_df)

        cat("Direct CpG matches found:", sum(direct_cpg_test), "\n")
        cat("Gene matches found:", sum(gene_test), "\n")

        if ("rbl2" %in% tolower(tmp_df)) {
          cat("Looking for RBL2...\n")
          rbl2_in_anno <- sum(tolower(anno$Gene) == "rbl2")
          cat("RBL2 entries in anno:", rbl2_in_anno, "\n")
        }

        available_cpgs <- anno$CpG[direct_cpg_test]
        gene_matched_cpgs <- anno$CpG[gene_test]
        available <- unique(c(available_cpgs, gene_matched_cpgs))

        cat("Final available count:", length(available), "\n")
        cat("========================\n")

        validate(
          need(length(available) > 0,
               "Note: The entered CpGs/genes are not available in COMPUTE.")
        )

        return(available)
      }

      # ---- Case 3: Nothing provided ----
      validate(
        need(FALSE, "Note: Please type CpGs/genes or upload a CSV.")
      )
    })

    cmrData <- reactive({
      input_cpgs <- validated_input()

      # ---- Filter annotation based on validated CpGs ----
      anno[, c("CpG", "MapInfo", "Chr", "Gene", "CMR_Cord", "CMR_Placenta")] %>%
        dplyr::filter(tolower(CpG) %in% tolower(input_cpgs))
    })


    output$cmrtable <- DT::renderDT({
      df <- cmrData()

      # Handle empty or NULL data to avoid errors
      if (is.null(df) || nrow(df) == 0) {
        empty_df <- data.frame(
          CpG = character(0),
          MapInfo = character(0),
          Chr = character(0),
          Gene = character(0),
          CMR_Cord = character(0),
          CMR_Placenta = character(0),
          stringsAsFactors = FALSE
        )

        display_names <- gsub("_", " ", colnames(empty_df))
        display_names[display_names == "MapInfo"] <- "Genomic Location"


        return(
          DT::datatable(
            empty_df,
            colnames = display_names,
            options = list(
              pageLength = 10,
              lengthMenu = c(5, 10, 20),
              stripe = TRUE,
              autoWidth = TRUE,
              scrollX = TRUE,
              dom = 'lrtip',
              columnDefs = list(
                list(className = 'dt-center', targets = "_all")
              )
            ),
            rownames = FALSE,
            filter = 'none'
          )
        )
      }

      display_names <- gsub("_", " ", colnames(df))
      display_names[display_names == "MapInfo"] <- "Genomic Location"

      DT::datatable(
        df,
        colnames = display_names,
        options = list(
          pageLength = 10,
          lengthMenu = c(5, 10, 20),
          stripe = TRUE,
          autoWidth = FALSE,
          scrollX = TRUE,
          dom = 'lrtip',
          columnDefs = list(
            list(className = 'dt-center', targets = "_all")
          ),
          initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'color': '#158cba'});",
            "}"
          )
        ),
        rownames = FALSE,
        filter = 'none'
      ) %>%
        DT::formatStyle(
          'CpG',
          fontStyle = 'italic'
        ) %>%
        DT::formatStyle(
          'CMR_Cord',
          fontStyle = 'italic'
        )  %>%
        DT::formatStyle(
          'CMR_Placenta',
          fontStyle = 'italic'
        )
    })



    output$cmrdownloadData <- downloadHandler(
      filename = function() {
        paste('cmr-data-', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(cmrData(), con)
      })



  # ------------------------- CREATING LINE PLOTS - TAB 1 -------------------------

    # Initialize pagination reactive values for line plots
    values_line <- reactiveValues(
      current_page_line = 1,
      total_pages_line = 1
    )

    # Update pagination info when validated input changes
    observeEvent(validated_input(), {
      tmp_df_line <- validated_input()

      page_size_line <- 4
      values_line$total_pages_line <- if (length(tmp_df_line) > 0) ceiling(length(tmp_df_line) / page_size_line) else 1
      values_line$current_page_line <- 1
    })

    # Pagination button observers
    observeEvent(input$next_page_line, {
      values_line$current_page_line <- min(values_line$current_page_line + 1, values_line$total_pages_line)
    })
    observeEvent(input$prev_page_line, {
      values_line$current_page_line <- max(values_line$current_page_line - 1, 1)
    })

    # Reactive: CpGs for the current page
    plots_for_page_line <- reactive({
      req(validated_input())
      tmp_df_line <- validated_input()
      page_size_line <- 4
      start_i <- (values_line$current_page_line - 1) * page_size_line + 1
      end_i <- min(values_line$current_page_line * page_size_line, length(tmp_df_line))
      tmp_df_line[start_i:end_i]
    })

    # UI output: line plots for current page
    output$line_plot_ui <- renderUI({
      cpgs_line <- plots_for_page_line()
      if (length(cpgs_line) == 0) return(NULL)

      plot_output_list_line <- lapply(seq_along(cpgs_line), function(i) {
        plotname_line <- paste0("line_plot_", i)
        cpg_name <- cpgs_line[i]

        div(
          style = "display: inline-block; vertical-align: top; margin-right: 20px;",
          div(
            style = "text-align: center; font-style: italic; font-size: 16px; font-family: serif; margin-bottom: 10px; font-weight: normal;",
            cpg_name
          ),
          plotOutput(plotname_line, height = "400px", width = "500px")
        )
      })

      div(
        style = "width: 100%; overflow-x: auto; white-space: nowrap; display: flex; gap: 0px;",
        plot_output_list_line
      )
    })

    # Render plots individually
    observe({
      cpgs_line <- plots_for_page_line()

      cord_data <- cord()
      plac_data <- plac()

      for (i in seq_along(cpgs_line)) {
        local({
          my_i <- i
          plotname_line <- paste0("line_plot_", my_i)
          cpg_line <- cpgs_line[my_i]
          output[[plotname_line]] <- renderPlot({
            data_prep(cpg_line, cord_data, plac_data)
          }, bg = "transparent")
        })
      }
    })

    # Pagination info text
    output$line_pagination_info <- renderText({
      paste("Page", values_line$current_page_line, "of", values_line$total_pages_line)
    })

    output$downloadLinePlots <- downloadHandler(
      filename = function() {
        paste0('line-plots-page-', values_line$current_page_line, '-', Sys.Date(), '.png')
      },
      content = function(file) {
        cpgs_line <- plots_for_page_line()
        n_plots <- length(cpgs_line)

        if (n_plots == 0) return(NULL)

        cord_data <- cord()
        plac_data <- plac()

        # Create PNG with plots side by side
        png(file, width = 500 * n_plots, height = 500, res = 100)
        par(mfrow = c(1, n_plots), mar = c(4, 4, 3, 2))

        for (i in seq_along(cpgs_line)) {
          cpg_line <- cpgs_line[i]
          data_prep(cpg_line, cord_data, plac_data)
          title(main = cpg_line, font.main = 3, cex.main = 1.2)
        }

        dev.off()
      }
    )
  # ------------------------- CREATING ANNOTATION TABLE - TAB 1 -------------------------
    # Load chromHMM annotations
    chromHMM <- readRDS("Data/chromHmm_placenta_cordblood.rds")

    dynamicData <- reactive({
      tmp_df <- validated_input()

      # Base anno selection
      base_df <- anno %>%
        select(CpG, Chr, Gene, Gene_Region, CpG_Island_Class, Is_Variable) %>%
        filter(tolower(CpG) %in% tmp_df)

      # Ensure CpG column type matches for merging
      base_df$CpG <- as.character(base_df$CpG)
      chromHMM$CpG <- as.character(chromHMM$CpG)

      # Merge chromHMM annotations
      base_df <- base_df %>%
        left_join(chromHMM, by = "CpG")

      cord_data <- cord()
      plac_data <- plac()

      # Get correlation stats
      cor_stats <- get_correlation_stats(tmp_df, cord_data, plac_data)
      cor_stats$CpG <- as.character(cor_stats$CpG)

      # Merge correlations
      final_df <- base_df %>%
        left_join(cor_stats, by = "CpG") %>%
        select(CpG, Chr, Gene, Correlation, Is_Variable,
               Gene_Region, CpG_Island_Class, everything())

      return(final_df)
    })


    output$testtable <- DT::renderDT({
      df <- dynamicData()

      # Return empty table if no data
      if (is.null(df) || nrow(df) == 0) {
        empty_df <- data.frame(
          CpG = character(0),
          Chr = character(0),
          Gene = character(0),
          Gene_Region = character(0),
          CpG_Island_Class = character(0),
          Correlation = character(0),
          `Variable In` = character(0),
          stringsAsFactors = FALSE
        )

        display_names <- gsub("_", " ", colnames(empty_df))

        return(DT::datatable(
          empty_df,
          colnames = display_names,
          options = list(
            pageLength = 10,
            lengthMenu = c(5, 10, 20),
            stripe = TRUE,
            autoWidth = TRUE,
            scrollX = TRUE,
            dom = 'lrtip',
            columnDefs = list(list(className = 'dt-center', targets = "_all"))
          ),
          rownames = FALSE,
          filter = 'none'
        ))
      }

      df <- df %>%
        rename("Variable In" = Is_Variable,
               "Spearman Correlation" = Correlation) %>%
        mutate(
          `Variable In` = case_when(
            grepl("cord", `Variable In`, ignore.case = TRUE) & grepl("placenta", `Variable In`, ignore.case = TRUE) ~ "Cord Blood and Placenta",
            grepl("cord", `Variable In`, ignore.case = TRUE) ~ "Cord Blood",
            grepl("placenta", `Variable In`, ignore.case = TRUE) ~ "Placenta",
            TRUE ~ NA_character_
          )
        )

      # Display names
      display_names <- gsub("_", " ", colnames(df))

      DT::datatable(
        df,
        colnames = display_names,
        escape = FALSE,
        options = list(
          pageLength = 10,
          lengthMenu = c(5, 10, 20),
          stripe = TRUE,
          autoWidth = FALSE,
          scrollX = TRUE,
          dom = 'lrtip',
          columnDefs = list(list(className = 'dt-center', targets = "_all")),
          initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'color': '#158cba'});",
            "}"
          )
        ),
        rownames = FALSE,
        filter = 'none'
      ) %>%
        DT::formatStyle('CpG', fontStyle = 'italic') %>%
        DT::formatStyle('Spearman Correlation', `white-space` = "pre-line")
    })


    output$dynamicDownloadData <- downloadHandler(
      filename = function() {
        paste('cpg-annotation-data-', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        df <- dynamicData()

        # Apply the same transformations as in the display table
        if (!is.null(df) && nrow(df) > 0) {
          df <- df %>%
            rename("Variable In" = Is_Variable,
                   "Spearman Correlation" = Correlation) %>%
            mutate(
              `Variable In` = case_when(
                grepl("cord", `Variable In`, ignore.case = TRUE) & grepl("placenta", `Variable In`, ignore.case = TRUE) ~ "Cord Blood and Placenta",
                grepl("cord", `Variable In`, ignore.case = TRUE) ~ "Cord Blood",
                grepl("placenta", `Variable In`, ignore.case = TRUE) ~ "Placenta",
                TRUE ~ NA_character_
              )
            )
        }

        write.csv(df, con, row.names = FALSE)
      }
    )

  ################################# mqtl table
#     mqtl <- read.fst("Data/mqtl_compressed.fst")
#     # Remove the first column
#     mqtl <- mqtl[, -1]
#
#     mqtlData <- reactive({
#       req(validated_input())
#       # Make a temporary copy and rename 'gene' back to 'cpg' if needed
#       df <- mqtl
#       if ("gene" %in% colnames(df)) {
#         df <- df %>% dplyr::rename(cpg = gene)
#       }
#       # Get the validated input (which is already processed and contains CpG names)
#       input_cpgs <- validated_input()
#       # Handle different types of validated input
#       if (is.data.frame(input_cpgs)) {
#         # Case 1: CSV upload
#         if ("cpg_or_gene_name" %in% colnames(input_cpgs)) {
#           input_cpgs <- tolower(trimws(input_cpgs$cpg_or_gene_name))
#         } else {
#           input_cpgs <- tolower(trimws(input_cpgs[[1]]))
#         }
#       } else {
#         # Case 2: Text input
#         input_cpgs <- tolower(input_cpgs)
#       }
#       input_cpgs <- unique(input_cpgs)
#       # Filter rows where 'cpg' matches input, case insensitive
#       filtered <- df %>%
#         filter(tolower(cpg) %in% input_cpgs)
#       if (nrow(filtered) == 0) {
#         # Return empty dataframe with correct columns if no matches
#         return(data.frame(
#           cpg = character(0),
#           SNP_cord_cis = character(0),
#           SNP_cord_trans = character(0),
#           SNP_plac_cis = character(0),
#           SNP_plac_trans = character(0),
#           stringsAsFactors = FALSE
#         ))
#       }
#       filtered
#     })
#     output$mqtlTable <- DT::renderDT({
#       df <- mqtlData()
#       # Map internal column names to display names
#       display_names <- colnames(df)
#       display_names[display_names == "cpg"] <- "CpG"
#       display_names[display_names == "SNP_cord_cis"] <- "Cord Blood cis mQTL"
#       display_names[display_names == "SNP_cord_trans"] <- "Cord Blood trans mQTL"
#       display_names[display_names == "SNP_plac_cis"] <- "Placenta cis mQTL"
#       display_names[display_names == "SNP_plac_trans"] <- "Placenta trans mQTL"
#       # Function to format SNP columns with line breaks every 5 entries
#       formatSNPColumn <- function(text) {
#         if (is.na(text) || text == "" || is.null(text)) return("")
#         entries <- trimws(strsplit(as.character(text), ",")[[1]])
#         if (length(entries) <= 5) return(paste(entries, collapse = ", "))
#         chunks <- split(entries, ceiling(seq_along(entries)/5))
#         formatted_chunks <- sapply(chunks, function(chunk) paste(chunk, collapse = ", "))
#         paste(formatted_chunks, collapse = "<br>")
#       }
#       # Find SNP columns
#       snp_columns <- which(colnames(df) %in% c("SNP_cord_cis", "SNP_cord_trans", "SNP_plac_cis", "SNP_plac_trans"))
#       # Apply formatting to SNP columns
#       for (col in snp_columns) {
#         df[[col]] <- sapply(df[[col]], formatSNPColumn, USE.NAMES = FALSE)
#       }
#       DT::datatable(
#         df,
#         colnames = display_names,
#         escape = FALSE,
#         options = list(
#           pageLength = 10,
#           lengthMenu = c(5, 10, 20),
#           stripe = TRUE,
#           autoWidth = TRUE,
#           scrollX = TRUE,
#           dom = 'lrtip',
#           language = list(emptyTable = "No results found for the entered CpG(s)"),
#           columnDefs = list(
#             list(className = 'dt-center', targets = "_all"),
#             list(targets = snp_columns - 1, className = 'dt-left')
#           ),
#           initComplete = JS(
#             "function(settings, json) {",
#             "$(this.api().table().header()).css({'color': '#158cba'});",
#             "}"
#           )
#         ),
#         rownames = FALSE,
#         filter = 'none'
#       ) %>%
#         DT::formatStyle('cpg', fontStyle = 'italic') %>%
#         DT::formatStyle(columns = colnames(df),
#                         `vertical-align` = 'top')
#     })
#
# output$mqtlDownloadData <- downloadHandler(
#   filename = function() {
#     paste('mqtl-data-', Sys.Date(), '.csv', sep = '')
#   },
#   content = function(con) {
#     write.csv(mqtlData(), con, row.names = FALSE)
#   }
# )


    # Load both mQTL files
    mqtl_index <- read.fst("Data/CANDLE_mqtls_index.fst")
    mqtl_long <- read.fst("Data/CANDLE_mqtls_long.fst")

    mqtl_index <- mqtl_index %>% rename(cpg = CpG)
    mqtl_long  <- mqtl_long  %>% rename(cpg = CpG)


    # Reactive for index table
    mqtlIndexData <- reactive({
      req(validated_input())

      # Get validated input
      input_cpgs <- validated_input()

      if (is.data.frame(input_cpgs)) {
        if ("cpg_or_gene_name" %in% colnames(input_cpgs)) {
          input_cpgs <- tolower(trimws(input_cpgs$cpg_or_gene_name))
        } else {
          input_cpgs <- tolower(trimws(input_cpgs[[1]]))
        }
      } else {
        input_cpgs <- tolower(input_cpgs)
      }
      input_cpgs <- unique(input_cpgs)

      # Filter index table
      filtered <- mqtl_index %>%
        filter(tolower(cpg) %in% input_cpgs)

      if (nrow(filtered) == 0) {
        return(data.frame(
          cpg = character(0),
          n_cord_cis = integer(0),
          n_cord_trans = integer(0),
          n_plac_cis = integer(0),
          n_plac_trans = integer(0),
          stringsAsFactors = FALSE
        ))
      }

      filtered
    })

    mqtlLongData <- reactive({
      req(validated_input())

      input_cpgs <- validated_input()

      if (is.data.frame(input_cpgs)) {
        if ("cpg_or_gene_name" %in% colnames(input_cpgs)) {
          input_cpgs <- tolower(trimws(input_cpgs$cpg_or_gene_name))
        } else {
          input_cpgs <- tolower(trimws(input_cpgs[[1]]))
        }
      } else {
        input_cpgs <- tolower(input_cpgs)
      }
      input_cpgs <- unique(input_cpgs)

      filtered <- mqtl_long %>%
        filter(tolower(cpg) %in% input_cpgs) %>%
        mutate(
          coef = round(coef, 4),
          pval = formatC(pval, format = "e", digits = 2)
          ) %>%
        rename(
          `p-value` = pval
        )

      if (nrow(filtered) == 0) {
        return(data.frame(
          cpg = character(0),
          snp = character(0),
          tissue = character(0),
          type = character(0),
          `p-value` = character(0),
          stringsAsFactors = FALSE
        ))
      }

      filtered
    })


    output$mqtlIndexTable <- DT::renderDT({
      df <- mqtlIndexData()

      # Map column names to display names
      display_names <- colnames(df)
      display_names[display_names == "cpg"] <- "CpG"
      display_names[display_names == "n_cord_cis"] <- "Cord Blood cis mQTL"
      display_names[display_names == "n_cord_trans"] <- "Cord Blood trans mQTL"
      display_names[display_names == "n_plac_cis"] <- "Placenta cis mQTL"
      display_names[display_names == "n_plac_trans"] <- "Placenta trans mQTL"

      DT::datatable(
        df,
        colnames = display_names,
        escape = TRUE,
        options = list(
          pageLength = 20,
          lengthMenu = c(5, 10, 20),
          stripe = TRUE,
          autoWidth = FALSE,
          scrollX = TRUE,
          dom = 'lrtip',
          language = list(emptyTable = "No results found for the entered CpG(s)"),
          columnDefs = list(
            list(className = 'dt-center', targets = "_all")
          ),
          initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'color': '#158cba'});",
            "}"
          )
        ),
        rownames = FALSE,
        filter = 'none'
      ) %>%
        DT::formatStyle('cpg', fontStyle = 'italic')
    })

    output$mqtlLongTable <- DT::renderDT({
      df <- mqtlLongData()

      display_names <- colnames(df)
      display_names[display_names == "cpg"] <- "CpG"
      display_names[display_names == "snp"] <- "SNP"
      display_names[display_names == "tissue"] <- "Tissue"
      display_names[display_names == "type"] <- "Type"
      display_names[display_names == "pval"] <- "p-value"
      display_names[display_names == "coef"] <- "Coef"



      DT::datatable(
        df,
        colnames = display_names,
        escape = TRUE,
        options = list(
          pageLength = 20,
          lengthMenu = c(10, 20, 50, 100),
          stripe = TRUE,
          autoWidth = FALSE,
          scrollX = TRUE,
          dom = 'lrtip',
          language = list(emptyTable = "No results found for the entered CpG(s)"),
          columnDefs = list(
            list(className = 'dt-center', targets = "_all")
          ),
          initComplete = JS(
            "function(settings, json) {",
            "$(this.api().table().header()).css({'color': '#158cba'});",
            "}"
          )
        ),
        rownames = FALSE,
        filter = 'none'
      ) %>%
        DT::formatStyle('cpg', fontStyle = 'italic')
    })

    output$mqtlIndexDownloadData <- downloadHandler(
      filename = function() {
        paste('mqtl-summary-', Sys.Date(), '.csv', sep = '')
      },
      content = function(con) {
        write.csv(mqtlIndexData(), con, row.names = FALSE)
      }
    )

    output$mqtlLongDownloadData <- downloadHandler(
      filename = function() {
        paste('mqtl-detailed-', Sys.Date(), '.csv', sep = '')
      },
      content = function(con) {
        write.csv(mqtlLongData(), con, row.names = FALSE)
      }
    )

################################# EWAS Catalog

#tmp_df <- c("cg00000029","cg00000103","cg00000109","cg00000289","cg00029284")
#
# ewas <- reactive({
#   req(input$input_cpgs_genes)
#
#   # Split by comma or whitespace and remove extra spaces
#   tmp_df <- unlist(strsplit(input$input_cpgs_genes, "[,\\s]+"))
#   tmp_df <- trimws(tmp_df)
#   tmp_df <- unique(tmp_df)
#
#   results_list <- lapply(tmp_df, function(cpg) {
#     tryCatch(
#       ewascatalog(cpg, type = "cpg"),
#       error = function(e) NULL
#     )
#   })
#
#   # Combine into one data frame
#   results <- do.call(rbind, results_list)
#
#   if (is.null(results) || nrow(results) == 0) {
#     return(NULL)
#   }
#
#   # Cleaning up the results table
#   results <- results %>%
#     filter(str_detect(outcome, regex("dna methylation", ignore_case = TRUE))) %>%
#     filter(str_detect(tissue, regex("placenta|cord", ignore_case = TRUE))) %>%
#     select(-c(1,12,13,16,22:26,30,31))
#
#   # Reordering columns
#   results <- results[, c(c(17,12,2,3,1,4:11,13:16,18:20))]
#
#   results
# })
#
#
# output$ewasTable <- DT::renderDT({
#   df <- ewas()
#
#   # If df is NULL or empty, create an empty dataframe with the same columns as normal results
#   if (is.null(df) || nrow(df) == 0) {
#     df <- data.frame(matrix(ncol = 20, nrow = 0))
#     colnames(df) <- c("cpg", "efo", "tissue", "pmid", "date", "consortium",
#                       "trait", "analysis", "source", "outcome", "exposure",
#                       "covariates", "array", "n", "n_studies", "age", "sex",
#                       "beta", "se", "p")
#   }
#
#   # Rename columns for display
#   colnames(df)[colnames(df) == "cpg"] <- "CpG"
#   colnames(df)[colnames(df) == "efo"] <- "Experimental Factor Ontology (EFO)"
#   colnames(df)[colnames(df) == "tissue"] <- "Tissue"
#   colnames(df)[colnames(df) == "pmid"] <- "PMID"
#   colnames(df)[colnames(df) == "date"] <- "Date"
#   colnames(df)[colnames(df) == "consortium"] <- "Consortium"
#   colnames(df)[colnames(df) == "trait"] <- "Trait"
#   colnames(df)[colnames(df) == "analysis"] <- "Analysis"
#   colnames(df)[colnames(df) == "source"] <- "Source"
#   colnames(df)[colnames(df) == "outcome"] <- "Outcome"
#   colnames(df)[colnames(df) == "exposure"] <- "Exposure"
#   colnames(df)[colnames(df) == "covariates"] <- "Covariates"
#   colnames(df)[colnames(df) == "array"] <- "Array"
#   colnames(df)[colnames(df) == "n_studies"] <- "n Studies"
#   colnames(df)[colnames(df) == "age"] <- "Age"
#   colnames(df)[colnames(df) == "sex"] <- "Sex"
#   colnames(df)[colnames(df) == "beta"] <- "Beta"
#   colnames(df)[colnames(df) == "se"] <- "SE"
#
#   display_names <- gsub("_", " ", colnames(df))
#
#   DT::datatable(
#     df,
#     colnames = display_names,
#     escape = TRUE,
#     options = list(
#       pageLength = 10,
#       lengthMenu = c(5, 10, 20),
#       stripe = TRUE,
#       autoWidth = TRUE,
#       scrollX = TRUE,
#       dom = 'lrtip',
#       language = list(
#         emptyTable = "No results found for the selected CpG(s)"
#       ),
#       columnDefs = list(
#         list(className = 'dt-center', targets = "_all"),
#         list(targets = which(colnames(df) == "CpG") - 1, className = 'dt-left')
#       ),
#       initComplete = JS(
#         "function(settings, json) {",
#         "$(this.api().table().header()).css({'color': '#158cba'});",
#         "}"
#       )
#     ),
#     rownames = FALSE,
#     filter = 'none'
#   ) %>%
#     DT::formatStyle('CpG', fontStyle = 'italic')
# })
#

  ################################# Table - Tab 4
  uploaded_data <- reactive({
    req(input$upload_userdata)
    user_file <- input$upload_userdata$datapath
    user_data <- read.csv(user_file, row.names = 1)
    return(user_data)
  })

  filtered_annotations <- reactive({
    user_data <- uploaded_data()
    anno <- call_anno()$anno
    cpg_ids <- rownames(user_data)

    annotations <- anno %>% filter(CpG %in% cpg_ids)

    return(annotations)
  })

  output$user_table <- renderDataTable({
    filtered_annotations()
  })

  output$downloadUserData <- downloadHandler(
    filename = function() {
      paste('filtered_user_data-', Sys.Date(), '.csv', sep = '')
    },
    content = function(con) {
      write.csv(filtered_annotations(), con)
    }
  )


  # ============================================================
  # DATA SECURITY
  # ============================================================

    data_security <- reactiveValues(
    last_activity = Sys.time(),
    upload_temp_file = NULL,
    session_active = TRUE,
    data_cleared = FALSE  # Track if data has been manually cleared
  )

  # Track user activity
  observe({
    # Update last activity time on any input change
    input_values <- reactiveValuesToList(input)
    data_security$last_activity <- Sys.time()
  })

  # Function to clear all user data AND plots
  clear_user_data <- function(reason = "manual") {
    cat("SECURITY: Clearing user uploaded data (reason:", reason, ")...\n")

    temp_file <- isolate(data_security$upload_temp_file)

    # 1. Remove temporary uploaded file if it exists
    if (!is.null(temp_file) && file.exists(temp_file)) {
      tryCatch({
        file.remove(temp_file)
        cat("SECURITY: Removed temporary file:", temp_file, "\n")
      }, error = function(e) {
        cat("SECURITY: Error removing temp file:", e$message, "\n")
      })
    }

    # 2. Clear file input
    shinyjs::reset("upload_userdata")

    # 3. Clear reactive values that store user data
    isolate({
      data_security$upload_temp_file <- NULL
      data_security$data_cleared <- TRUE
    })

    # 4. Reset pagination for violin plots
    isolate({
      values$current_page <- 1
      values$total_pages <- 0
      values$common_cpgs <- NULL
    })

    # 5. Force garbage collection to free memory
    gc()

    cat("SECURITY: User data and plots cleared at", as.character(Sys.time()), "\n")
  }

  # Store temp file path when user uploads
  observeEvent(input$upload_userdata, {
    if (!is.null(input$upload_userdata)) {
      data_security$upload_temp_file <- input$upload_userdata$datapath
      data_security$last_activity <- Sys.time()
      data_security$data_cleared <- FALSE
      cat("SECURITY: File uploaded at", as.character(Sys.time()), "\n")
      cat("SECURITY: Tracking temp file:", data_security$upload_temp_file, "\n")
    }
  })

  # Monitor for 1 hour inactivity
  observe({
    invalidateLater(60000)  # Check every 60 seconds

    if (!is.null(isolate(data_security$upload_temp_file)) && !isolate(data_security$data_cleared)) {
      time_elapsed <- as.numeric(difftime(Sys.time(), isolate(data_security$last_activity), units = "hours"))

      if (time_elapsed >= 1) {
        cat("SECURITY: 1 hour of inactivity detected. Time elapsed:", time_elapsed, "hours\n")
        clear_user_data(reason = "inactivity")

        # Show notification to user
        showNotification(
          "Your uploaded data has been automatically cleared due to inactivity (>1 hour).",
          type = "warning",
          duration = 10
        )
      }
    }
  })

  # Clear data on session end (when user closes browser/tab)
  session$onSessionEnded(function() {
    cat("SECURITY: Session ended at", as.character(Sys.time()), "\n")

    # Access reactive values safely with isolate
    temp_file <- isolate(data_security$upload_temp_file)

    if (!is.null(temp_file) && file.exists(temp_file)) {
      tryCatch({
        file.remove(temp_file)
        cat("SECURITY: Removed temporary file:", temp_file, "\n")
      }, error = function(e) {
        cat("SECURITY: Error removing temp file:", e$message, "\n")
      })
    }

    cat("SECURITY: Session cleanup complete\n")
  })

  # Manual clear button handler
  observeEvent(input$clear_data_btn, {
    clear_user_data(reason = "manual")
    showNotification(
      "Your uploaded data and visualizations have been cleared.",
      type = "message",
      duration = 5
    )
  })

  ############################## Beta Plots
  user_plot_methylation_data <- reactive({
    req(input$upload_userdata)

    if (data_security$data_cleared) {
      return(NULL)
    }

    data_security$last_activity <- Sys.time()

    user_data <- tryCatch({
      read.csv(input$upload_userdata$datapath, stringsAsFactors = FALSE)
    }, error = function(e) {
      return(NULL) # Return NULL if read fails
    })

    validate(need(!is.null(user_data), "Failed to read uploaded CSV file."))

    cpg_pattern <- "^cg[0-9]{8}$"
    cpg_columns <- grep(cpg_pattern, colnames(user_data), value = TRUE)

    validate(need(length(cpg_columns) > 0, "No valid CpG columns found in the uploaded file."))

    cpgs_to_load <- cpg_columns

    withProgress(message = 'Loading reference data for comparison...', value = 0.4, {
      meth_data <- get_methylation_data_chunked(cpgs_to_load, anno)

      validate(
        need(!is.null(meth_data$cord) && !is.null(meth_data$plac),
             "Could not load reference data for the CpGs found in your file.")
      )

      return(meth_data)
    })
  })


  plot_cord_plac_user <- function(cpgs, cord, plac, user_data, user_tissue_label = "User", page = 1, plots_per_page = 4) {

    # Helper function to prepare dataframes
    prep_df <- function(df) {
      df <- as.data.frame(df, stringsAsFactors = FALSE)
      if (!"CpG_ID" %in% colnames(df)) {
        df$CpG_ID <- rownames(df)
        rownames(df) <- NULL
      }
      df[!duplicated(df$CpG_ID), ]
    }

    # Helper function to extract CpG values
    extract_cpg_values <- function(df, cpg, exclude_cols = "CpG_ID") {
      if (!(cpg %in% df$CpG_ID)) return(numeric(0))
      row_idx <- which(df$CpG_ID == cpg)[1]
      vals <- as.numeric(as.matrix(df[row_idx, setdiff(names(df), exclude_cols)]))

      # Check for and handle values outside [0,1] range
      if (any(vals > 1, na.rm = TRUE)) {
        warning(paste("Found beta values > 1 for CpG", cpg, "- clamping to 1"))
        vals[vals > 1] <- 1
      }
      if (any(vals < 0, na.rm = TRUE)) {
        warning(paste("Found beta values < 0 for CpG", cpg, "- clamping to 0"))
        vals[vals < 0] <- 0
      }

      return(vals)
    }

    # Custom function to calculate p-value and Cohen's d
    custom_stat_function <- function(data, group1, group2) {
      data1 <- data$BetaValue[data$Tissue == group1]
      data2 <- data$BetaValue[data$Tissue == group2]

      # Remove NA values
      data1 <- data1[!is.na(data1)]
      data2 <- data2[!is.na(data2)]

      if (length(data1) > 1 && length(data2) > 1) {
        # Calculate t-test
        t_result <- t.test(data1, data2)
        p_val <- t_result$p.value

        # Calculate Cohen's d
        pooled_sd <- sqrt(((length(data1) - 1) * var(data1) + (length(data2) - 1) * var(data2)) /
                            (length(data1) + length(data2) - 2))
        cohens_d <- abs(mean(data1) - mean(data2)) / pooled_sd

        # Format p-value
        if (p_val < 0.001) {
          p_text <- "p<0.001"
        } else if (p_val < 0.01) {
          p_text <- sprintf("p=%.3f", p_val)
        } else {
          p_text <- sprintf("p=%.2f", p_val)
        }

        # Format Cohen's d
        d_text <- sprintf("d=%.2f", cohens_d)

        return(paste(p_text, d_text, sep = ", "))
      } else {
        return("Insufficient data")
      }
    }

    # Prepare data
    cord_df <- prep_df(cord)
    plac_df <- prep_df(plac)
    user_cord <- user_data[user_data$TissueType == "Cord Blood", ]
    user_plac <- user_data[user_data$TissueType == "Placenta", ]
    user_cord <- user_cord[!duplicated(user_cord$CpG_ID), ]
    user_plac <- user_plac[!duplicated(user_plac$CpG_ID), ]

    # Build plot data for all CpGs
    all_plot_data <- lapply(cpgs, function(cpg) {
      # Skip if CpG not found in any dataset
      if (!(cpg %in% c(cord_df$CpG_ID, plac_df$CpG_ID, user_cord$CpG_ID, user_plac$CpG_ID))) return(NULL)

      # Extract values
      vals <- list(
        cord = extract_cpg_values(cord_df, cpg),
        plac = extract_cpg_values(plac_df, cpg),
        user_cord = extract_cpg_values(user_cord, cpg, c("CpG_ID", "TissueType")),
        user_plac = extract_cpg_values(user_plac, cpg, c("CpG_ID", "TissueType"))
      )

      # Apply sampling limits
      vals$cord <- if (length(vals$cord) > 500) sample(vals$cord, 500) else vals$cord
      vals$plac <- if (length(vals$plac) > 500) sample(vals$plac, 500) else vals$plac
      vals$user_cord <- if (length(vals$user_cord) > 500) sample(vals$user_cord, 500) else vals$user_cord
      vals$user_plac <- if (length(vals$user_plac) > 500) sample(vals$user_plac, 500) else vals$user_plac

      # Skip if no data
      if (sum(lengths(vals)) == 0) return(NULL)

      # Create data frame for this CpG
      data.frame(
        BetaValue = unlist(vals),
        Tissue = c(rep("CANDLE_CordBlood", length(vals$cord)),
                   rep("CANDLE_Placenta", length(vals$plac)),
                   rep(paste0(user_tissue_label, "_CordBlood"), length(vals$user_cord)),
                   rep(paste0(user_tissue_label, "_Placenta"), length(vals$user_plac))),
        CpG = cpg,
        stringsAsFactors = FALSE
      )
    })

    # Remove NULL entries
    all_plot_data <- all_plot_data[!sapply(all_plot_data, is.null)]
    names(all_plot_data) <- sapply(all_plot_data, function(x) x$CpG[1])

    if (length(all_plot_data) == 0) {
      return(list(ggplot() + geom_text(aes(x = 0, y = 0), label = "No matching CpG data found", size = 5) + theme_void()))
    }

    # Pagination
    unique_cpgs <- names(all_plot_data)
    total_pages <- ceiling(length(unique_cpgs) / plots_per_page)
    page <- max(1, min(page, total_pages))

    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, length(unique_cpgs))
    page_cpgs <- unique_cpgs[start_idx:end_idx]

    # Create individual plots for each CpG on this page
    fill_colors <- setNames(
      c("#FADA7A", "#87cefa", "#bfaf3b", "#1E6091"),
      c("CANDLE_CordBlood", "CANDLE_Placenta",
        paste0(user_tissue_label, "_CordBlood"), paste0(user_tissue_label, "_Placenta"))
    )

    plots <- lapply(page_cpgs, function(cpg) {
      cpg_data <- all_plot_data[[cpg]]

      # Create comparisons list for statistical tests
      comparisons <- list()
      custom_labels <- character()

      # Check if both CANDLE and User data exist for Cord Blood
      if (any(cpg_data$Tissue == "CANDLE_CordBlood") &&
          any(cpg_data$Tissue == paste0(user_tissue_label, "_CordBlood"))) {
        comparisons <- append(comparisons, list(c("CANDLE_CordBlood", paste0(user_tissue_label, "_CordBlood"))))
        custom_labels <- c(custom_labels, custom_stat_function(cpg_data, "CANDLE_CordBlood", paste0(user_tissue_label, "_CordBlood")))
      }

      # Check if both CANDLE and User data exist for Placenta
      if (any(cpg_data$Tissue == "CANDLE_Placenta") &&
          any(cpg_data$Tissue == paste0(user_tissue_label, "_Placenta"))) {
        comparisons <- append(comparisons, list(c("CANDLE_Placenta", paste0(user_tissue_label, "_Placenta"))))
        custom_labels <- c(custom_labels, custom_stat_function(cpg_data, "CANDLE_Placenta", paste0(user_tissue_label, "_Placenta")))
      }

      # Base plot with constrained violin plots
       p <- ggplot(cpg_data, aes(x = Tissue, y = BetaValue, fill = Tissue)) +
        geom_violin(trim = FALSE, adjust = 1.5, alpha = 0.7) +
        geom_jitter(aes(color = Tissue), width = 0.05, alpha = 0.3, size = 0.8) +
        scale_fill_manual(values = fill_colors) +
        scale_color_manual(values = fill_colors) +
        scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1, 0.2)) +
        labs(title = cpg, y = "Beta Value") +
        theme_minimal() +
        theme(
          text = element_text(size = 14),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          plot.title = element_text(face = "italic", size = 14, hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "pt")
        )

      # Add statistical comparisons using geom_signif
      if (length(comparisons) > 0) {

        # Add each comparison separately with explicit positioning
        if (length(comparisons) >= 1) {
          p <- p + geom_signif(
            comparisons = list(comparisons[[1]]),
            annotations = custom_labels[1],
            y_position = 1.08,
            textsize = 3.0,
            vjust = -0.2
          )
        }

        if (length(comparisons) >= 2) {
          p <- p + geom_signif(
            comparisons = list(comparisons[[2]]),
            annotations = custom_labels[2],
            y_position = 1.18,
            textsize = 3.0,
            vjust = -0.2
          )
        }
      }

      return(p)
    })

    return(plots)
  }

  values <- reactiveValues(common_cpgs = NULL, current_page = 1, total_pages = 0)

  observeEvent(input$upload_userdata, { values$current_page <- 1 })
  observeEvent(input$prev_page, { values$current_page <- max(1, values$current_page - 1) })
  observeEvent(input$next_page, { values$current_page <- min(values$total_pages, values$current_page + 1) })

  output$user_violin_plot <- renderPlot({
    req(input$upload_userdata)

    withProgress(message = "Generating plots...", value = 0, {

      reference_data <- user_plot_methylation_data()
      cord_loaded <- reference_data$cord
      plac_loaded <- reference_data$plac

      # Safely try to read the uploaded CSV file
      user_data <- tryCatch({
        read.csv(input$upload_userdata$datapath, stringsAsFactors = FALSE)
      }, error = function(e) {
        return(NULL)
      })

      incProgress(0.1, detail = "Validating input...")

      # Function to display error in a neutral box, vertically centered
      error_plot <- function(msg) {
        validate(
          need(FALSE, msg)
        )
      }



      # Error Handling for Beta Values Violin Plot
      if (is.null(user_data)) {
        return(error_plot("Note: Failed to read uploaded CSV file."))
      }

      # Check for required columns
      required_cols <- c("TissueType")
      if (!all(required_cols %in% colnames(user_data))) {
        return(error_plot("Note: Uploaded file must contain 'TissueType' column."))
      }

      incProgress(0.2, detail = "Processing data format...")

      # Identify CpG columns (assuming they start with "cg" followed by 8 digits)
      cpg_pattern <- "^cg[0-9]{8}$"
      cpg_columns <- grep(cpg_pattern, colnames(user_data), value = TRUE)

      if (length(cpg_columns) == 0) {
        return(error_plot("Note: CpG not available in COMPUTE."))
      }

      # Validate CpG values - Enhanced validation
      cpg_data <- user_data[, cpg_columns, drop = FALSE]

      # Check for non-numeric values first
      non_numeric_cols <- sapply(cpg_data, function(x) !is.numeric(x))
      if (any(non_numeric_cols)) {
        problem_cols <- names(non_numeric_cols)[non_numeric_cols]
        return(error_plot(paste("Note: Non-numeric values found in CpG columns:", paste(problem_cols[1:min(3, length(problem_cols))], collapse = ", "))))
      }

      # Check for values outside 0-1 range (ignore NAs completely)
      range_violations <- sapply(cpg_data, function(x) {
        valid_values <- x[!is.na(x)]
        if (length(valid_values) == 0) return(FALSE)
        any(valid_values < 0 | valid_values > 1)
      })

      if (any(range_violations)) {
        problem_cols <- names(range_violations)[range_violations]
        return(error_plot(paste("Note: Values outside 0-1 range found in CpG columns:", paste(problem_cols[1:min(3, length(problem_cols))], collapse = ", "))))
      }

      # Identify rows that have at least one non-NA CpG value
      valid_cpg_rows <- apply(cpg_data, 1, function(row) {
        any(!is.na(row))
      })

      num_valid_cpg_samples <- sum(valid_cpg_rows)

      user_data_with_cpg <- user_data[valid_cpg_rows, ]

      tissue_counts <- table(user_data_with_cpg$TissueType)
      if (any(tissue_counts == 0)) {
        return(error_plot("Note: TissueType column contains empty values."))
      }

      tissue_sample_counts <- table(user_data_with_cpg$TissueType)
      if (any(tissue_sample_counts < 1)) {
        return(error_plot("Note: Each tissue type must have at least one sample."))
      }

      # Check for reasonable number of CpG samples (not total rows)
      # if (num_valid_cpg_samples > 50) {
      #   return(error_plot(paste("Note: Too many CpG samples. Found", num_valid_cpg_samples, "samples with CpG data. Maximum 50 supported.")))
      # }
      #

      if (num_valid_cpg_samples == 0) {
        return(error_plot("Note: No samples found with valid CpG data."))
      }

      # # Check for reasonable number of CpGs (optional upper limit)
      # if (length(cpg_columns) > 50) {
      #   return(error_plot("Note: Too many CpG columns. Maximum 50 CpGs supported."))
      # }

      incProgress(0.3, detail = "Transforming data format...")

      # Transform data: Convert from wide format (CpGs as columns) to long format (CpGs as rows)
      user_data_filtered <- user_data_with_cpg

      # Create sample identifiers if they don't exist
      if (!"Sample_ID" %in% colnames(user_data_filtered)) {
        user_data_filtered$Sample_ID <- paste0("Sample_", seq_len(nrow(user_data_filtered)))
      }

      user_long <- reshape2::melt(user_data_filtered,
                                  id.vars = c("Sample_ID", "TissueType"),
                                  measure.vars = cpg_columns,
                                  variable.name = "CpG_ID",
                                  value.name = "BetaValue")

      # Convert to the expected format for plot_cord_plac_user function
      user_transformed <- reshape2::dcast(user_long,
                                          CpG_ID + TissueType ~ Sample_ID,
                                          value.var = "BetaValue")

      # Check if any user CpGs are present in reference tissues
      prep_df <- function(df) {
        df <- as.data.frame(df, stringsAsFactors = FALSE)
        if (!"CpG_ID" %in% colnames(df)) {
          df$CpG_ID <- rownames(df)
          rownames(df) <- NULL
        }
        df$CpG_ID <- as.character(df$CpG_ID)
        df
      }

      cord_df <- prep_df(cord_loaded)
      plac_df <- prep_df(plac_loaded)

      cpgs_in_user <- unique(user_transformed$CpG_ID)
      cpgs_in_cord <- unique(cord_df$CpG_ID)
      cpgs_in_plac <- unique(plac_df$CpG_ID)
      cpgs_in_ref <- union(cpgs_in_cord, cpgs_in_plac)

      # Identify which CpGs are missing from reference data
      missing_cpgs <- setdiff(cpgs_in_user, cpgs_in_ref)
      common_cpgs <- intersect(cpgs_in_user, cpgs_in_ref)

      if (length(common_cpgs) == 0) {
        if (length(missing_cpgs) <= 5) {
          return(error_plot(paste("Note: None of the uploaded CpG IDs are available in COMPUTE reference data. Missing CpGs:", paste(missing_cpgs, collapse = ", "))))
        } else {
          return(error_plot(paste("Note: None of the uploaded CpG IDs are available in COMPUTE reference data. Missing", length(missing_cpgs), "CpGs including:", paste(missing_cpgs[1:3], collapse = ", "), "...")))
        }
      }

      # Handle partially missing CpGs - provide user feedback
      if (length(missing_cpgs) > 0) {
        warning_msg <- if (length(missing_cpgs) <= 5) {
          paste("Note:", length(missing_cpgs), "CpG(s) not found in reference data and will be skipped:", paste(missing_cpgs, collapse = ", "))
        } else {
          paste("Note:", length(missing_cpgs), "CpG(s) not found in reference data and will be skipped. First few:", paste(missing_cpgs[1:3], collapse = ", "), "...")
        }


        cat("Warning:", warning_msg, "\n")
      }

      incProgress(0.5, detail = "Preparing plot data...")

      # Separate user data by tissue type
      user_cord <- user_transformed[user_transformed$TissueType == "Cord Blood", ]
      user_plac <- user_transformed[user_transformed$TissueType == "Placenta", ]

      has_user_cord <- nrow(user_cord) > 0
      has_user_plac <- nrow(user_plac) > 0

      if (!has_user_cord && !has_user_plac) {
        return(error_plot("No valid tissue types found in user data."))
      }

      # Filter to common CpGs only
      if (has_user_cord) {
        user_cord <- user_cord[user_cord$CpG_ID %in% common_cpgs, ]
      }
      if (has_user_plac) {
        user_plac <- user_plac[user_plac$CpG_ID %in% common_cpgs, ]
      }

      # Combine user data back
      user_final <- rbind(user_cord, user_plac)

      values$common_cpgs <- common_cpgs
      values$total_pages <- ceiling(length(common_cpgs) / 4)

      incProgress(0.7, detail = "Building plot...")

      plots_list <- plot_cord_plac_user(common_cpgs, cord_df, plac_df, user_final, page = values$current_page, plots_per_page = 4)

      # Arrange plots in a single row (4 columns, 1 row)
      if (length(plots_list) == 1) {
        plots_list[[1]]
      } else {
        do.call(gridExtra::grid.arrange, c(plots_list, nrow = 1))
      }

    })
  }, bg = "transparent")

  output$pagination_info <- renderText({
    req(values$total_pages > 0)
    paste("Page", values$current_page, "of", values$total_pages)
  })


  ##################### PCA plots
  plot_pca_comparison <- function(common_cpgs, cord, plac, user_data, user_tissue_label = "User") {

    # Helper function to prepare dataframes (same as violin plot)
    prep_df <- function(df) {
      df <- as.data.frame(df, stringsAsFactors = FALSE)
      if (!"CpG_ID" %in% colnames(df)) {
        df$CpG_ID <- rownames(df)
        rownames(df) <- NULL
      }
      df[!duplicated(df$CpG_ID), ]
    }

    # Prepare reference data
    cord_df <- prep_df(cord)
    plac_df <- prep_df(plac)

    # Prepare user data - FIXED: Remove TissueType column before processing
    user_cord <- user_data[user_data$TissueType == "Cord Blood", ]
    user_plac <- user_data[user_data$TissueType == "Placenta", ]
    user_cord <- user_cord[!duplicated(user_cord$CpG_ID), ]
    user_plac <- user_plac[!duplicated(user_plac$CpG_ID), ]

    if (nrow(user_cord) > 0 && "TissueType" %in% colnames(user_cord)) {
      user_cord$TissueType <- NULL
    }
    if (nrow(user_plac) > 0 && "TissueType" %in% colnames(user_plac)) {
      user_plac$TissueType <- NULL
    }

    # Filter to common CpGs only
    cord_df <- cord_df[cord_df$CpG_ID %in% common_cpgs, ]
    plac_df <- plac_df[plac_df$CpG_ID %in% common_cpgs, ]
    user_cord <- user_cord[user_cord$CpG_ID %in% common_cpgs, ]
    user_plac <- user_plac[user_plac$CpG_ID %in% common_cpgs, ]

    # Check if we have enough data
    if (nrow(cord_df) == 0 && nrow(plac_df) == 0) {
      return(ggplot() +
               geom_text(aes(x = 0, y = 0), label = "No CANDLE reference data available for PCA", size = 5) +
               theme_void())
    }

    if (nrow(user_cord) == 0 && nrow(user_plac) == 0) {
      return(ggplot() +
               geom_text(aes(x = 0, y = 0), label = "No user data available for PCA", size = 5) +
               theme_void())
    }

    # Create matrices for each dataset
    create_matrix <- function(df, prefix, tissue_type) {
      if (nrow(df) == 0) return(NULL)

      # Get sample columns (exclude CpG_ID and TissueType)
      sample_cols <- setdiff(names(df), c("CpG_ID", "TissueType"))

      if (length(sample_cols) == 0) return(NULL)

      # Create matrix with CpGs as rows, samples as columns
      mat <- as.matrix(df[, sample_cols, drop = FALSE])
      rownames(mat) <- df$CpG_ID

      colnames(mat) <- paste0(prefix, "_", colnames(mat))

      # Remove columns with all NA values
      mat <- mat[, !apply(is.na(mat), 2, all), drop = FALSE]

      if (ncol(mat) > 0) {
        cat("Created", prefix, "matrix:", nrow(mat), "CpGs x", ncol(mat), "samples\n")
      }

      return(mat)
    }

    # Create matrices for each group with explicit tissue type labels
    candle_cord_mat <- create_matrix(cord_df, "CANDLE_Cord", "Cord")
    candle_plac_mat <- create_matrix(plac_df, "CANDLE_Placenta", "Placenta")
    user_cord_mat <- create_matrix(user_cord, paste0(user_tissue_label, "_Cord"), "Cord")
    user_plac_mat <- create_matrix(user_plac, paste0(user_tissue_label, "_Placenta"), "Placenta")

    # Combine all matrices
    all_matrices <- list(candle_cord_mat, candle_plac_mat, user_cord_mat, user_plac_mat)
    all_matrices <- all_matrices[!sapply(all_matrices, is.null)]

    if (length(all_matrices) == 0) {
      return(ggplot() +
               geom_text(aes(x = 0, y = 0), label = "No valid data for PCA analysis", size = 5) +
               theme_void())
    }

    # Find common CpGs across all matrices
    common_cpgs_matrices <- Reduce(intersect, lapply(all_matrices, rownames))

    if (length(common_cpgs_matrices) < 2) {
      return(ggplot() +
               geom_text(aes(x = 0, y = 0),
                         label = paste("Insufficient CpGs for PCA analysis.\nNeed at least 2 CpGs, found:", length(common_cpgs_matrices)),
                         size = 4) +
               theme_void())
    }

    # Filter matrices to common CpGs and combine
    all_matrices_filtered <- lapply(all_matrices, function(mat) {
      if (is.null(mat)) return(NULL)
      mat[common_cpgs_matrices, , drop = FALSE]
    })

    # Combine into single matrix
    beta_combined <- do.call(cbind, all_matrices_filtered)

    cat("Combined matrix sample names (first 5):", head(colnames(beta_combined), 5), "\n")
    cat("Total samples in combined matrix:", ncol(beta_combined), "\n")

    # Remove any remaining NA values by column
    complete_cols <- apply(beta_combined, 2, function(x) !any(is.na(x)))
    beta_combined <- beta_combined[, complete_cols, drop = FALSE]

    if (ncol(beta_combined) < 3) {
      return(ggplot() +
               geom_text(aes(x = 0, y = 0), label = "Insufficient samples for PCA analysis", size = 5) +
               theme_void())
    }

    # Transpose for PCA (samples as rows, CpGs as columns)
    beta_t <- t(beta_combined)

    samples <- rownames(beta_t)
    meta_df <- data.frame(Sample = samples, stringsAsFactors = FALSE) %>%
      mutate(
        Tissue = case_when(
          grepl("_Cord_", Sample, fixed = FALSE) ~ "Cord",
          grepl("_Placenta_", Sample, fixed = FALSE) ~ "Placenta",
          TRUE ~ "Unknown"
        ),
        Source = case_when(
          grepl("^CANDLE_", Sample) ~ "CANDLE",
          grepl(paste0("^", user_tissue_label, "_"), Sample) ~ user_tissue_label,
          TRUE ~ "Unknown"
        ),
        Group = paste(Source, Tissue, sep = "_")
      )

    unknown_samples <- meta_df %>% filter(Tissue == "Unknown" | Source == "Unknown")
    if (nrow(unknown_samples) > 0) {
      cat("WARNING: Found samples with unknown classification:\n")
      print(unknown_samples)
    }

    cat("\nSample distribution by group:\n")
    print(table(meta_df$Group))

    # Perform PCA
    pca_result <- tryCatch({
      prcomp(beta_t, center = TRUE, scale. = TRUE)
    }, error = function(e) {
      cat("PCA error:", e$message, "\n")
      return(NULL)
    })

    if (is.null(pca_result)) {
      return(ggplot() +
               geom_text(aes(x = 0, y = 0), label = "PCA analysis failed", size = 5) +
               theme_void())
    }

    # Create PCA dataframe
    pca_df <- as.data.frame(pca_result$x)
    pca_df$Sample <- rownames(pca_df)
    pca_df <- left_join(pca_df, meta_df, by = "Sample")

    # PERMANOVA analysis
    dist_mat <- dist(beta_t, method = "euclidean")

    # Global PERMANOVA
    adonis_result <- adonis2(dist_mat ~ Group, data = meta_df)
    r2_val <- round(adonis_result$R2[1], 3)
    p_val <- format.pval(adonis_result$`Pr(>F)`[1], digits = 2, eps = .001)

    # Pairwise PERMANOVA for matched tissues
    pairwise_matched_adonis <- function(dist.mat, grouping) {
      results <- data.frame(Group1 = character(), Group2 = character(),
                            R2 = numeric(), p.value = numeric(), stringsAsFactors = FALSE)

      # Define pairs based on available groups
      available_groups <- unique(grouping)
      pairs <- list()

      if (all(c("CANDLE_Cord", paste0(user_tissue_label, "_Cord")) %in% available_groups)) {
        pairs <- append(pairs, list(c("CANDLE_Cord", paste0(user_tissue_label, "_Cord"))))
      }
      if (all(c("CANDLE_Placenta", paste0(user_tissue_label, "_Placenta")) %in% available_groups)) {
        pairs <- append(pairs, list(c("CANDLE_Placenta", paste0(user_tissue_label, "_Placenta"))))
      }

      for (pair in pairs) {
        sel <- grouping %in% pair
        if (sum(sel) < 3) next  # Need at least 3 samples for PERMANOVA

        sub.dist <- as.dist(as.matrix(dist.mat)[sel, sel])
        sub.group <- grouping[sel]

        tryCatch({
          ad <- adonis2(sub.dist ~ sub.group)
          results <- rbind(results,
                           data.frame(Group1 = pair[1], Group2 = pair[2],
                                      R2 = round(ad$R2[1], 3),
                                      p.value = ad$`Pr(>F)`[1]))
        }, error = function(e) {
          # Skip failed comparisons
        })
      }

      if (nrow(results) > 0) {
        results$p.adj <- p.adjust(results$p.value, method = "BH")
      }

      return(results)
    }

    pairwise_results <- pairwise_matched_adonis(dist_mat, meta_df$Group)

    # Calculate centroids
    centroids <- pca_df %>%
      group_by(Group) %>%
      summarise(across(c(PC1, PC2), mean, na.rm = TRUE), .groups = 'drop')

    pca_df <- pca_df %>%
      left_join(centroids, by = "Group", suffix = c("", "_centroid"))

    # Define colors (matching violin plot colors)
    group_colors <- c(
      "CANDLE_Cord" = "#FADA7A",
      "CANDLE_Placenta" = "#87CEFA"
    )

    # Add dynamic names for user tissue labels
    group_colors[paste0(user_tissue_label, "_Cord")] <- "#BFAF3B"
    group_colors[paste0(user_tissue_label, "_Placenta")] <- "#1E6091"

    # Filter colors to only available groups
    available_groups <- unique(pca_df$Group)
    group_colors <- group_colors[names(group_colors) %in% available_groups]

    # Prepare annotation text
    global_annot <- paste0("Global PERMANOVA:\nR² = ", r2_val, ", p = ", p_val)

    annot_lines <- c(global_annot)

    if (nrow(pairwise_results) > 0) {
      pairwise_annot_lines <- pairwise_results %>%
        mutate(
          txt = paste0(Group1, " vs ", Group2, ": R²=", R2, ", p=", format.pval(p.adj, digits = 2, eps = .001))
        ) %>%
        pull(txt)

      annot_lines <- c(annot_lines, "\nPairwise PERMANOVA:", pairwise_annot_lines)
    }

    annot_text <- paste(annot_lines, collapse = "\n")

    # Create annotation positioning
    annot_label_df <- data.frame(
      x = max(pca_df$PC1, na.rm = TRUE) * 0.95,
      y = max(pca_df$PC2, na.rm = TRUE) * 0.95,
      label = annot_text
    )

    # Create the plot
    pca_df$Dataset <- paste0("PCA of CANDLE and ", user_tissue_label, " samples using ", length(common_cpgs_matrices), " CpGs")

    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 3, alpha = 0.85) +
      geom_segment(aes(xend = PC1_centroid, yend = PC2_centroid), alpha = 0.3, linetype = "dashed") +
      geom_point(data = centroids, aes(x = PC1, y = PC2, color = Group), size = 6, shape = 4, stroke = 2) +
      labs(
        x = paste0("PC1 (", round(100 * summary(pca_result)$importance[2, 1], 1), "% variance)"),
        y = paste0("PC2 (", round(100 * summary(pca_result)$importance[2, 2], 1), "% variance)")
      ) +
      scale_color_manual(values = group_colors) +
      geom_label(data = annot_label_df, aes(x = x, y = y, label = label),
                 hjust = 1, vjust = 1,
                 fill = NA,
                 label.size = 0.5,
                 color = "black",
                 size = 3.5,
                 inherit.aes = FALSE) +
      facet_wrap(~ Dataset) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_blank(),
        strip.text = element_text(size = 13, hjust = 0.5, color = "black"),
        strip.background = element_rect(fill = "lightgray", color = NA),
        panel.background = element_rect(fill = NA, color = NA),
        plot.background = element_rect(fill = NA, color = NA),
        panel.grid = element_blank(),
        legend.background = element_rect(fill = NA, color = NA),
        legend.key = element_rect(fill = NA, color = NA)
      )


    return(p)
  }

  # Add PCA plot generation
  output$pca_plot <- renderPlot({
    req(input$upload_userdata)

    withProgress(message = "Generating PCA plot...", value = 0, {

      reference_data <- user_plot_methylation_data()
      cord_loaded <- reference_data$cord
      plac_loaded <- reference_data$plac

      # Safely try to read the uploaded CSV file
      user_data <- tryCatch({
        read.csv(input$upload_userdata$datapath, stringsAsFactors = FALSE)
      }, error = function(e) {
        return(NULL)
      })

      incProgress(0.1, detail = "Validating input...")

      # Function to display error in a neutral box, vertically centered
      error_plot <- function(msg) {
        validate(
          need(FALSE, msg)
        )
      }

      # Error Handling
      if (is.null(user_data)) {
        return(error_plot("Note: Failed to read uploaded CSV file."))
      }

      # Check for required columns
      required_cols <- c("TissueType")
      if (!all(required_cols %in% colnames(user_data))) {
        return(error_plot("Note: Uploaded file must contain 'TissueType' column."))
      }

      incProgress(0.2, detail = "Processing data format...")

      # Identify CpG columns (assuming they start with "cg" followed by 8 digits)
      cpg_pattern <- "^cg[0-9]{8}$"
      cpg_columns <- grep(cpg_pattern, colnames(user_data), value = TRUE)

      if (length(cpg_columns) == 0) {
        return(error_plot("Note: CpG not available in COMPUTE."))
      }

      # Validate CpG values - Enhanced validation
      cpg_data <- user_data[, cpg_columns, drop = FALSE]

      # Check for non-numeric values first
      non_numeric_cols <- sapply(cpg_data, function(x) !is.numeric(x))
      if (any(non_numeric_cols)) {
        problem_cols <- names(non_numeric_cols)[non_numeric_cols]
        return(error_plot(paste("Note: Non-numeric values found in CpG columns:", paste(problem_cols[1:min(3, length(problem_cols))], collapse = ", "))))
      }

      # Check for values outside 0-1 range (ignore NAs completely)
      range_violations <- sapply(cpg_data, function(x) {
        valid_values <- x[!is.na(x)]  # Only check non-NA values
        if (length(valid_values) == 0) return(FALSE)  # If all NA, no violation
        any(valid_values < 0 | valid_values > 1)
      })

      if (any(range_violations)) {
        problem_cols <- names(range_violations)[range_violations]
        return(error_plot(paste("Note: Values outside 0-1 range found in CpG columns:", paste(problem_cols[1:min(3, length(problem_cols))], collapse = ", "))))
      }

      # Count only rows with valid CpG data
      # Identify rows that have at least one non-NA CpG value
      valid_cpg_rows <- apply(cpg_data, 1, function(row) {
        any(!is.na(row))
      })

      # Count only rows with CpG data
      num_valid_cpg_samples <- sum(valid_cpg_rows)

      # Filter user_data to only rows with CpG data for tissue validation
      user_data_with_cpg <- user_data[valid_cpg_rows, ]

      # Check if TissueType column has valid entries (only for rows with CpG data)
      tissue_counts <- table(user_data_with_cpg$TissueType)
      if (any(tissue_counts == 0)) {
        return(error_plot("Note: TissueType column contains empty values."))
      }

      # Check minimum data requirements (only for rows with CpG data)
      tissue_sample_counts <- table(user_data_with_cpg$TissueType)
      if (any(tissue_sample_counts < 1)) {
        return(error_plot("Note: Each tissue type must have at least one sample."))
      }

      if (num_valid_cpg_samples == 0) {
        return(error_plot("Note: No samples found with valid CpG data."))
      }

      incProgress(0.3, detail = "Transforming data format...")

      # Transform data: Convert from wide format (CpGs as columns) to long format (CpGs as rows)
      user_data_filtered <- user_data_with_cpg

      # Create sample identifiers if they don't exist
      if (!"Sample_ID" %in% colnames(user_data_filtered)) {
        user_data_filtered$Sample_ID <- paste0("Sample_", seq_len(nrow(user_data_filtered)))
      }


      user_long <- reshape2::melt(user_data_filtered,
                                  id.vars = c("Sample_ID", "TissueType"),
                                  measure.vars = cpg_columns,
                                  variable.name = "CpG_ID",
                                  value.name = "BetaValue")

      # Convert to the expected format for plot functions
      user_transformed <- reshape2::dcast(user_long,
                                          CpG_ID + TissueType ~ Sample_ID,
                                          value.var = "BetaValue")

      # Check if any user CpGs are present in reference tissues - Enhanced validation
      prep_df <- function(df) {
        df <- as.data.frame(df, stringsAsFactors = FALSE)
        if (!"CpG_ID" %in% colnames(df)) {
          df$CpG_ID <- rownames(df)
          rownames(df) <- NULL
        }
        df$CpG_ID <- as.character(df$CpG_ID)
        df
      }

      cord_df <- prep_df(cord_loaded)
      plac_df <- prep_df(plac_loaded)

      cpgs_in_user <- unique(user_transformed$CpG_ID)
      cpgs_in_cord <- unique(cord_df$CpG_ID)
      cpgs_in_plac <- unique(plac_df$CpG_ID)
      cpgs_in_ref <- union(cpgs_in_cord, cpgs_in_plac)

      # Identify which CpGs are missing from reference data
      missing_cpgs <- setdiff(cpgs_in_user, cpgs_in_ref)
      common_cpgs <- intersect(cpgs_in_user, cpgs_in_ref)

      if (length(common_cpgs) == 0) {
        if (length(missing_cpgs) <= 5) {
          return(error_plot(paste("Note: None of the uploaded CpG IDs are available in COMPUTE reference data. Missing CpGs:", paste(missing_cpgs, collapse = ", "))))
        } else {
          return(error_plot(paste("Note: None of the uploaded CpG IDs are available in COMPUTE reference data. Missing", length(missing_cpgs), "CpGs including:", paste(missing_cpgs[1:3], collapse = ", "), "...")))
        }
      }

      # Handle partially missing CpGs - provide user feedback
      if (length(missing_cpgs) > 0) {
        warning_msg <- if (length(missing_cpgs) <= 5) {
          paste("Note:", length(missing_cpgs), "CpG(s) not found in reference data and will be skipped:", paste(missing_cpgs, collapse = ", "))
        } else {
          paste("Note:", length(missing_cpgs), "CpG(s) not found in reference data and will be skipped. First few:", paste(missing_cpgs[1:3], collapse = ", "), "...")
        }

        # For console logging
        cat("Warning:", warning_msg, "\n")
      }

      incProgress(0.5, detail = "Preparing PCA data...")

      user_cord <- user_transformed[user_transformed$TissueType == "Cord Blood", ]
      user_plac <- user_transformed[user_transformed$TissueType == "Placenta", ]

      if (nrow(user_cord) > 0 && "TissueType" %in% colnames(user_cord)) {
        user_cord$TissueType <- NULL
      }
      if (nrow(user_plac) > 0 && "TissueType" %in% colnames(user_plac)) {
        user_plac$TissueType <- NULL
      }

      has_user_cord <- nrow(user_cord) > 0
      has_user_plac <- nrow(user_plac) > 0

      cat("User Cord samples:", if(has_user_cord) ncol(user_cord) - 1 else 0, "\n")
      cat("User Placenta samples:", if(has_user_plac) ncol(user_plac) - 1 else 0, "\n")

      if (!has_user_cord && !has_user_plac) {
        return(error_plot("No valid tissue types found in user data."))
      }

      # Filter to common CpGs only
      if (has_user_cord) {
        user_cord <- user_cord[user_cord$CpG_ID %in% common_cpgs, ]
      }
      if (has_user_plac) {
        user_plac <- user_plac[user_plac$CpG_ID %in% common_cpgs, ]
      }

      # Combine user data back
      user_final <- rbind(user_cord, user_plac)

      if (has_user_cord && nrow(user_cord) > 0) {
        user_cord$TissueType <- "Cord Blood"
      }
      if (has_user_plac && nrow(user_plac) > 0) {
        user_plac$TissueType <- "Placenta"
      }

      user_final <- rbind(
        if(has_user_cord && nrow(user_cord) > 0) user_cord else NULL,
        if(has_user_plac && nrow(user_plac) > 0) user_plac else NULL
      )

      if (length(common_cpgs) < 3) {
        return(error_plot(paste("Insufficient CpGs for PCA analysis.\nNeed at least 3 CpGs, found:", length(common_cpgs))))
      }

      incProgress(0.7, detail = "Running PCA analysis...")

      # Generate PCA plot using the same transformed data
      pca_plot <- plot_pca_comparison(common_cpgs, cord_df, plac_df, user_final, user_tissue_label = "User")

      incProgress(1.0, detail = "Complete!")

      return(pca_plot)

    })
  }, bg = "transparent")

  ##################### Epigenetic Clocks

  CANDLE_clocks$GestAge_numeric <- as.numeric(as.character(CANDLE_clocks$GestAge))
  CANDLE_celltypes <- readRDS("Data/CANDLE_celltypes.rds")

timing <- reactiveValues(start_time = NULL)

observeEvent(input$upload_userdata, {
  timing$start_time <- Sys.time()
  # Send start time to JavaScript
  session$sendCustomMessage("startTiming", as.character(timing$start_time))
})

# Receive the timing back from JavaScript
observeEvent(input$total_timing, {
  cat("FULL END-TO-END TIME (JavaScript measured):", input$total_timing, "seconds\n")
})

  # Basic validation
  user_combined_data <- reactive({
    req(input$upload_userdata)

    if (data_security$data_cleared) {
      return(NULL)
    }

   data_security$last_activity <- Sys.time()


    user_data <- tryCatch({
      read.csv(input$upload_userdata$datapath)
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(user_data)) {
      return(NULL)
    }

    # Check TissueType
    if (!"TissueType" %in% colnames(user_data)) {
      return(NULL)
    }

    # Check valid TissueType values
    if (!all(user_data$TissueType %in% c("Cord Blood", "Placenta"))) {
      return(NULL)
    }

    user_data
  })

  # ========= helper functions =========
  prepare_violin_data <- function(data, pattern, source) {
    filtered_data <- data[, grep(pattern, colnames(data), value = TRUE), drop = FALSE]
    filtered_data <- filtered_data[, !grepl("sample_id|Chip_Position", colnames(filtered_data))]
    long_data <- reshape2::melt(
      filtered_data,
      variable.name = "Cell_Type",
      value.name = "Cell_Proportions"
    )
    long_data$Cell_Type <- gsub("_cord|_plac", "", long_data$Cell_Type)
    long_data$Source <- source
    return(long_data)
  }

  merge_data_for_plotting <- function(candle_data, user_data, pattern, source_name) {
    user_long <- prepare_violin_data(user_data, pattern, source = source_name)
    matching_cell_types <- intersect(candle_data$Cell_Type, user_long$Cell_Type)
    candle_subset <- candle_data[candle_data$Cell_Type %in% matching_cell_types, ]
    user_subset <- user_long[user_long$Cell_Type %in% matching_cell_types, ]
    missing_cols_in_user <- setdiff(names(candle_subset), names(user_subset))
    missing_cols_in_candle <- setdiff(names(user_subset), names(candle_subset))
    for (col in missing_cols_in_user) user_subset[[col]] <- NA
    for (col in missing_cols_in_candle) candle_subset[[col]] <- NA
    user_subset <- user_subset[, names(candle_subset)]
    combined_data <- rbind(candle_subset, user_subset)
    return(combined_data)
  }

  render_violin_plot <- function(data, title, x_label, y_label, legend_labels = NULL, missing_note = NULL) {
    data$Source_Tissue <- paste(data$Source, data$TissueType, sep = "_")

    all_colors <- c(
      "CANDLE_Cord Blood" = "#FADA7A",
      "CANDLE_Placenta"  = "#87cefa",
      "User_Cord Blood"  = "#bfaf3b",
      "User_Placenta"    = "#1E6091"
    )

    used_levels <- unique(data$Source_Tissue)
    used_colors <- all_colors[used_levels]

    # Use legend_labels as passed
    if (is.null(legend_labels)) {
      legend_labels <- names(used_colors)
    }

    # Base plot
    base_plot <- ggplot(data, aes(x = Cell_Type, y = Cell_Proportions, fill = Source_Tissue)) +
      geom_violin(trim = FALSE, alpha = 0.7, position = position_dodge(width = 0.8)) +
      geom_jitter(aes(color = Source_Tissue), alpha = 0.3, size = 0.8,
                  position = position_jitterdodge(dodge.width = 0.8)) +
      scale_fill_manual(
        values = used_colors,
        breaks = names(used_colors),
        labels = legend_labels,
        name = "Source"
      ) +
      scale_color_manual(
        values = used_colors,
        breaks = names(used_colors),
        labels = legend_labels,
        name = "Source"
      ) +
      ylim(0, 1) +
      ylab(y_label) +
      xlab(x_label) +
      ggtitle(title) +
      theme_minimal() +
      theme(
        text = element_text(size = 12),
        panel.spacing = unit(0.5, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0)
      )

    # Add missing note as caption if present
    if (!is.null(missing_note) && nchar(missing_note) > 0) {
      base_plot <- base_plot + labs(caption = missing_note) +
        theme(plot.caption = element_text(hjust = 0.5, color = "#6c757d", face = "italic", size = 10))
    }

    return(base_plot)
  }

  # ========= 1) Render scatter_both (clocks) =========

  compute_stats <- function(df) {
    r <- cor(df$GestAge, df$EpigenAge, method = "pearson", use = "complete.obs")
    p_val <- cor.test(df$GestAge, df$EpigenAge)$p.value
    MAE <- median(abs(df$GestAge - df$EpigenAge), na.rm = TRUE)
    label <- paste0("r = ", round(r, 3),
                    ", p = ", format.pval(p_val, digits = 2, eps = .001),
                    ", MdAE = ", round(MAE, 3), " weeks")
    return(label)
  }

  output$scatter_both <- renderPlot({

    # Check if user data exists first
    user_age_data <- user_combined_data()


    error_plot <- function(msg) {
      ggplot() +
        xlim(0, 1) +
        ylim(0, 1) +
        geom_text(aes(x = 0.05, y = 0.95),
                  label = msg,
                  size = 4,
                  color = "gray60",
                  hjust = 0,
                  vjust = 1) +
        theme_void() +
        theme(
          plot.background = element_rect(fill = NA, color = NA),
          panel.background = element_rect(fill = NA, color = NA),
          plot.margin = margin(20, 20, 20, 20)
        )
    }


    if (is.null(user_age_data)) {
      return(error_plot("Note: Please upload a valid CSV file with required columns:\nGestAge, EpigeneticAge, TissueType"))
    }

    required_cols <- c("GestAge", "EpigeneticAge", "TissueType")
    missing_cols <- setdiff(required_cols, colnames(user_age_data))
    if (length(missing_cols) > 0) {
      return(error_plot(paste("Note: Missing required columns:", paste(missing_cols, collapse = ", "))))
    }

    valid_tissues <- c("Cord Blood", "Placenta")
    invalid_tissues <- setdiff(unique(user_age_data$TissueType), valid_tissues)
    if (length(invalid_tissues) > 0) {
      return(error_plot(paste("Note: Invalid TissueType values found:", paste(invalid_tissues, collapse = ", "),
                              "\nAllowed values: Cord Blood, Placenta")))
    }

    user_age_data$GestAge_numeric <- suppressWarnings(as.numeric(as.character(user_age_data$GestAge)))
    user_age_data$EpigeneticAge_numeric <- suppressWarnings(as.numeric(as.character(user_age_data$EpigeneticAge)))

    if (all(is.na(user_age_data$GestAge_numeric))) {
      return(error_plot("Note: GestAge column contains no valid numeric values."))
    }

    if (all(is.na(user_age_data$EpigeneticAge_numeric))) {
      return(error_plot("Note: EpigeneticAge column contains no valid numeric values."))
    }

    user_age_data <- user_age_data[!is.na(user_age_data$GestAge_numeric) &
                                     !is.na(user_age_data$EpigeneticAge_numeric), ]
    if (nrow(user_age_data) == 0) {
      return(error_plot("Note: No valid data rows after removing missing values."))
    }


    # Subset data by tissue type
    user_cord <- subset(user_age_data, TissueType == "Cord Blood")
    user_placenta <- subset(user_age_data, TissueType == "Placenta")

    # Calculate axis limits
    min_x_y <- min(
      min(CANDLE_clocks$GestAge_numeric, na.rm = TRUE),
      min(CANDLE_clocks$DNAme_GA_Bohlin, na.rm = TRUE),
      min(CANDLE_clocks$DNAme_GA_RPC, na.rm = TRUE),
      min(user_age_data$GestAge_numeric, na.rm = TRUE),
      min(user_age_data$EpigeneticAge_numeric, na.rm = TRUE)
    )

    max_x_y <- max(
      max(CANDLE_clocks$GestAge_numeric, na.rm = TRUE),
      max(CANDLE_clocks$DNAme_GA_Bohlin, na.rm = TRUE),
      max(CANDLE_clocks$DNAme_GA_RPC, na.rm = TRUE),
      max(user_age_data$GestAge_numeric, na.rm = TRUE),
      max(user_age_data$EpigeneticAge_numeric, na.rm = TRUE)
    )

    axis_breaks <- seq(floor(min_x_y), ceiling(max_x_y), by = 1)
    padding <- 0.5
    x_limits <- c(floor(min_x_y) - padding, ceiling(max_x_y) + padding)
    y_limits <- c(floor(min_x_y) - padding, ceiling(max_x_y) + padding)

    # Define theme
    clean_theme <- theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        text = element_text(size = 13),
        axis.line = element_line(color = "black", size = 0.6),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5)
      )

    # Calculate sample sizes for legend
    candle_cord_n <- sum(!is.na(CANDLE_clocks$DNAme_GA_Bohlin))
    user_cord_n <- nrow(user_cord)
    candle_placenta_n <- sum(!is.na(CANDLE_clocks$DNAme_GA_RPC))
    user_placenta_n <- nrow(user_placenta)

    # Define colors
    candle_cord_color <- "#FADA7A"     # Amber
    candle_placenta_color <- "#87cefa"  # Light blue
    user_cord_clock <- "#bfaf3b"       # Dark amber
    user_placenta_clock <- "#1E6091"   # Dark blue

    # Create cord blood plot
    if (user_cord_n > 0) {
      cord_data <- rbind(
        data.frame(
          Source = "CANDLE",
          GestAge = CANDLE_clocks$GestAge_numeric,
          EpigenAge = CANDLE_clocks$DNAme_GA_Bohlin
        ),
        data.frame(
          Source = "User",
          GestAge = user_cord$GestAge_numeric,
          EpigenAge = user_cord$EpigeneticAge_numeric
        )
      )
      cord_data <- subset(cord_data, !is.na(EpigenAge))

      cord_plot <- ggplot(cord_data, aes(x = GestAge, y = EpigenAge, color = Source)) +
        geom_hline(yintercept = axis_breaks, color = "gray90", size = 0.4) +
        geom_vline(xintercept = axis_breaks, color = "gray90", size = 0.4) +
        geom_point(size = 3, alpha = 0.6) +
        geom_smooth(aes(color = Source), method = "lm", se = FALSE, size = 1) +  # regression line
        labs(title = "Cord Blood",
             x = "Gestational Age (Weeks)",
             y = "Epigenetic Gestational Age (Weeks)") +
        scale_x_continuous(breaks = axis_breaks, expand = c(0, 0)) +
        scale_y_continuous(breaks = axis_breaks, expand = c(0, 0)) +
        geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
        coord_cartesian(xlim = x_limits, ylim = y_limits) +
        scale_color_manual(
          values = c("CANDLE" = candle_cord_color, "User" = user_cord_clock),
          labels = c(
            CANDLE = paste0("CANDLE (n = ", candle_cord_n, ")"),
            User = paste0("User (n = ", user_cord_n, ")")
          ),
          name = "Source"
        ) +
        clean_theme +
        theme(legend.position = "right")
    } else {
      cord_plot <- error_plot("Note: No Cord Blood samples found in your data")
    }

    # Create placenta plot
    if (user_placenta_n > 0) {
      placenta_data <- rbind(
        data.frame(
          Source = "CANDLE",
          GestAge = CANDLE_clocks$GestAge_numeric,
          EpigenAge = CANDLE_clocks$DNAme_GA_RPC
        ),
        data.frame(
          Source = "User",
          GestAge = user_placenta$GestAge_numeric,
          EpigenAge = user_placenta$EpigeneticAge_numeric
        )
      )
      placenta_data <- subset(placenta_data, !is.na(EpigenAge))

      placenta_plot <- ggplot(placenta_data, aes(x = GestAge, y = EpigenAge, color = Source)) +
        geom_hline(yintercept = axis_breaks, color = "gray90", size = 0.4) +
        geom_vline(xintercept = axis_breaks, color = "gray90", size = 0.4) +
        geom_point(size = 3, alpha = 0.6) +
        geom_smooth(aes(color = Source), method = "lm", se = FALSE, size = 1) +  # regression line
        labs(title = "Placenta",
             x = "Gestational Age (Weeks)",
             y = "Epigenetic Gestational Age (Weeks)") +
        scale_x_continuous(breaks = axis_breaks, expand = c(0, 0)) +
        scale_y_continuous(breaks = axis_breaks, expand = c(0, 0)) +
        geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
        coord_cartesian(xlim = x_limits, ylim = y_limits) +
        scale_color_manual(
          values = c("CANDLE" = candle_placenta_color, "User" = user_placenta_clock),
          labels = c(
            CANDLE = paste0("CANDLE (n = ", candle_placenta_n, ")"),
            User = paste0("User (n = ", user_placenta_n, ")")
          ),
          name = "Source"
        ) +
        clean_theme +
        theme(legend.position = "right")
    } else {
      placenta_plot <- error_plot("Note: No Placenta samples found in your data")
    }

    # For cord blood plot stats:
    if (user_cord_n > 0) {
      cord_candle_label <- compute_stats(subset(cord_data, Source == "CANDLE"))
      cord_user_label   <- compute_stats(subset(cord_data, Source == "User"))

      # Position text in top-left corner with better spacing
      text_x_pos <- x_limits[1] + 0.2
      text_y_top <- y_limits[2] - 0.3
      text_y_spacing <- 0.8

      cord_plot <- cord_plot +
        annotate("text",
                 x = text_x_pos,
                 y = text_y_top,
                 label = paste0("CANDLE: ", cord_candle_label),
                 hjust = 0, vjust = 1, size = 3.5,
                 color = candle_cord_color,
                 fontface = "bold") +
        annotate("text",
                 x = text_x_pos,
                 y = text_y_top - text_y_spacing,
                 label = paste0("User: ", cord_user_label),
                 hjust = 0, vjust = 1, size = 3.5,
                 color = user_cord_clock,
                 fontface = "bold")
    }

    # For placenta plot stats:
    if (user_placenta_n > 0) {
      placenta_candle_label <- compute_stats(subset(placenta_data, Source == "CANDLE"))
      placenta_user_label   <- compute_stats(subset(placenta_data, Source == "User"))

      # Position text in top-left corner with better spacing
      text_x_pos <- x_limits[1] + 0.2
      text_y_top <- y_limits[2] - 0.3
      text_y_spacing <- 0.8

      placenta_plot <- placenta_plot +
        annotate("text",
                 x = text_x_pos,
                 y = text_y_top,
                 label = paste0("CANDLE: ", placenta_candle_label),
                 hjust = 0, vjust = 1, size = 3.5,
                 color = candle_placenta_color,
                 fontface = "bold") +
        annotate("text",
                 x = text_x_pos,
                 y = text_y_top - text_y_spacing,
                 label = paste0("User: ", placenta_user_label),
                 hjust = 0, vjust = 1, size = 3.5,
                 color = user_placenta_clock,
                 fontface = "bold")
    }



    # Create spacer
    spacer <- grid::nullGrob()

    # Arrange plots
    tryCatch({
      gridExtra::grid.arrange(
        cord_plot, spacer, placenta_plot,
        nrow = 3,
        heights = c(1, 0.1, 1)
      )
    }, error = function(e) {
      if (user_cord_n > 0) {
        return(cord_plot)
      } else if (user_placenta_n > 0) {
        return(placenta_plot)
      } else {
        return(ggplot() +
                 geom_text(aes(x = 0, y = 0), label = "Error arranging plots") +
                 theme_void())
      }
    })

  }, bg = "transparent")


  # ========= 2) Render Violin Plots =========
  check_missing_cols <- function(data, required_cols) {
    missing <- required_cols[
      !required_cols %in% colnames(data) |
        sapply(required_cols, function(c) {
          c %in% colnames(data) && all(is.na(data[[c]]))
        })
      ]
    return(missing)
  }

  # ---- Cord Blood violin ----
  output$cord_violin_plot <- renderPlot({
    req(user_combined_data())
    user_data <- user_combined_data()

    error_plot <- function(msg) {
      validate(
        need(FALSE, msg)
      )
    }


    # Check for samples
    user_cord_samples_n <- nrow(subset(user_data, TissueType == "Cord Blood"))
    if (user_cord_samples_n == 0) {
      return(error_plot("Note: No Cord Blood samples found in your data"))
    }

    required_cord_cols <- c(
      "CD4T_cord", "CD8T_cord", "Mono_cord", "Gran_cord",
      "NK_cord", "nRBC_cord", "Bcell_cord"
    )
    missing_cols <- check_missing_cols(user_data, required_cord_cols)

    # If ALL missing, throw error
    if (length(missing_cols) == length(required_cord_cols)) {
      return(error_plot("Note: Missing required cell type columns for Cord Blood plot."))
    }

    missing_note <- NULL
    if (length(missing_cols) > 0) {
      missing_note <- paste(
        "Note: The following Cord Blood cell types are missing from your data:",
        paste(gsub("_cord$", "", missing_cols), collapse = ", ")
      )
    }



    # Merge with CANDLE reference
    combined_cord_data <- merge_data_for_plotting(
      candle_data = prepare_violin_data(CANDLE_celltypes, "_cord$", "CANDLE"),
      user_data = user_data,
      pattern = "_cord$",
      source_name = "User"
    )
    combined_cord_data$TissueType <- "Cord Blood"

    # Legend labels
    candle_sample_n <- 551
    legend_labels <- c(
      "CANDLE_Cord Blood" = paste0("CANDLE (n = ", candle_sample_n, ")"),
      "User_Cord Blood" = paste0("User (n = ", user_cord_samples_n, ")")
    )

    render_violin_plot(
      combined_cord_data,
      title = "Cord Blood",
      x_label = "",
      y_label = "Estimated Cell Proportions",
      legend_labels = legend_labels,
      missing_note = missing_note
    )
  }, bg = "transparent")


  # ---- Placenta violin ----
  output$plac_violin_plot <- renderPlot({
    req(user_combined_data())
    user_data <- user_combined_data()

    error_plot <- function(msg) {
      validate(
        need(FALSE, msg)
      )
    }


    # Check for samples
    user_placenta_samples_n <- nrow(subset(user_data, TissueType == "Placenta"))
    if (user_placenta_samples_n == 0) {
      return(error_plot("Note: No Placenta samples found in your data"))
    }

    # Required columns
    required_plac_cols <- c(
      "Trophoblasts_plac", "Stromal_plac", "Hofbauer_plac",
      "Endothelial_plac", "nRBC_plac", "Syncytiotrophoblast_plac"
    )
    missing_cols <- check_missing_cols(user_data, required_plac_cols)

    # If ALL missing, throw error
    if (length(missing_cols) == length(required_plac_cols)) {
      return(error_plot("Note: Missing required cell type columns for Placenta plot."))
    }

    # Otherwise, make a note if partially missing
    missing_note <- NULL
    if (length(missing_cols) > 0) {
      missing_note <- paste(
        "Note: The following Placenta cell types are missing from your data:",
        paste(gsub("_plac$", "", missing_cols), collapse = ", ")
      )
    }


    # Merge with CANDLE reference
    combined_plac_data <- merge_data_for_plotting(
      candle_data = prepare_violin_data(CANDLE_celltypes, "_plac$", "CANDLE"),
      user_data = user_data,
      pattern = "_plac$",
      source_name = "User"
    )
    combined_plac_data$TissueType <- "Placenta"

    # Legend labels
    candle_sample_n <- 551
    legend_labels <- c(
      "CANDLE_Placenta" = paste0("CANDLE (n = ", candle_sample_n, ")"),
      "User_Placenta" = paste0("User (n = ", user_placenta_samples_n, ")")
    )

    render_violin_plot(
      combined_plac_data,
      title = "Placenta",
      x_label = "",
      y_label = "Estimated Cell Proportions",
      legend_labels = legend_labels,
      missing_note = missing_note
    )
  }, bg = "transparent")


}

