library(tictoc)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggExtra)
library(ggpubr)
library(ggtext)
library(cowplot)
library(qs)

call_anno <- function() {
  library(qs)

  CANDLE_corr_anno <- readRDS("Data/CANDLE_corr_anno.rds")
  anno <- CANDLE_corr_anno %>% dplyr::select(CpG,Chr,MapInfo,Gene,Gene_Region,CpG_Island_Class,
                                             Cross_Hybridizing,Correlation, Is_Variable,
                                             Is_in_CMR,CMR_Cord,CMR_Length_Cord,CMR_Placenta,CMR_Length_Placenta)

  return(list("anno" = anno))
}

# Load methylation data from chunks
get_methylation_data_chunked <- function(cpgs, anno) {
  library(qs)

  # Find which chromosomes these CpGs are on
  needed_chrs <- unique(anno$Chr[anno$CpG %in% cpgs])

  if (length(needed_chrs) == 0) {
    return(list(cord = NULL, plac = NULL))
  }

  # Initialize results
  cord_results <- list()
  plac_results <- list()

  # Process one chromosome at a time
  for (chr in needed_chrs) {
    # Check if chunk file exists
    cord_file <- paste0("Data/chunks/cord_chr", chr, ".qs")
    plac_file <- paste0("Data/chunks/plac_chr", chr, ".qs")

    if (!file.exists(cord_file) || !file.exists(plac_file)) {
      warning(paste("Chunk files not found for chromosome", chr))
      next
    }

    # ============ CORD PROCESSING ============
    cord_chunk <- qread(cord_file)

    cpgs_in_chunk <- intersect(cpgs, rownames(cord_chunk))

    if (length(cpgs_in_chunk) > 0) {
      cord_extract <- cord_chunk[cpgs_in_chunk, , drop = FALSE]

      # Save metadata BEFORE any conversion
      saved_rownames <- rownames(cord_extract)
      saved_colnames <- colnames(cord_extract)
      n_rows <- nrow(cord_extract)
      n_cols <- ncol(cord_extract)

      # Convert from float32
      if (inherits(cord_extract, "float32")) {
        if (requireNamespace("float", quietly = TRUE)) {
          cord_extract <- as.numeric(float::dbl(cord_extract))
          cord_extract <- matrix(cord_extract, nrow = n_rows, ncol = n_cols)
          rownames(cord_extract) <- saved_rownames
          colnames(cord_extract) <- saved_colnames
        }
      } else if (!is.matrix(cord_extract)) {
        cord_extract <- as.matrix(cord_extract)
      }

      cord_results[[chr]] <- cord_extract
      rm(cord_extract)
    }

    rm(cord_chunk, cpgs_in_chunk, saved_rownames, saved_colnames, n_rows, n_cols)
    gc(full = TRUE, reset = TRUE)  # Full GC with heap reset

    # Small delay to ensure GC completes
    Sys.sleep(0.05)

    # ============ PLAC PROCESSING ============
    plac_chunk <- qread(plac_file)

    cpgs_in_chunk <- intersect(cpgs, rownames(plac_chunk))

    if (length(cpgs_in_chunk) > 0) {
      plac_extract <- plac_chunk[cpgs_in_chunk, , drop = FALSE]

      # Save metadata BEFORE any conversion
      saved_rownames <- rownames(plac_extract)
      saved_colnames <- colnames(plac_extract)
      n_rows <- nrow(plac_extract)
      n_cols <- ncol(plac_extract)

      # Convert from float32
      if (inherits(plac_extract, "float32")) {
        if (requireNamespace("float", quietly = TRUE)) {
          plac_extract <- as.numeric(float::dbl(plac_extract))
          plac_extract <- matrix(plac_extract, nrow = n_rows, ncol = n_cols)
          rownames(plac_extract) <- saved_rownames
          colnames(plac_extract) <- saved_colnames
        }
      } else if (!is.matrix(plac_extract)) {
        plac_extract <- as.matrix(plac_extract)
      }

      plac_results[[chr]] <- plac_extract
      rm(plac_extract)
    }

    rm(plac_chunk, cpgs_in_chunk, saved_rownames, saved_colnames, n_rows, n_cols)
    gc(full = TRUE, reset = TRUE)

    Sys.sleep(0.05)
  }

  # Combine results
  cord <- if (length(cord_results) > 0) do.call(rbind, cord_results) else NULL
  plac <- if (length(plac_results) > 0) do.call(rbind, plac_results) else NULL

  # Final cleanup of intermediate results
  rm(cord_results, plac_results)
  gc(full = TRUE)

  return(list(cord = cord, plac = plac))
}
#################################################################

data_prep <- function(tmp_df, cord, plac) {
  # Subset
  cord_subset <- cord[rownames(cord) %in% tmp_df, , drop = FALSE]
  plac_subset <- plac[rownames(plac) %in% tmp_df, , drop = FALSE]

  # Transpose and prepare
  cord_long <- as.data.frame(t(cord_subset))
  plac_long <- as.data.frame(t(plac_subset))
  cord_long$Subject <- rownames(cord_long)
  plac_long$Subject <- rownames(plac_long)

  merged <- merge(cord_long, plac_long, by = "Subject", suffixes = c("_cord", "_plac"))

  long_data <- merged %>%
    pivot_longer(cols = -Subject,
                 names_to = c("CpG", "Tissue"),
                 names_sep = "_",
                 values_to = "Beta") %>%
    pivot_wider(names_from = Tissue, values_from = Beta) %>%
    drop_na(cord, plac)

  num_cpgs <- length(unique(long_data$CpG))
  plot_color <- "#4D4D4D"

  clean_theme <- theme_minimal(base_size = 15) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(color = "black", size = 0.6),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      axis.text.x = element_text(angle = 0),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14),
      strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "italic", size = 15)
    )

  base_plot <- ggplot(long_data, aes(x = cord, y = plac)) +
    geom_hline(yintercept = seq(0, 1, by = 0.1), color = "gray90", size = 0.4) +
    geom_vline(xintercept = seq(0, 1, by = 0.1), color = "gray90", size = 0.4) +
    geom_point(color = plot_color, alpha = 0.4, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    xlim(0, 1) + ylim(0, 1) +
    xlab("Cord Blood Beta Value") + ylab("Placenta Beta Value") +
    clean_theme

  # Correlation stats
  cor_result <- cor.test(long_data$cord, long_data$plac, method = "pearson")
  r_val <- cor_result$estimate
  r_fmt <- formatC(r_val, format = "f", digits = 3)
  r_squared <- r_val^2
  r2_fmt <- formatC(r_squared, format = "f", digits = 3)
  p_val <- cor_result$p.value
  p_fmt <- if (p_val < 0.001) format(p_val, scientific = TRUE, digits = 2) else round(p_val, 3)

  # Add correlation annotation to plot
  base_plot <- base_plot +
    annotate("text", x = 0.05, y = 0.95,
             label = paste0("ρ = ", r_fmt, "\nR² = ", r2_fmt, "\np = ", p_fmt),
             hjust = 0, vjust = 1, size = 5, color = "black")

  # Create marginal density plot
  marginal_plot <- ggMarginal(
    base_plot,
    type = "density",
    margins = "both",
    size = 5,
    alpha = 0.4,
    fill = plot_color,
    color = plot_color
  )

  # Return the marginal plot without any title (title is now handled in the UI)
  return(marginal_plot)
}


get_correlation_stats <- function(tmp_df, cord, plac) {
  # Subset
  cord_subset <- cord[rownames(cord) %in% tmp_df, , drop = FALSE]
  plac_subset <- plac[rownames(plac) %in% tmp_df, , drop = FALSE]

  cord_long <- as.data.frame(t(cord_subset))
  plac_long <- as.data.frame(t(plac_subset))
  cord_long$Subject <- rownames(cord_long)
  plac_long$Subject <- rownames(plac_long)

  merged <- merge(cord_long, plac_long, by = "Subject", suffixes = c("_cord", "_plac"))

  long_data <- merged %>%
    pivot_longer(cols = -Subject,
                 names_to = c("CpG", "Tissue"),
                 names_sep = "_",
                 values_to = "Beta") %>%
    pivot_wider(names_from = Tissue, values_from = Beta) %>%
    drop_na(cord, plac)

  long_data %>%
    group_by(CpG) %>%
    summarise(
      r = cor(cord, plac, method = "pearson"),
      r_squared = r^2,
      p_value = cor.test(cord, plac, method = "pearson")$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      Correlation = paste0(
        "ρ = ", formatC(r, digits = 3, format = "f"), "\n",
        "R² = ", formatC(r_squared, digits = 3, format = "f"), "\n",
        "p = ", ifelse(p_value < 0.001,
                       format(p_value, scientific = TRUE, digits = 2),
                       round(p_value, 3))
      )
    ) %>%
    select(CpG, Correlation)
}

#############################################################

cmr_prep <- function(anno, forplot, colour, tissue) {
  # Merge with genomic annotations
  forplot <- merge(forplot, anno[, c("CpG", "MapInfo", "Chr")],
                   by.x = "CpG", by.y = "CpG")

  forplot <- reshape2::melt(forplot, id = c("CpG", "MapInfo", "Chr"))
  colnames(forplot) <- c("CpG", "Genomic_Location", "Chr", "Subject", "Beta")

  forplot$Tissue <- tissue

  # Order by genomic location
  forplot <- forplot[order(forplot$Genomic_Location), ]

  # Title for facet
  forplot$title <- paste0(
    "Chr ", unique(forplot$Chr), ": ",
    min(forplot$Genomic_Location), "-", max(forplot$Genomic_Location)
  )

  # Plot
  plot <- ggplot(forplot, aes(x = Genomic_Location, y = Beta, color = Tissue)) +
    geom_point() +
    stat_summary(fun = "mean", aes(group = Tissue), geom = "line", size = 1, show.legend = FALSE) +
    scale_color_manual(values = setNames(colour, tissue), name = "Tissue") +
    theme_bw(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 12)
    ) +
    xlab("Genomic Location") +
    ylab("Beta Value") +
    ylim(0, 1) +
    facet_grid(. ~ title)


  return(plot)
}




cmr_prep_common <- function(anno, forplot) {
  # Merge plot data with genomic location annotations
  forplot <- merge(forplot, anno[, c("CpG", "MapInfo", "Chr")],
                   by.x = "CpG", by.y = "CpG")
  forplot <- reshape2::melt(forplot, id = c("CpG", "MapInfo", "Chr", "Tissue"))
  colnames(forplot) <- c("CpG", "Genomic_Location", "Chr", "Tissue", "Subject", "Beta")

  # Order by genomic location
  forplot <- forplot[order(forplot$Genomic_Location), ]

  forplot$Tissue <- factor(
    forplot$Tissue,
    levels = c("CordBlood", "Placenta"),
    labels = c("Cord Blood", "Placenta")
  )

  forplot$title <- paste0(
    "Chr ", unique(forplot$Chr), ": ",
    min(forplot$Genomic_Location), "-", max(forplot$Genomic_Location)
  )

  # Plot with genomic location on x-axis
  plot2 <- ggplot(forplot, aes(x = Genomic_Location, y = Beta, color = Tissue)) +
    geom_point() +
    stat_summary(fun = "mean", aes(group = Tissue), geom = "line", size = 1, show.legend = FALSE) +
    scale_color_manual(values = c(
      "Cord Blood" = '#FADA7A',
      "Placenta" = '#87cefa'
    )) +
    theme_bw(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    xlab("Genomic Location") +
    labs(y = "Beta Value") +
    ylim(0, 1) +
    facet_grid(. ~ title)

  return(plot2)
}
