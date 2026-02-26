library(shiny)
library(shinycssloaders)
library(ComplexHeatmap)
library(tidyverse)

# ==== define ui ====

ui <- fluidPage(
  titlePanel("Newton/Giembycz Lab RNAseq Explorer"),
  sidebarLayout(
    sidebarPanel(
      
      selectInput(
        "dataset", "Dataset:",
        choices = list(
          "A549 datasets" = c(
            "A549 IL1B Bud",
            "A549 IL1B Bud (transcripts)",
            "A549 Ad-DUSP1 IL1B",
            "A549 Ad-IkBa IL1B Dex"
          ),
          "NEW/EXPERIMENTAL" = c(
            "Basal expression - single gene",
            "Basal expression - heatmap"
          ),
          "ALI datasets" = c(
            "ALI IL1B Bud",
            "ALI Bud Form",
            "ALI Cig IL17A Bud",
            "ALI CMV",
            "ALI HRV (Bai 2015)"
          ),
          "BEAS-2B datasets" = c(
            "BEAS-2B AKO/BKO Form Bud",
            "BEAS-2B AKO/BKO Form Bud (transcripts)",
            "BEAS-2B dKO Form",
            "BEAS-2B dKO Form (transcripts)",
            "BEAS-2B Mom Ind",
            "BEAS-2B Tap Bay Vil",
            "BEAS-2B RNO ONO Vil (GSE267218)"
          ),
          "HBE datasets" = c(
            "HBE IL1B IFNg Dex",
            "HBE TNF Form",
            "HBE CMV timecourse (Parkins Lab)"
          )
        )
      ),
      
      # textInput("gene", "Official gene symbol:", value = "FKBP5"),
      uiOutput("gene_input"),
      
      radioButtons(
        "metric",
        "Y-axis:",
        choices = c("log2_tpm", "log2_fold", "tpm", "fold"),
        selected = "log2_tpm"
      ),
      
      actionButton("button", "Generate", icon = icon("redo")),
      
      helpText(HTML("
      v2.2, February 2025<br/>
      - Added experimental heatmap option for basal expression data<br/>
      <br/>
      v2.1, December 2025<br/>
      - Plots now update only when Generate is clicked<br/>
      - Added transcript-level datasets<br/>
      - Added log2fold, tpm, and fold options<br/>
      - Added loading indicator<br/>
      - Datasets cleaned and consolidated (https://tinyurl.com/2ae5ajap)<br/>
      - Large datasets chunked to avoid OOM crash"
      )
      )
      
    ),
    mainPanel(
      tabsetPanel(type = "tabs", 
                  tabPanel("Plot", plotOutput("expression_plot") %>% 
                             withSpinner(
                               type = 6, 
                               color = "#555555", 
                               caption = "Working...")), 
                  tabPanel("Table", tableOutput("table"))
      )
    )
  )
)

single_gene_ui <- textInput(
  "gene",
  "Official gene symbol:",
  value = "FKBP5"
)

multi_gene_ui <- textAreaInput(
  "genes",
  "Gene list (one per line or comma-separated):",
  value = "FKBP5\nDUSP1\nNFKBIA",
  rows = 6
)

# ==== define server ====
server <- function(input, output) {
  
  ## ==== visual elements ====
  
  # MM's colors
  my_cols = c("NS" = "grey10", "IL1B" = "#FC4E07", "Bud" = "#E7B800", "IL1B + Bud" = "#2E9FDF",
              "IL17A" = "#FC4E07", "IL17A + Bud" = "#2E9FDF",
              "TNF" = "#FC4E07", "Dex" = "#E7B800", "IL1B + Dex" = "#2E9FDF",
              "Indacaterol" = "#CC79A7", "Mometazone" = "#E7B800", "Mometazone + Indacaterol" = "#0072B2",
              "Form" = "#CC79A7", "TNF + Form" = "#0072B2", "Fsk" = "grey50", "Bud + Form" = "#0072B2",
              "naive" = "grey10", "Ad-GFP" = "#2E9FDF", "Ad-DUSP1" = "firebrick", "Ad-IkBa" = "firebrick", 
              "Control" = "#0072B2", "CMV" = "firebrick", "Vehicle" = "#0072B2", "HRV16" = "firebrick",
              "WT" = "grey50", "AKO" = "firebrick1", "BKO" = "steelblue1", "DKO" = "#CC79A7",
              "Tap" = "firebrick", "Bay" = "#2E9FDF", "Vil" = "#E7B800")
  
  # MM's theme
  my_theme = theme_bw(base_size = 12) + 
    theme(legend.position= "right", 
          legend.justification = "top", 
          legend.title = element_text(colour = "grey10", size = 12, face = 2),
          legend.text = element_text(colour = "grey10", size = 12, margin = margin(r = 5, unit = "pt")),
          legend.key.size = unit(15, 'points'),
          axis.text = element_text(colour = "grey10", size = 12),
          panel.border = element_rect(colour = "grey10", linewidth = 0.6),
          panel.grid = element_line(linewidth = 0.25),
          plot.background = element_rect(fill = 'transparent', color = "transparent"),
          plot.margin = unit(c(0, 2, 0, 2), "points"),
          strip.text = element_text(size = 14),
          strip.background = element_rect(fill = "transparent", colour = "transparent", linewidth = 0.5)) 
  
  ## ==== registries ====
  
  # dataset type registry
  dataset_meta <- list(
    "A549 IL1B Bud"                          = "single_gene",
    "A549 IL1B Bud (transcripts)"            = "single_gene",
    "HBE IL1B IFNg Dex"                      = "single_gene",
    "HBE TNF Form"                           = "single_gene",
    "ALI IL1B Bud"                           = "single_gene",
    "ALI Bud Form"                           = "single_gene",
    "ALI Cig IL17A Bud"                      = "single_gene",
    "ALI CMV"                                = "single_gene",
    "ALI HRV (Bai 2015)"                     = "single_gene",
    "BEAS-2B Mom Ind"                        = "single_gene",
    "A549 Ad-DUSP1 IL1B"                     = "single_gene",
    "A549 Ad-IkBa IL1B Dex"                  = "single_gene",
    "BEAS-2B AKO/BKO Form Bud"               = "single_gene",
    "BEAS-2B AKO/BKO Form Bud (transcripts)" = "single_gene",
    "BEAS-2B dKO Form"                       = "single_gene",
    "BEAS-2B dKO Form (transcripts)"         = "single_gene",
    "BEAS-2B Tap Bay Vil"                    = "single_gene",
    "HBE CMV timecourse (Parkins Lab)"       = "single_gene",
    "Basal expression - single gene"         = "single_gene",
    "Basal expression - heatmap"             = "multi_gene",
    "BEAS-2B RNO ONO Vil (GSE267218)"        = "single_gene"
  )
  
  # dataset file registry
  dataset_map <- list(
    "A549 IL1B Bud"                          = "data/a549_ib clean.rds",
    "A549 IL1B Bud (transcripts)"            = NULL, # handled via chunks
    "HBE IL1B IFNg Dex"                      = "data/hbe_iid clean.rds",
    "HBE TNF Form"                           = "data/hbe_tf clean.rds",
    "ALI IL1B Bud"                           = "data/ali_ib clean.rds",
    "ALI Bud Form"                           = "data/ali_bf clean.rds",
    "ALI Cig IL17A Bud"                      = "data/ali_cig clean.rds",
    "ALI CMV"                                = "data/ali_cmv clean.rds",
    "ALI HRV (Bai 2015)"                     = "data/ali_hrv clean.rds",
    "BEAS-2B Mom Ind"                        = "data/b2b_mi clean.rds",
    "A549 Ad-DUSP1 IL1B"                     = "data/dusp clean.rds",
    "A549 Ad-IkBa IL1B Dex"                  = "data/ikba clean.rds",
    "BEAS-2B AKO/BKO Form Bud"               = "data/b2b_pka clean.rds",
    "BEAS-2B AKO/BKO Form Bud (transcripts)" = NULL, # handled via chunks
    "BEAS-2B dKO Form"                       = "data/b2b_dko clean.rds",
    "BEAS-2B dKO Form (transcripts)"         = NULL, # handled via chunks
    "BEAS-2B Tap Bay Vil"                    = "data/b2b_tap_bay_vil clean.rds",
    "HBE CMV timecourse (Parkins Lab)"       = "data/hbe_cmv clean.rds",
    "Basal expression - single gene"         = "data/all_cells clean.rds",
    "Basal expression - heatmap"             = "data/all_cells clean.rds",
    "BEAS-2B RNO ONO Vil (GSE267218)"        = "data/b2b_rno_ono_vil clean.rds"
  )
  
  # chunk mappings
  chunk_map <- list(
    "A-E" = LETTERS[1:5],
    "F-J" = LETTERS[6:10],
    "K-O" = LETTERS[11:15],
    "P-T" = LETTERS[16:20],
    "U-Z" = LETTERS[21:26]
  )
  
  ## ==== functions ====
  
  # function to pick the correct dataset type
  dataset_type <- function(dataset) {
    dataset_meta[[dataset]]
  }
  
  # function to pick the correct chunk
  get_transcript_chunk <- function(gene, dataset_prefix) {
    
    # Determine chunk based on first letter
    first_letter <- toupper(substr(gene, 1, 1))
    chunk <- switch(first_letter,
                    "A" = "A-E", "B" = "A-E", "C" = "A-E", "D" = "A-E", "E" = "A-E",
                    "F" = "F-J", "G" = "F-J", "H" = "F-J", "I" = "F-J", "J" = "F-J",
                    "K" = "K-O", "L" = "K-O", "M" = "K-O", "N" = "K-O", "O" = "K-O",
                    "P" = "P-T", "Q" = "P-T", "R" = "P-T", "S" = "P-T", "T" = "P-T",
                    "U" = "U-Z", "V" = "U-Z", "W" = "U-Z", "X" = "U-Z", "Y" = "U-Z", "Z" = "U-Z",
                    stop("Gene not recognized"))
    
    # Build full path to RDS
    path <- paste0(dataset_prefix, " ", chunk, ".rds")
    
    if (!file.exists(path)) stop("File does not exist: ", path)
    
    # Read and return
    readRDS(path)
  }
  
  # function to parse genes
  parse_genes <- function(x) {
    x %>%
      strsplit("[,\n]") %>%
      unlist() %>%
      trimws() %>%
      toupper() %>%
      (\(.) .[. != "" & !is.na(.)])()
  }
  
  # function to build a heatmap (avoids monster if_else tree)
  
  render_heatmap <- function(dat, metric) {
    
    mat <- dat %>% 
      select("Gene", "sample", metric) %>%
      mutate(Gene = factor(Gene, levels = unique(state()$genes))) %>%
      arrange(Gene) %>%
      pivot_wider(names_from = "sample", values_from = metric) %>%
      column_to_rownames("Gene") %>%
      as.matrix()
    
    colsplit = factor(
      c(
        rep("A549", 4),
        rep("ALI", 5),
        rep("B2B", 4),
        rep("Brushings", 3),
        rep("HBE", 5)
      ),
      levels = c("A549", "ALI", "B2B", "Brushings", "HBE")
    )
    
    vals <- as.numeric(mat)
    vals <- vals[is.finite(vals)]
    
    lo  <- quantile(vals, 0.02)
    hi  <- quantile(vals, 0.98)
    mid <- median(vals)
    
    if (lo < 0 && hi > 0) {
      # diverging
      max_abs <- max(abs(c(lo, hi)))
      
      col_fun <- circlize::colorRamp2(
        c(-max_abs, -0.25 * max_abs, 0,
          0.25 * max_abs,  max_abs),
        c("steelblue4",
          "steelblue2",
          "white",
          "firebrick2",
          "firebrick4")
      )
      
    } else if (hi <= 0) {
      # all negative
      col_fun <- circlize::colorRamp2(
        c(lo, mid, hi),
        c("steelblue4", "steelblue2", "white")
      )
      
    } else {
      # all positive
      col_fun <- circlize::colorRamp2(
        c(lo, mid, hi),
        c("white", "firebrick2", "firebrick4")
      )
    }
    
    hm <- Heatmap(mat,
                  cluster_rows = F, 
                  cluster_columns = F,
                  column_split = colsplit,
                  rect_gp = grid::gpar(col = "black", lwd = 0.5),
                  show_row_names = TRUE,
                  show_column_names = FALSE,
                  row_names_side = "left",
                  column_title_rot = 45,
                  width  = grid::unit(0.75 * ncol(mat), "cm"),
                  height = grid::unit(0.75 * nrow(mat), "cm"),
                  name = metric,
                  col = col_fun
    )
    
    draw(hm)
      
  }
  
  ## ==== snapshot choices and define state ====
  
  # state <- eventReactive(input$button, {
  #   list(
  #     dataset = input$dataset,
  #     gene    = toupper(input$gene),
  #     metric  = input$metric
  #   )
  # })
  
  state <- eventReactive(input$button, {
    type <- dataset_type(input$dataset)
    
    genes <- if (type == "multi_gene") {
      parse_genes(input$genes)
    } else {
      toupper(input$gene)
    }
    
    list(
      dataset = input$dataset,
      genes   = genes,
      metric  = input$metric,
      type    = type
    )
  })
  
  ## ==== input ui loader ====
  
  output$gene_input <- renderUI({
    req(input$dataset)
    
    if (dataset_type(input$dataset) == "multi_gene") {
      multi_gene_ui
    } else {
      single_gene_ui
    }
  })
  
  ## ==== data loader ====
  
  x <- reactive({
    
    req(state())
    
    # gene <- toupper(state()$gene)
    genes <- state()$genes
    dataset <- state()$dataset
    
    # Check if dataset is transcript-level and needs chunked loading
    if (dataset %in% c("A549 IL1B Bud (transcripts)",
                       "BEAS-2B AKO/BKO Form Bud (transcripts)",
                       "BEAS-2B dKO Form (transcripts)")) {
      
      # Determine dataset prefix (for path) from dataset name
      prefix <- switch(dataset,
                       "A549 IL1B Bud (transcripts)"            = "data/a549_ib_transcripts clean",
                       "BEAS-2B AKO/BKO Form Bud (transcripts)" = "data/b2b_pka_transcripts clean",
                       "BEAS-2B dKO Form (transcripts)"         = "data/b2b_dko_transcripts clean")
      
      # Load correct chunk
      gene <- genes[[1]]
      dat <- get_transcript_chunk(gene, prefix)
    } 
    
    else {
      # Regular single-file dataset
      rds_path <- dataset_map[[dataset]]
      req(rds_path)
      dat <- readRDS(rds_path)
    }
    
    # Filter for the selected gene
    # dat %>% filter(Gene == gene)
    dat %>% filter(Gene %in% genes)
  })
  
  y_label <- reactive({
    switch(state()$metric,
           "log2_tpm"  = expression(log[2]*" TPM"),
           "tpm"       = "TPM",
           "log2_fold" = expression(log[2]*" fold"),
           "fold"      = "fold",
           state()$metric)  # fallback
  })
  
  ## ==== plot and table loader ====
  
  output$expression_plot <- renderPlot({
    
    req(state(), x())
    
    # defer to a heatmap function if a multi-gene option is picked
    if (state()$type == "multi_gene") {
      req(length(state()$genes) > 0)
      return(render_heatmap(x(), metric = state()$metric))
    }
    
    # otherwise use MM's monster tree
    if (state()$dataset == "A549 IL1B Bud") {
      x() %>% 
        ggplot(aes(time, .data[[state()$metric]], group = treatment, color = treatment)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1)+
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), size = 0.5, show.legend = F) + 
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.5, show.legend = F, geom = "errorbar", width = 0.2) +  
        scale_color_manual(values = my_cols) +  
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        #expand_limits(y = c(0)) +
        facet_wrap(~Gene)+
        labs(y = y_label(), x = "time (h)", color = "Treatment") + 
        my_theme + theme(aspect.ratio = 1)
    }
    else if (state()$dataset == "A549 IL1B Bud (transcripts)") {
      x() %>% 
        ggplot(aes(time, .data[[state()$metric]], group = treatment, color = treatment)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1)+
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), size = 0.5, show.legend = F) + 
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.5, show.legend = F, geom = "errorbar", width = 0.2) +  
        scale_color_manual(values = my_cols) +  
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        #expand_limits(y = c(0)) +
        facet_wrap(~transcript, ncol = 2)+
        labs(y = y_label(), x = "time (h)", color = "Treatment", title = paste(state()$genes, collapse = ", ")) + 
        my_theme + theme(aspect.ratio = 1, legend.position = "top")
    } 
    else if (state()$dataset == "BEAS-2B Mom Ind") {
      x() %>% 
        ggplot(aes(time, .data[[state()$metric]], group = treatment, color = treatment)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1)+
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), size = 0.5, show.legend = F) + 
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.5, show.legend = F, geom = "errorbar", width = 0.2) +  
        scale_color_manual(values = my_cols) +  
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        #expand_limits(y = c(0)) +
        facet_wrap(~Gene)+
        labs(y = y_label(), x = "time (h)", color = "Treatment") + 
        my_theme + theme(aspect.ratio = 1)
    }
    else if (state()$dataset == "HBE IL1B IFNg Dex") {
      x() %>% 
        ggplot(aes(time, .data[[state()$metric]], group = ttt, color = ttt)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1)+
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), size = 0.5, show.legend = F) + 
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.5, show.legend = F, geom = "errorbar", width = 0.2) +  
        scale_color_manual(values = my_cols) +  
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        #expand_limits(y = c(0)) +
        facet_grid(Gene~ifn)+
        labs(y = y_label(), x = "time (h)", color = "Treatment") + 
        my_theme + theme(aspect.ratio = 1.5)
    }
    else if (state()$dataset == "HBE TNF Form") {
      x() %>% 
        ggplot(aes(time, .data[[state()$metric]], group = treatment, color = treatment)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1)+
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), size = 0.5, show.legend = F) + 
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.5, show.legend = F, geom = "errorbar", width = 0.2) +  
        scale_color_manual(values = my_cols) +  
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        #expand_limits(y = c(0)) +
        facet_wrap(~Gene)+
        labs(y = y_label(), x = "time (h)", color = "Treatment") + 
        my_theme + theme(aspect.ratio = 1)
    }
    else if (state()$dataset == "ALI IL1B Bud") {
      x() %>% 
        ggplot(aes(treatment, .data[[state()$metric]], fill = treatment))+
        geom_boxplot(coef = 3, outlier.size = 0, outlier.stroke = 0, size = 0.5, show.legend = T)  + 
        labs(y = y_label(), fill = "Treatment", x= NULL) + 
        scale_fill_manual(values = my_cols)+
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        facet_wrap(~Gene)+
        my_theme + theme(aspect.ratio = 1.5, axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    else if (state()$dataset == "ALI Bud Form") {
      x() %>% 
        ggplot(aes(time, .data[[state()$metric]], group = treatment, color = treatment)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1)+
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), size = 0.5, show.legend = F) + 
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.5, show.legend = F, geom = "errorbar", width = 0.2) +  
        scale_color_manual(values = my_cols) +  
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        #expand_limits(y = c(0)) +
        facet_wrap(~Gene)+
        labs(y = y_label(), x = "time (h)", color = "Treatment") + 
        my_theme + theme(aspect.ratio = 2)
    }
    else if (state()$dataset == "ALI Cig IL17A Bud") {
      x() %>% 
        ggplot(aes(treatment, .data[[state()$metric]], group = treatment, fill = treatment)) +
        geom_boxplot(aes(linetype = condition),color = "grey10", linewidth = 0.4, outlier.shape = NULL, 
                     outlier.size = 0, outlier.stroke = 0, coef = 3) +
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        scale_fill_manual(values = my_cols)+
        facet_grid(Gene~cond, scales = "free_y")+
        labs(y = y_label(), x = NULL, fill = NULL) + 
        my_theme + theme(aspect.ratio = 2,axis.ticks.x = element_blank(), axis.text.x = element_blank()) 
    }
    else if (state()$dataset == "ALI CMV") {
      x() %>% 
        ggplot(aes(condition, .data[[state()$metric]], fill = condition)) +
        geom_boxplot(color = "grey10", linewidth = 0.4, outlier.shape = NULL, outlier.size = 0, outlier.stroke = 0, coef = 3) +
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        scale_fill_manual(values = my_cols)+
        facet_wrap(Gene~.)+
        labs(y = y_label(), x = NULL, fill = NULL) + 
        my_theme + theme(aspect.ratio = 2,axis.ticks.x = element_blank(), axis.text.x = element_blank()) 
    }
    else if (state()$dataset == "ALI HRV (Bai 2015)") {
      x() %>% 
        ggplot(aes(condition, .data[[state()$metric]], fill = condition)) +
        geom_boxplot(aes(linetype = disease_state), color = "grey10", linewidth = 0.4, outlier.shape = NULL, 
                     outlier.size = 0, outlier.stroke = 0, coef = 3) +
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        scale_fill_manual(values = my_cols)+
        facet_grid(Gene~disease_state)+
        labs(y = y_label(), x = NULL, fill = NULL, linetype = NULL) + 
        my_theme + theme(aspect.ratio = 1.5, strip.text.x = element_text(size = 9),
                         axis.ticks.x = element_blank(), axis.text.x = element_blank()) 
    }
    else if (state()$dataset == "A549 Ad-DUSP1 IL1B") {
      x() %>% 
        ggplot(aes(time, .data[[state()$metric]], group = interaction(treatment,condition) , color = treatment)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1, alpha = 0.8) +
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), stroke = 0.5, size = 0.4, show.legend = F, alpha = 0.8) +  
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.3, show.legend = F, geom = "errorbar", width = 0.2, alpha = 0.8) +  
        scale_color_manual(values = my_cols) + 
        facet_grid(Gene~condition)+
        labs(y = y_label(), x= "time (h)" , color = "Treatment") + 
        my_theme + theme(aspect.ratio = 2) 
      
    }
    else if (state()$dataset == "A549 Ad-IkBa IL1B Dex") {
      x() %>% 
        ggplot(aes(condition, .data[[state()$metric]], group = condition , fill = condition)) +
        geom_hline(yintercept = 0, linewidth = 0.75, alpha = 0.5, linetype = 2) +
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5)  + 
        scale_fill_manual(values = my_cols) + 
        facet_grid(Gene~treatment)+
        labs(y = y_label(), x= NULL , fill = "Condition") + 
        my_theme + theme(aspect.ratio = 2, axis.ticks.x = element_blank(), axis.text.x = element_blank()) 
    }
    else if (state()$dataset == "BEAS-2B AKO/BKO Form Bud") {
      x() %>% 
        ggplot(aes(x = condition, y = .data[[state()$metric]], fill = condition)) +
        facet_grid(Gene~treatment) +
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5) +
        geom_point() +
        labs(y = y_label(), x = NULL, fill = "Condition") +
        scale_fill_manual(values = my_cols) +
        my_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 2)
    }
    else if (state()$dataset == "BEAS-2B AKO/BKO Form Bud (transcripts)") {
      x() %>% 
        mutate(condition = factor(condition, levels = c("WT", "AKO", "BKO"))) %>% #faceting error, corrected here
        ggplot(aes(x = condition, y = .data[[state()$metric]], fill = condition)) +
        facet_grid(transcript~treatment) +
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5) +
        geom_point() +
        labs(y = y_label(), x = NULL, fill = "Condition", title = paste(state()$genes, collapse = ", ")) +
        scale_fill_manual(values = my_cols) +
        my_theme + 
        theme(axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              aspect.ratio = 1,
              strip.text.y = element_text(size = 8))
    }
    else if (state()$dataset == "BEAS-2B dKO Form") {
      x() %>% 
        ggplot(aes(x = condition, y = .data[[state()$metric]], fill = condition)) +
        facet_grid(Gene~treatment) +
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5) +
        geom_point() +
        labs(y = y_label(), x = NULL, fill = "Condition") +
        scale_fill_manual(values = my_cols) +
        my_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 2)
    }
    else if (state()$dataset == "BEAS-2B dKO Form (transcripts)") {
      x() %>% 
        ggplot(aes(x = condition, y = .data[[state()$metric]], fill = condition)) +
        facet_grid(transcript~treatment) +
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5) +
        geom_point() +
        labs(y = y_label(), x = NULL, fill = "Condition", title = paste(state()$genes, collapse = ", ")) +
        scale_fill_manual(values = my_cols) +
        my_theme +
        theme(axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              aspect.ratio = 1,
              strip.text.y = element_text(size = 8))
    }
    else if (state()$dataset == "BEAS-2B Tap Bay Vil") {
      x() %>% 
        ggplot(aes(x = treatment, y = .data[[state()$metric]], fill = treatment)) +
        facet_wrap(~Gene) +
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5) +
        geom_point() +
        labs(y = y_label(), x = NULL, fill = "Treatment") +
        scale_fill_manual(values = c("NS" = "grey50", "Tap" = "firebrick", "Bay" = "#2E9FDF", "Vil" = "#E7B800")) +
        theme_bw() +
        my_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 1.5)
    }
    else if (state()$dataset == "HBE CMV timecourse (Parkins Lab)") {
      x() %>% 
        ggplot(aes(x = time, y = .data[[state()$metric]], color = condition, group = condition)) +
        stat_summary(geom = "point", fun = "mean") +
        stat_summary(geom = "line", fun = "mean") +
        stat_summary(geom = "errorbar", fun.data = "mean_se", fun.args = list(mult = 1),
                     width = 0.2) +
        facet_wrap(~Gene) +
        labs(x = "time (h)", y = y_label(), color = "Condition") +
        scale_color_manual(values = my_cols, labels = c("CON" = "control", "CMV" = "CMV")) +
        my_theme + theme(aspect.ratio = 1)
    } 
    else if (state()$dataset == "Basal expression - single gene") {
      x() %>% 
        ggplot(aes(cells, .data[[state()$metric]]))+
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5)  + 
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        labs(y = y_label(), x = NULL) + 
        expand_limits(y = 0) +
        facet_wrap(~Gene)+
        my_theme + 
        theme(aspect.ratio = 0.8 , 
              legend.position="none", 
              axis.text.x = element_text(angle = 45, hjust = 1))
    }
    else if (state()$dataset == "BEAS-2B RNO ONO Vil (GSE267218)") {
      x() %>%
        ggplot(aes(x = treatment, y = .data[[state()$metric]])) +
        facet_wrap(~Gene) +
        geom_boxplot(outliers = F) +
        geom_point() +
        my_theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              aspect.ratio = 0.5) +
        labs(x = NULL)
    }
    
    
    else {output$text = renderPrint("select a valid dataset")} 
  }, 
  
  res = 120,
  
  height = function() {
    dat <- x()
    
    if ("transcript" %in% colnames(dat)) {
      n <- dplyr::n_distinct(dat$transcript)
      max(800, n * 180)
    } 
    
    else if (length(unique(dat$Gene)) > 1){
      800
    }
    
    else {
      400
    }
  }
  )
  
  output$table <- renderTable({x()})
  

  }

# ==== create shiny app ====
shinyApp(ui = ui, server = server)





















