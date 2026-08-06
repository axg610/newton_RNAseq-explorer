library(shiny)
library(shinycssloaders)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(magick)
library(tidyverse)

# When adding a new dataset:
# GENE LOOKUP MODE:
# - add dataset name to lookup ui
# - declare data type (single_gene or multiple) in registries
# - declare data location (data/...) in registries
# - add plot to MM's monster if-else tree
# INTERACTIVE HEATMAP MODE:
# - add dataset name to lookup ui
# - map name to file location in the interactive heatmap loader
# - add experimental details as html in interactive heatmap loader


# ==== define ui ====

ui <- fluidPage(
  
  titlePanel("Newton/Giembycz Lab RNAseq Explorer"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      # select mode
      tabsetPanel(
        id = "app_mode",
        type = "tabs",
        
        tabPanel(
          "Gene search mode",
          value = "lookup"
        ),
        
        tabPanel(
          "Interactive heatmap mode",
          value = "interactive"
        )
      ),
      
      # boot lookup ui
      conditionalPanel(
        condition = "input.app_mode == 'lookup'",
        
        selectInput(
          "dataset", "Dataset:",
          choices = list(
            "A549 datasets" = c(
              "A549 IL1B Bud - single gene",
              "A549 IL1B Bud - heatmap",
              # "A549 IL1B Bud (transcripts) - single gene",
              "A549 Ad-DUSP1 IL1B",
              "A549 Ad-IkBa IL1B Dex 6h",
              "A549 Ad-IkBa IL1B Bud timecourse",
              "A549 IRF1-KO IL1B 6h",
              "A549 IRF1 siRNA IL1B 6h",
              "A549 TAK1-KO IL1B 2h"
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
              # "BEAS-2B AKO/BKO Form Bud (transcripts)",
              "BEAS-2B dKO Form",
              # "BEAS-2B dKO Form (transcripts)",
              "BEAS-2B Mom Ind",
              "BEAS-2B Tap Bay Vil",
              "BEAS-2B RNO ONO Vil (GSE267218)",
              "BEAS-2B RNO Salm (GSE126981)",
              "BEAS-2B CRTC-KO Form 2h"
            ),
            "HBE datasets" = c(
              "HBE IL1B IFNg Dex",
              "HBE TNF Form",
              "HBE CMV timecourse (Parkins Lab)"
            ),
            "Datasets with 2+ celltypes" = c(
              "Basal expression - single gene",
              "Basal expression - heatmap",
              "ALI vs HBE vs BEAS-2B Fsk - single gene",
              "ALI vs HBE vs BEAS-2B Fsk - heatmap"
            )
          )
        ),
        
        uiOutput("gene_input"),
        
        radioButtons(
          "metric",
          "Y-axis:",
          choices = c("log2_tpm", "log2_fold", "tpm", "fold"),
          selected = "log2_tpm"
        ),
        
        actionButton("button", "Generate", icon = icon("redo"))
        
      ),
      
      # boot interactive heatmaps ui
      conditionalPanel(
        condition = "input.app_mode == 'interactive'",
        
        selectInput(
          "interactive_heatmap",
          "Interactive heatmap:",
          choices = c(
            "IL1B-regulated genes: A549/ALI/HBE @ 6 h",
            "IL1B-regulated genes: naive/Ad-GFP/Ad-IKBA A549 timecourse",
            "IL1B-regulated genes: IRF1 KO @ 6 h", 
            "IL1B-regulated genes: IRF1 siRNA @ 6 h",
            "IL1B-regulated genes: TAK1 KO and rescue @ 2 h",
            "Form-regulated genes: CRTC KOs @ 2 h"
          )
        ),
        
        actionButton(
          "launch_heatmap",
          "Launch interactive heatmap",
          icon = icon("th")
        ),
        
        helpText(
          HTML("Heatmaps may take ~30 seconds to launch.")
        )
      ),
      
      # static helptext
      helpText(
        HTML(
          '<br/>
          <a href="https://drive.google.com/drive/folders/1iumRvf95fkqg2FvuGiZunExZcJpQr5lN?usp=sharing" target="_blank"><b>Data repository</b></a><br/>
     <a href="https://github.com/axg610/newton_RNAseq-explorer" target="_blank"><b>Code repository</b></a>
          <br/>'
        )
      ),
      
      wellPanel(
        tags$details(
          tags$summary(tags$b("Change log")),
          HTML(
            "
          
<b>v3.0.2, Aug 6 2026</b><br/>
  - Added IRF1 KO dataset and heatmap <br/>
  - Added IRF1 siRNA dataset and heatmap <br/>
  - Added TAK1 KO dataset and heatmap <br/>
  - Added CRTC KO dataset and heatmap <br/>
  - Archived transcript-level browsers for room <br/>

<b>v3.0.1, May 21 2026</b><br/>
  - Added Ad-IKBA dataset and heatmap <br/>

<b>v3.0, May 11 2026</b><br/>
  - Added interactive heatmaps <br/>

<b>v2.3a, May 5 2026</b><br/>
  - Added heatmap option to A549 IL1B Bud <br/>
  - Upgrades to internal heatmap functions <br/>

<b>v2.3, May 3 2026</b><br/>
  - Added ALI/HBE/B2B forskolin comparison <br/>

<b>v2.2b, April 21 2026</b><br/>
  - Synonym genes no longer break search (search SIK1 on HBE TNF Form) <br/>
  - Set unpublished data to private <br/>

<b>v2.2a, February 25 2026</b><br/>
  - Added Radhika's data (GSE267218 and GSE126981) <br/>

<b>v2.2, February 7 2026</b><br/>
  - Added experimental heatmap option for basal expression data <br/>

<b>v2.1, December 2025</b><br/>
  - Plots now update only when Generate is clicked <br/>
  - Added transcript-level datasets <br/>
  - Added log2fold, tpm, and fold options <br/>
  - Added loading indicator <br/>
  - Datasets cleaned and consolidated <br/>
  - Large datasets chunked to avoid OOM crash

    "
            
          )
        )
      )
    ),
    
    mainPanel(
      
      # ==== gene lookup mode ====
      conditionalPanel(
        condition = "input.app_mode == 'lookup'",
        
        tabsetPanel(
          type = "tabs", 
          
          tabPanel(
            "Plot",
            
            plotOutput("expression_plot") %>% 
              withSpinner(
                type = 6, 
                color = "#555555", 
                caption = "Working..."
              )
          ), 
          
          tabPanel(
            "Table",
            tableOutput("table")
          )
        )
      ),
      
      # ==== interactive heatmap mode ====
      
      # should get a spinner working here
      conditionalPanel(
        condition = "input.app_mode == 'interactive'",
        
        tagList(
          
          uiOutput("heatmap_note"),
          
          br(),
          
          uiOutput("interactive_heatmap_container")
        )
      ),
      
      conditionalPanel(
        condition = "input.app_mode == 'interactive' && input.launch_heatmap > 0",
        
        InteractiveComplexHeatmapOutput(
          output_ui = "heatmap_ui",
          response_ui = "heatmap_response",
          width1 = 500,
          height1 = 1400
        )
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
  value = paste(
    readLines("default_multi_genes.txt"),
    collapse = "\n"
  ),
  rows = 6
)

# ==== define server ====
server <- function(input, output, session) {
  
  ## ==== visual elements ====
  
  # MM's colors
  my_cols = c(
    "NS" = "grey10", "IL1B" = "#FC4E07", 
    "Bud" = "#E7B800", "IL1B + Bud" = "#2E9FDF", "IB" = "#2E9FDF",
    "IL17A" = "#FC4E07", "IL17A + Bud" = "#2E9FDF",
    "TNF" = "#FC4E07", "Dex" = "#E7B800", "IL1B + Dex" = "#2E9FDF",
    "Indacaterol" = "#CC79A7", "Mometazone" = "#E7B800", 
    "Mometazone + Indacaterol" = "#0072B2",
    "Form" = "#CC79A7", "TNF + Form" = "#0072B2", 
    "Fsk" = "grey50", "Bud + Form" = "#0072B2",
    "naive" = "grey10", "Ad-GFP" = "#2E9FDF", 
    "Ad-DUSP1" = "firebrick", "Ad-IkBa" = "firebrick", 
    "Control" = "#0072B2", "CMV" = "firebrick", 
    "Vehicle" = "#0072B2", "HRV16" = "firebrick",
    "WT" = "grey50", "AKO" = "firebrick1", 
    "BKO" = "steelblue1", "DKO" = "#CC79A7",
    "Tap" = "firebrick", "Bay" = "#2E9FDF", 
    "Vil" = "#E7B800"
  )
  
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
    "A549 IL1B Bud - single gene"               = "single_gene",
    "A549 IL1B Bud - heatmap"                   = "multi_gene",
    # "A549 IL1B Bud (transcripts) - single gene" = "single_gene",
    "HBE IL1B IFNg Dex"                         = "single_gene",
    "HBE TNF Form"                              = "single_gene",
    "ALI IL1B Bud"                              = "single_gene",
    "ALI Bud Form"                              = "single_gene",
    "ALI Cig IL17A Bud"                         = "single_gene",
    "ALI CMV"                                   = "single_gene",
    "ALI HRV (Bai 2015)"                        = "single_gene",
    "BEAS-2B Mom Ind"                           = "single_gene",
    "A549 Ad-DUSP1 IL1B"                        = "single_gene",
    "A549 Ad-IkBa IL1B Dex 6h"                  = "single_gene",
    "BEAS-2B AKO/BKO Form Bud"                  = "single_gene",
    # "BEAS-2B AKO/BKO Form Bud (transcripts)"    = "single_gene",
    "BEAS-2B dKO Form"                          = "single_gene",
    # "BEAS-2B dKO Form (transcripts)"            = "single_gene",
    "BEAS-2B Tap Bay Vil"                       = "single_gene",
    "HBE CMV timecourse (Parkins Lab)"          = "single_gene",
    "Basal expression - single gene"            = "single_gene",
    "Basal expression - heatmap"                = "multi_gene",
    "BEAS-2B RNO ONO Vil (GSE267218)"           = "single_gene",
    "BEAS-2B RNO Salm (GSE126981)"              = "single_gene",
    "ALI vs HBE vs BEAS-2B Fsk - single gene"   = "single_gene",
    "ALI vs HBE vs BEAS-2B Fsk - heatmap"       = "multi_gene",
    "A549 Ad-IkBa IL1B Bud timecourse"          = "single_gene",
    "A549 IRF1-KO IL1B 6h"                      = "single_gene",
    "A549 IRF1 siRNA IL1B 6h"                   = "single_gene",
    "A549 TAK1-KO IL1B 2h"                      = "single_gene",
    "BEAS-2B CRTC-KO Form 2h"                   = "single_gene"
  )
  
  # dataset file registry
  dataset_map <- list(
    "A549 IL1B Bud - single gene"                = "data/a549_ib clean.rds",
    "A549 IL1B Bud - heatmap"                    = "data/a549_ib clean.rds",
    # "A549 IL1B Bud (transcripts) - single gene"  = NULL, # handled via chunks
    "HBE IL1B IFNg Dex"                          = "data/hbe_iid clean.rds",
    "HBE TNF Form"                               = "data/hbe_tf clean.rds",
    "ALI IL1B Bud"                               = "data/ali_ib clean.rds",
    "ALI Bud Form"                               = "data/ali_bf clean.rds",
    "ALI Cig IL17A Bud"                          = "data/ali_cig clean.rds",
    "ALI CMV"                                    = "data/ali_cmv clean.rds",
    "ALI HRV (Bai 2015)"                         = "data/ali_hrv clean.rds",
    "BEAS-2B Mom Ind"                            = "data/b2b_mi clean.rds",
    "A549 Ad-DUSP1 IL1B"                         = "data/dusp clean.rds",
    "A549 Ad-IkBa IL1B Dex 6h"                   = "data/ikba clean.rds",
    "BEAS-2B AKO/BKO Form Bud"                   = "data/b2b_pka clean.rds",
    # "BEAS-2B AKO/BKO Form Bud (transcripts)"     = NULL, # handled via chunks
    "BEAS-2B dKO Form"                           = "data/b2b_dko clean.rds",
    # "BEAS-2B dKO Form (transcripts)"             = NULL, # handled via chunks
    "BEAS-2B Tap Bay Vil"                        = "data/b2b_tap_bay_vil clean.rds",
    "HBE CMV timecourse (Parkins Lab)"           = "data/hbe_cmv clean.rds",
    "Basal expression - single gene"             = "data/all_cells clean.rds",
    "Basal expression - heatmap"                 = "data/all_cells clean.rds",
    "BEAS-2B RNO ONO Vil (GSE267218)"            = "data/b2b_rno_ono_vil clean.rds",
    "BEAS-2B RNO Salm (GSE126981)"               = "data/b2b_rno_salm clean.rds",
    "ALI vs HBE vs BEAS-2B Fsk - single gene"    = "data/ali_b2b_hbe_fsk clean.rds",
    "ALI vs HBE vs BEAS-2B Fsk - heatmap"        = "data/ali_b2b_hbe_fsk clean.rds",
    "A549 Ad-IkBa IL1B Bud timecourse"           = "data/ikba_timecourse clean.rds",
    "A549 IRF1-KO IL1B 6h"                       = "data/a549_irf1ko clean.rds",
    "A549 IRF1 siRNA IL1B 6h"                    = "data/a549_irf1_siRNA clean.rds",
    "A549 TAK1-KO IL1B 2h"                       = "data/a549_tak1ko clean.rds",
    "BEAS-2B CRTC-KO Form 2h"                    = "data/b2b_crtc_ko clean.rds"
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
  
  # function to build basal cell types heatmap (avoids monster if_else tree)
  
  render_heatmap <- function(
    dat, 
    colsplit, 
    metric, 
    vertical.scale = 1,
    horizontal.scale = 1
    ) {
    
    mat <- dat %>% 
      select("Gene", "sample", metric) %>%
      mutate(Gene = factor(Gene, levels = unique(state()$genes))) %>%
      arrange(Gene) %>%
      pivot_wider(names_from = "sample", values_from = metric) %>%
      column_to_rownames("Gene") %>%
      as.matrix()
    
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
    
    hm <- Heatmap(
      mat,
      cluster_rows = F, 
      cluster_columns = F,
      column_split = colsplit,
      rect_gp = grid::gpar(col = "black", lwd = 0.5),
      show_row_names = TRUE,
      show_column_names = FALSE,
      row_names_side = "left",
      column_title_rot = 45,
      width  = grid::unit(horizontal.scale * 0.75 * ncol(mat), "cm"),
      height = grid::unit(vertical.scale * 0.75 * nrow(mat), "cm"),
      name = metric,
      col = col_fun
    )
    
    draw(hm)
    
  }
  
  ## ==== snapshot choices and define state ====
  
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
    if (dataset %in% c("A549 IL1B Bud (transcripts) - single gene",
                       "BEAS-2B AKO/BKO Form Bud (transcripts)",
                       "BEAS-2B dKO Form (transcripts)")) {
      
      # Determine dataset prefix (for path) from dataset name
      prefix <- switch(dataset,
                       "A549 IL1B Bud (transcripts) - single gene" = "data/a549_ib_transcripts clean",
                       "BEAS-2B AKO/BKO Form Bud (transcripts)"    = "data/b2b_pka_transcripts clean",
                       "BEAS-2B dKO Form (transcripts)"            = "data/b2b_dko_transcripts clean")
      
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
    # fuck!!! i guess we doing if_else trees again
    
    if (state()$type == "multi_gene") {
      
      req(length(state()$genes) > 0)
      
      if (state()$dataset == "Basal expression - heatmap"){
        return(
          render_heatmap(
            x(),
            metric = state()$metric,
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
          )
        )
      }
      
      else if (state()$dataset == "ALI vs HBE vs BEAS-2B Fsk - heatmap"){
        return(
          render_heatmap(
            x(),
            metric = state()$metric,
            vertical.scale = 0.6,
            horizontal.scale = 0.4,
            colsplit = factor(
              c(
                rep("ALI NS", 10),
                rep("ALI Fsk", 10),
                rep("HBE NS", 5),
                rep("HBE Fsk", 5),
                rep("B2B NS", 5),
                rep("B2B Fsk", 5)
              ),
              levels = c(
                "ALI NS", "ALI Fsk",
                "HBE NS", "HBE Fsk",
                "B2B NS", "B2B Fsk"
              )
            )
          )
        )
      }
      
      else if (state()$dataset == "A549 IL1B Bud - heatmap"){
        return(
          render_heatmap(
            x(),
            metric = state()$metric,
            vertical.scale = 0.5,
            horizontal.scale = 0.3,
            colsplit = factor(
              c(
                rep("NS 1h", 4), rep("NS 2h", 4), rep("NS 6h", 4), 
                rep("NS 12h", 4), rep("NS 24h", 4),
                
                rep("IL1B 1h", 4), rep("IL1B 2h", 4), rep("IL1B 6h", 4), 
                rep("IL1B 12h", 4), rep("IL1B 24h", 4),
                
                rep("Bud 1h", 4), rep("Bud 2h", 4), rep("Bud 6h", 4), 
                rep("Bud 12h", 4), rep("Bud 24h", 4),
                
                rep("I+B 1h", 4), rep("I+B 2h", 4), rep("I+B 6h", 4), 
                rep("I+B 12h", 3), # <- intentional
                rep("I+B 24h", 4)
              ),
              levels = c(
                "NS 1h", "NS 2h", "NS 6h", "NS 12h", "NS 24h",
                "IL1B 1h", "IL1B 2h", "IL1B 6h", "IL1B 12h", "IL1B 24h",
                "Bud 1h", "Bud 2h", "Bud 6h", "Bud 12h", "Bud 24h",
                "I+B 1h", "I+B 2h", "I+B 6h", "I+B 12h", "I+B 24h"
              )
            )
          )
        )
      }
      
      
    }
    
    # otherwise use MM's tree
    if (state()$dataset == "A549 IL1B Bud - single gene") {
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
    else if (state()$dataset == "A549 IL1B Bud (transcripts) - single gene") {
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
    else if (state()$dataset == "A549 Ad-IkBa IL1B Dex 6h") {
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
    else if (state()$dataset %in% c(
      "BEAS-2B RNO ONO Vil (GSE267218)",
      "BEAS-2B RNO Salm (GSE126981)" 
    )) {
      x() %>%
        ggplot(aes(x = treatment, y = .data[[state()$metric]])) +
        facet_wrap(~Gene) +
        geom_boxplot(outliers = F) +
        geom_point() +
        my_theme +
        theme(axis.text.x = element_text(angle = 315, hjust = 0),
              aspect.ratio = 0.5) +
        labs(x = NULL, y = y_label())
    }
    else if (state()$dataset == "ALI vs HBE vs BEAS-2B Fsk - single gene"){
      x() %>%
        ggplot(aes(x = treatment, y = .data[[state()$metric]])) +
        facet_grid(Gene~celltype) +
        geom_boxplot(outliers = F) +
        geom_point() +
        my_theme +
        theme(aspect.ratio = 2) +
        labs(
          x = NULL, y = y_label(), 
          caption = "***Note: this figure is a composite of three\nseparate sequencing runs (ALI, HBE, and BEAS-2B)"
        )
    }
    else if (state()$dataset == "A549 Ad-IkBa IL1B Bud timecourse"){
      x() %>%
        ggplot(
          aes(
            x = as.character(time), 
            y = .data[[state()$metric]], 
            color = treatment, 
            group = treatment
          )
        ) +
        facet_grid(
          Gene~condition,
          labeller = labeller(condition = c("Ad-IKBA" = "Ad-IκBαΔN"))) +
        stat_summary(geom = "point", fun = "mean") +
        stat_summary(geom = "line", fun = "mean") +
        stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.2) +
        my_theme +
        theme(aspect.ratio = 1.5) +
      labs(
        x = "time (h)",
        y = y_label()
      ) +
        scale_color_manual(values = my_cols, labels = c("NS", "IL1B", "Bud", "IL1B + Bud"))
    }
    else if (state()$dataset == "A549 IRF1-KO IL1B 6h"){
      x() %>%
        ggplot(
          aes(
            x = treatment,
            y = .data[[state()$metric]],
          )
        ) +
        facet_grid(
          Gene~condition,
        ) +
        geom_boxplot(outliers = F, color = "black", lwd = 0.3) +
        geom_point(shape = 21, fill = "white", color = "black") +
        my_theme +
        theme(aspect.ratio = 1.5) +
        labs(
          x = NULL,
          y = y_label()
        )
    }
    else if (state()$dataset == "A549 IRF1 siRNA IL1B 6h"){
      x() %>%
        ggplot(
          aes(
            x = group, 
            y = .data[[state()$metric]]
            ) 
        ) +
        facet_wrap(~Gene)+
        geom_boxplot(outliers = F, color = "black", lwd = 0.3) +
        geom_point(shape = 21, fill = "white", color = "black") +
        my_theme +
        theme(
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 45, hjust = 1)
          ) +
        labs(
          x = NULL,
          y = y_label()
        )
    }
    else if (state()$dataset == "A549 TAK1-KO IL1B 2h"){
      x() %>%
        ggplot(
          aes(
            x = group,
            y = .data[[state()$metric]],
            fill = treatment,
            shape = treatment
          )
        ) +
        facet_wrap(~Gene) +
        geom_boxplot(
          position = position_dodge(width = 0.7),
          width = 0.7,
          outliers = F, 
          color = "black", 
          lwd = 0.3
          ) +
        geom_point(
          position = position_dodge(width = 0.7),
          fill = "white", 
          color = "black"
          ) +
        my_theme +
        theme(
          aspect.ratio = 0.75,
          axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        labs(x = NULL, y = y_label())  +
        scale_fill_manual(values = c("grey50", "#FC4E07")) +
        scale_shape_manual(values = c(21, 24))
    }
    else if (state()$dataset == "BEAS-2B CRTC-KO Form 2h"){
      x() %>%
        ggplot(
          aes(
            x = condition,
            y = .data[[state()$metric]],
            fill = treatment,
            shape = treatment
          )
        ) +
        facet_wrap(~Gene) +
        geom_boxplot(
          position = position_dodge(width = 0.7),
          width = 0.7,
          outliers = F, 
          color = "black", 
          lwd = 0.3
        ) +
        geom_point(
          position = position_dodge(width = 0.7),
          fill = "white", 
          color = "black"
        ) +
        my_theme +
        theme(
          aspect.ratio = 0.75,
          axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        labs(x = NULL, y = y_label())  +
        scale_fill_manual(values = c("grey50", "#CC79A7")) +
        scale_shape_manual(values = c(21, 24))
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
      500
    }
  }
  )
  
  output$table <- renderTable({x()})
  
  
  
  ## ==== interactive heatmap loader ====
  
  observeEvent(input$launch_heatmap, {
    
    req(input$app_mode == "interactive")
    
    hm <- switch(
      
      input$interactive_heatmap,
      
      "IL1B-regulated genes: A549/ALI/HBE @ 6 h" =
        readRDS("interactive-heatmaps/a549 ali hbe IL1B reg genes.rds"),
      
      "IL1B-regulated genes: naive/Ad-GFP/Ad-IKBA A549 timecourse" =
        readRDS("interactive-heatmaps/a549-ikbadn-ib-timecourse-prelim-heatmap.rds"),
      
      "IL1B-regulated genes: IRF1 KO @ 6 h" =
        readRDS("interactive-heatmaps/a549-irf1-ko-prelim-heatmap.rds"),
      
      "IL1B-regulated genes: IRF1 siRNA @ 6 h" =
        readRDS("interactive-heatmaps/a549-irf1-siRNA-prelim-heatmap.rds"),
      
      "IL1B-regulated genes: TAK1 KO and rescue @ 2 h" =
        readRDS("interactive-heatmaps/a549-tak1-ko-prelim-heatmap.rds"),
      
      "Form-regulated genes: CRTC KOs @ 2 h" =
        readRDS("interactive-heatmaps/b2b-crtc-ko-prelim-heatmap.rds")
      
    )
    
    ht <- draw(hm)
    
    makeInteractiveComplexHeatmap(
      input,
      output,
      session,
      ht
    )
  }
  )
  
  # ==== notes for interactive heatmaps
  output$heatmap_note <- renderUI({
    
    req(input$interactive_heatmap)
    
    switch(
      
      input$interactive_heatmap,
      
      "IL1B-regulated genes: A549/ALI/HBE @ 6 h" =
        HTML(
          "Created to aid exploration for the IL1B paper. Row titles indicate 
          which cell types genes are regulated in.<br/>
          e.g., A = A549 only, L = ALI only, H = HBE only. ALH = regulated in all
          three celltypes, etc."
        ),
      
      "IL1B-regulated genes: naive/Ad-GFP/Ad-IKBA A549 timecourse" =
        HTML(
          "Created to aid exploration for the IL1B paper. 1404 IL1B-regulated 
          genes (abs(log2fold) >= 1, FDR <= 0.05) are shown."
        ),
      
      "IL1B-regulated genes: IRF1 KO @ 6 h" =
        HTML(
          "For Amandah's IRF1 work. 2163 IL1B-regulated genes (abs(log2fold) 
          >= 1, FDR <= 0.05) are shown."
        ),
      
      "IL1B-regulated genes: IRF1 siRNA @ 6 h" =
        HTML(
          "For Amandah's IRF1 work. 3297 IL1B-regulated genes (abs(log2fold) >=
          1, FDR <= 0.05) are shown."
        ),
      
      "IL1B-regulated genes: TAK1 KO and rescue @ 2 h" =
        HTML(
          "For TAK1 paper work. 1303 IL1B-regulated genes (abs(log2fold) >= 1, 
          FDR <= 0.05) under WT, sham, or TAK1KO; +- 48 h of transfection using
          300 ng/cm2 of native TAK1 mRNA or D175 catalytically inactive TAK1 mRNA."
        ),
      
      "Form-regulated genes: CRTC KOs @ 2 h" =
        HTML(
          "For Priyanka's CRTC work. 382 Formoterol-regulated genes (fold >= 1.5
          or fold <= 0.67, and FDR <= 0.05) under CRTC1, CRTC2, CRTC3, or CRTC2+3
          knockout."
        )
      
    )
  })
}

# ==== create shiny app ====
shinyApp(ui = ui, server = server)





















