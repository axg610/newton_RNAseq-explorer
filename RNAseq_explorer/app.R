library(shiny)

# Define the UI
ui <- fluidPage(
  titlePanel("Newton Lab RNAseq Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Select a dataset:", choices=c("A549 IL1B Bud", 
                                                            "HBE IL1B IFNg Dex",
                                                            "HBE TNF Form",
                                                            "ALI IL1B Bud",
                                                            "ALI Bud Form",
                                                            "ALI Cig IL17A Bud",
                                                            "ALI CMV",
                                                            "ALI HRV (Bai 2015)",
                                                            "BEAS-2B Mom Ind",
                                                            "Basal expression - all cells" ,
                                                            "A549 Ad-DUSP1 IL1B - TPM",
                                                            "A549 Ad-DUSP1 IL1B - fold",
                                                            "A549 Ad-IkBa IL1B Dex - TPM",
                                                            "A549 Ad-IkBa IL1B Dex - fold",
                                                            "BEAS-2B AKO/BKO Form Bud (2h) - TPM",
                                                            "BEAS-2B AKO/BKO Form Bud (2h) - fold",
                                                            "BEAS-2B Tap Bay Vil - TPM",
                                                            "BEAS-2B Tap Bay Vil - fold")),
      textInput("gene", "Enter a gene name:",value = "FKBP5"),
      
      helpText("Note: This page takes one gene only as an input, and the official gene symbol must be used"),
      
      actionButton("button", "Show", icon = icon("redo"))
      
    ),
    mainPanel(
      tabsetPanel(type = "tabs", 
                  tabPanel("Plot", plotOutput("expression_plot")), 
                  tabPanel("Table", tableOutput("table"))
      )
    )
  )
)


# Define server
server <- function(input, output) {
  library(tidyverse)
  my_cols = c("NS" = "grey10", "IL1B" = "#FC4E07", "Bud" = "#E7B800", "IL1B + Bud" = "#2E9FDF",
              "IL17A" = "#FC4E07", "IL17A + Bud" = "#2E9FDF",
              "TNF" = "#FC4E07", "Dex" = "#E7B800", "IL1B + Dex" = "#2E9FDF",
              "Indacaterol" = "#CC79A7", "Mometazone" = "#E7B800", "Mometazone + Indacaterol" = "#0072B2",
              "Form" = "#CC79A7", "TNF + Form" = "#0072B2", "Fsk" = "grey50", "Bud + Form" = "#0072B2",
              "naive" = "grey10", "Ad-GFP" = "#2E9FDF", "Ad-DUSP1" = "firebrick", "Ad-IkBa" = "firebrick", 
              "Control" = "#0072B2", "CMV" = "firebrick", "Vehicle" = "#0072B2", "HRV16" = "firebrick",
              "WT" = "grey50", "AKO" = "firebrick1", "BKO" = "steelblue1",
              "Tap" = "firebrick", "Bay" = "#2E9FDF", "Vil" = "#E7B800")
  
  #standard theme for ggplot figures
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
  
  
  # Generate plot
  x = eventReactive({
    input$button
    input$dataset
  }, {     
    if (input$dataset == "A549 IL1B Bud") {
      read_tsv(paste0("processed_tables/a549_ib.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene)) %>% 
        mutate(treatment = factor(treatment,  levels = c("NS", "IL1B", "Bud", "IL1B + Bud")),
               time = factor(time, levels =  c("1", "2", "6", "12", "24"))) }
    
    else if (input$dataset == "BEAS-2B Mom Ind") {
      read_tsv(paste0("processed_tables/b2b_mi.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene)) %>% 
        mutate(treatment = factor(treatment, levels = unique(treatment)),
               time = factor(time, levels = unique(time)))}
    
    else if (input$dataset == "Basal expression - all cells") {
      read_tsv(paste0("processed_tables/all_cells.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene)) %>% 
        mutate(cells = factor(cells, levels = c("A549", "B2B", "HBE", "ALI", "Brushings"))) }
    
    else if (input$dataset == "HBE IL1B IFNg Dex") {
      read_tsv(paste0("processed_tables/hbe_iid.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene)) %>% 
        mutate(treatment = factor(treatment,  levels = unique(treatment)),
               ttt = factor(ttt, levels = unique(ttt)),
               time = factor(time, levels =  c("2", "6", "24"))) }
    
    else if (input$dataset == "HBE TNF Form") {
      read_tsv(paste0("processed_tables/hbe_tf.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene)) %>% 
        mutate(treatment = factor(treatment,  levels = c("NS", "TNF", "Form", "TNF + Form", "Fsk")),
               time = factor(time, levels =  c("1", "2", "6", "18"))) }
    
    else if (input$dataset == "ALI IL1B Bud") {
      read_tsv(paste0("processed_tables/ali_ib.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene)) %>% 
        mutate(treatment = factor(treatment,  levels = unique(treatment))) }
    
    else if (input$dataset == "ALI Bud Form") {
      read_tsv(paste0("processed_tables/ali_bf.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene)) %>% 
        mutate(treatment = factor(treatment,  levels = unique(treatment)),
               time = factor(time, levels =  c("2", "6"))) }
    
    else if (input$dataset == "ALI Cig IL17A Bud") {
      read_tsv(paste0("processed_tables/ali_cig.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene)) %>% 
        mutate(treatment = factor(treatment, levels= c("NS", "IL17A", "Bud", "IL17A + Bud")),
               condition = factor(condition, levels = c("Clean Air", "Cigarette Smoke")),
               cond = factor(cond, levels= c("clean", "cig")))}
    
    else if (input$dataset == "ALI CMV") {
      read_tsv(paste0("processed_tables/ali_cmv.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene)) %>% 
        mutate(condition = factor(condition, levels = c("Control", "CMV")))}
    
    else if (input$dataset == "ALI HRV (Bai 2015)") {
      read_tsv(paste0("processed_tables/ali_hrv.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene)) %>% 
        mutate(condition = factor(condition, levels = c("Vehicle", "HRV16")),
               disease_state = factor(disease_state, levels = c("non-asthmatic", "asthmatic")))}
    
    else if (input$dataset == "A549 Ad-DUSP1 IL1B - TPM") {
      read_tsv(paste0("processed_tables/dusp.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene)) %>% 
        mutate(treatment = factor(treatment,  levels = unique(treatment)),
               condition = factor(condition, levels = unique(condition)),
               time = factor(time, levels = unique(time))) }
    
    else if (input$dataset == "A549 Ad-DUSP1 IL1B - fold") {
      read_tsv(paste0("processed_tables/dusp.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene) & treatment != "NS") %>% 
        mutate(treatment = factor(treatment,  levels = unique(treatment)),
               condition = factor(condition, levels = unique(condition)),
               time = factor(time, levels = unique(time))) }
    
    else if (input$dataset == "A549 Ad-IkBa IL1B Dex - fold") {
      read_tsv(paste0("processed_tables/ikba.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene) & treatment != "NS") %>% 
        mutate(treatment = factor(treatment,  levels = unique(treatment)),
               condition = factor(condition, levels = unique(condition))) }
    
    else if (input$dataset == "A549 Ad-IkBa IL1B Dex - TPM") {
      read_tsv(paste0("processed_tables/ikba.txt"), show_col_types = F, progress = F) %>% 
        filter(Gene == toupper(input$gene)) %>% 
        mutate(treatment = factor(treatment,  levels = unique(treatment)),
               condition = factor(condition, levels = unique(condition))) }
    
    else if (input$dataset == "BEAS-2B AKO/BKO Form Bud (2h) - fold") {
      read_tsv(paste0("processed_tables/b2b_pka.txt"), show_col_types = F, progress = F) %>% 
        filter(gene == toupper(input$gene) & treatment != "NS") %>% 
        mutate(condition = factor(condition, levels = c("WT", "AKO", "BKO")),
               treatment = factor(treatment, levels = c("NS", "Form", "Bud", "F+B"))) }
    
    else if (input$dataset == "BEAS-2B AKO/BKO Form Bud (2h) - TPM") {
      read_tsv(paste0("processed_tables/b2b_pka.txt"), show_col_types = F, progress = F) %>% 
        filter(gene == toupper(input$gene)) %>% 
        mutate(condition = factor(condition, levels = c("WT", "AKO", "BKO")),
               treatment = factor(treatment, levels = c("NS", "Form", "Bud", "F+B"))) }
    
    else if (input$dataset == "BEAS-2B Tap Bay Vil - TPM") {
      read_tsv(paste0("processed_tables/b2b_tap_bay_vil.txt"), show_col_types = F, progress = F) %>% 
        filter(gene == toupper(input$gene)) %>% 
        mutate(treatment = factor(treatment, levels = c("NS", "Vil", "Bay", "Tap"))) }
    
    else if (input$dataset == "BEAS-2B Tap Bay Vil - fold") {
      read_tsv(paste0("processed_tables/b2b_tap_bay_vil.txt"), show_col_types = F, progress = F) %>% 
        filter(gene == toupper(input$gene)) %>% 
        mutate(treatment = factor(treatment, levels = c("NS", "Vil", "Bay", "Tap"))) %>%
        filter(treatment != "NS") }
    
    else {output$text = renderPrint("select a valid dataset")} 
  })
  
  output$expression_plot <- renderPlot({
    if (input$dataset == "A549 IL1B Bud") {
      x() %>% 
        ggplot(aes(time, log2_tpm, group = treatment, color = treatment)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1)+
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), size = 0.5, show.legend = F) + 
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.5, show.legend = F, geom = "errorbar", width = 0.2) +  
        scale_color_manual(values = my_cols) +  
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        #expand_limits(y = c(0)) +
        facet_wrap(~Gene)+
        labs(y = expression("log"[2]*" TPM"), x = "time (h)", color = "Treatment") + 
        my_theme + theme(aspect.ratio = 1)
    } 
    else if (input$dataset == "BEAS-2B Mom Ind") {
      x() %>% 
        ggplot(aes(time, log2_tpm, group = treatment, color = treatment)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1)+
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), size = 0.5, show.legend = F) + 
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.5, show.legend = F, geom = "errorbar", width = 0.2) +  
        scale_color_manual(values = my_cols) +  
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        #expand_limits(y = c(0)) +
        facet_wrap(~Gene)+
        labs(y = expression("log"[2]*" TPM"), x = "time (h)", color = "Treatment") + 
        my_theme + theme(aspect.ratio = 1)
    }
    
    else if (input$dataset == "Basal expression - all cells") {
      x() %>% 
        ggplot(aes(cells, log2_tpm))+
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5)  + 
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        labs(y = expression("log"[2]*" TPM"), x = NULL) + 
        expand_limits(y = 0) +
        facet_wrap(~Gene)+
        my_theme + 
        theme(aspect.ratio = 0.8 , 
              legend.position="none", 
              axis.text.x = element_text(angle = 45, hjust = 1))
    }
    else if (input$dataset == "HBE IL1B IFNg Dex") {
      x() %>% 
        ggplot(aes(time, log2_tpm, group = ttt, color = ttt)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1)+
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), size = 0.5, show.legend = F) + 
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.5, show.legend = F, geom = "errorbar", width = 0.2) +  
        scale_color_manual(values = my_cols) +  
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        #expand_limits(y = c(0)) +
        facet_grid(Gene~ifn)+
        labs(y = expression("log"[2]*" TPM"), x = "time (h)", color = "Treatment") + 
        my_theme + theme(aspect.ratio = 1.5)
    }
    else if (input$dataset == "HBE TNF Form") {
      x() %>% 
        ggplot(aes(time, log2_tpm, group = treatment, color = treatment)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1)+
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), size = 0.5, show.legend = F) + 
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.5, show.legend = F, geom = "errorbar", width = 0.2) +  
        scale_color_manual(values = my_cols) +  
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        #expand_limits(y = c(0)) +
        facet_wrap(~Gene)+
        labs(y = expression("log"[2]*" TPM"), x = "time (h)", color = "Treatment") + 
        my_theme + theme(aspect.ratio = 1)
    }
    else if (input$dataset == "ALI IL1B Bud") {
      x() %>% 
        ggplot(aes(treatment, log2_tpm, fill = treatment))+
        geom_boxplot(coef = 3, outlier.size = 0, outlier.stroke = 0, size = 0.5, show.legend = T)  + 
        labs(y = expression("log"[2]*" TPM"), fill = "Treatment", x= NULL) + 
        scale_fill_manual(values = my_cols)+
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        facet_wrap(~Gene)+
        my_theme + theme(aspect.ratio = 1.5, axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    else if (input$dataset == "ALI Bud Form") {
      x() %>% 
        ggplot(aes(time, log2_tpm, group = treatment, color = treatment)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1)+
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), size = 0.5, show.legend = F) + 
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.5, show.legend = F, geom = "errorbar", width = 0.2) +  
        scale_color_manual(values = my_cols) +  
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        #expand_limits(y = c(0)) +
        facet_wrap(~Gene)+
        labs(y = expression("log"[2]*" TPM"), x = "time (h)", color = "Treatment") + 
        my_theme + theme(aspect.ratio = 2)
    }
    else if (input$dataset == "ALI Cig IL17A Bud") {
      x() %>% 
        ggplot(aes(treatment, log2_tpm, group = treatment, fill = treatment)) +
        geom_boxplot(aes(linetype = condition),color = "grey10", linewidth = 0.4, outlier.shape = NULL, 
                     outlier.size = 0, outlier.stroke = 0, coef = 3) +
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        scale_fill_manual(values = my_cols)+
        facet_grid(Gene~cond, scales = "free_y")+
        labs(y = expression("log"[2]*" tpm"), x = NULL, fill = NULL) + 
        my_theme + theme(aspect.ratio = 2,axis.ticks.x = element_blank(), axis.text.x = element_blank()) 
    }
    else if (input$dataset == "ALI CMV") {
      x() %>% 
        ggplot(aes(condition, log2_tpm, fill = condition)) +
        geom_boxplot(color = "grey10", linewidth = 0.4, outlier.shape = NULL, outlier.size = 0, outlier.stroke = 0, coef = 3) +
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        scale_fill_manual(values = my_cols)+
        facet_wrap(Gene~.)+
        labs(y = expression("log"[2]*" tpm"), x = NULL, fill = NULL) + 
        my_theme + theme(aspect.ratio = 2,axis.ticks.x = element_blank(), axis.text.x = element_blank()) 
    }
    else if (input$dataset == "ALI HRV (Bai 2015)") {
      x() %>% 
        ggplot(aes(condition, log2_tpm, fill = condition)) +
        geom_boxplot(aes(linetype = disease_state), color = "grey10", linewidth = 0.4, outlier.shape = NULL, 
                     outlier.size = 0, outlier.stroke = 0, coef = 3) +
        scale_y_continuous(breaks = scales::breaks_pretty(min.n=4), expand = expansion(add = 0.5)) + 
        scale_fill_manual(values = my_cols)+
        facet_grid(Gene~disease_state)+
        labs(y = expression("log"[2]*" tpm"), x = NULL, fill = NULL, linetype = NULL) + 
        my_theme + theme(aspect.ratio = 1.5, strip.text.x = element_text(size = 9),
                         axis.ticks.x = element_blank(), axis.text.x = element_blank()) 
    }
    else if (input$dataset == "A549 Ad-DUSP1 IL1B - TPM") {
      x() %>% 
        ggplot(aes(time, log2_tpm, group = interaction(treatment,condition) , color = treatment)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1, alpha = 0.8) +
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), stroke = 0.5, size = 0.4, show.legend = F, alpha = 0.8) +  
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.3, show.legend = F, geom = "errorbar", width = 0.2, alpha = 0.8) +  
        scale_color_manual(values = my_cols) + 
        facet_grid(Gene~condition)+
        labs(y = expression("log"[2]*" TPM"), x= "time (h)" , color = "Treatment") + 
        my_theme + theme(aspect.ratio = 2) 
      
    }
    else if (input$dataset == "A549 Ad-DUSP1 IL1B - fold") {
      x() %>% 
        ggplot(aes(time, log2_fold, group = condition , color = condition)) +
        geom_hline(yintercept = 0, linewidth = 0.75, alpha = 0.5, linetype = 2) +
        stat_summary(fun = mean, geom = "line", linewidth = 1, alpha = 0.8) +
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), stroke = 0.5, size = 0.4, show.legend = F, alpha = 0.8) +  
        stat_summary(fun.data = "mean_se", fun.args = list(mult = 1), linewidth = 0.3, show.legend = F, geom = "errorbar", width = 0.2, alpha = 0.8) +  
        scale_color_manual(values = my_cols) + 
        facet_wrap(~Gene)+
        labs(y = expression("log"[2]*" fold"), x= "time (h)" , color = "Treatment") + 
        my_theme + theme(aspect.ratio = 2) 
      
    }
    else if (input$dataset == "A549 Ad-IkBa IL1B Dex - fold") {
      x() %>% 
        ggplot(aes(condition, log2_fold, group = condition , fill = condition)) +
        geom_hline(yintercept = 0, linewidth = 0.75, alpha = 0.5, linetype = 2) +
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5)  + 
        scale_fill_manual(values = my_cols) + 
        facet_grid(Gene~treatment)+
        labs(y = expression("log"[2]*" fold"), x= NULL , fill = "Condition") + 
        my_theme + theme(aspect.ratio = 2, axis.ticks.x = element_blank(), axis.text.x = element_blank()) 
    }
    else if (input$dataset == "A549 Ad-IkBa IL1B Dex - TPM") {
      x() %>% 
        ggplot(aes(condition, log2_tpm, group = condition , fill = condition)) +
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5)  + 
        scale_fill_manual(values = my_cols) + 
        facet_grid(Gene~treatment)+
        labs(y = expression("log"[2]*" fold"), x= NULL , fill = "Condition") + 
        my_theme + theme(aspect.ratio = 2, axis.ticks.x = element_blank(), axis.text.x = element_blank()) 
    }
    else if (input$dataset == "BEAS-2B AKO/BKO Form Bud (2h) - fold") {
      x() %>% 
        ggplot(aes(x = condition, y = fold, fill = condition)) +
        facet_grid(gene~treatment) +
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5) +
        geom_point() +
        labs(y = expression("log"[2]*" fold"), x = NULL, fill = "Condition") +
        scale_fill_manual(values = my_cols) +
        my_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 2)
    }
    else if (input$dataset == "BEAS-2B AKO/BKO Form Bud (2h) - TPM") {
      x() %>% 
        ggplot(aes(x = condition, y = tpm, fill = condition)) +
        facet_grid(gene~treatment) +
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5) +
        geom_point() +
        labs(y = expression("log"[2]*" TPM"), x = NULL, fill = "Condition") +
        scale_fill_manual(values = my_cols) +
        my_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 2)
    }
    else if (input$dataset == "BEAS-2B Tap Bay Vil - fold") {
      x() %>% 
        ggplot(aes(x = treatment, y = fold, fill = treatment)) +
        facet_wrap(~gene) +
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5) +
        geom_point() +
        labs(y = expression("log"[2]*" fold"), x = NULL, fill = "Treatment") +
        scale_fill_manual(values = c("NS" = "grey50", "Tap" = "firebrick", "Bay" = "#2E9FDF", "Vil" = "#E7B800")) +
        theme_bw() +
        my_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 1.5)
    }
    else if (input$dataset == "BEAS-2B Tap Bay Vil - TPM") {
      x() %>% 
        ggplot(aes(x = treatment, y = tpm, fill = treatment)) +
        facet_wrap(~gene) +
        geom_boxplot(coef = 3, outlier.size = 0, linewidth = 0.5) +
        geom_point() +
        labs(y = expression("log"[2]*" TPM"), x = NULL, fill = "Treatment") +
        scale_fill_manual(values = c("NS" = "grey50", "Tap" = "firebrick", "Bay" = "#2E9FDF", "Vil" = "#E7B800")) +
        theme_bw() +
        my_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 1.5)
    }
    else {output$text = renderPrint("select a valid dataset")} 
  }, res = 120,)
  
  output$table <- renderTable({x()})
  
}

# Create Shiny app
shinyApp(ui = ui, server = server)