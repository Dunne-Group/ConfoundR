#### Load required packages ####
library(shiny)
library(shinydashboard)
library(dashboardthemes)
library(shinyFeedback)
library(shinybusy)
library(msigdbr)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(ggpubr)
library(ComplexHeatmap)
library(limma)
library(DESeq2)
library(clusterProfiler)
library(fgsea)
library(enrichplot)
library(RColorBrewer)


#### shinyServer ####
shinyServer(function(input, output, session) {
    
    #### INTRO TAB ####
    
    ## Switch tab based on buttons embedded in introduction text
    
    # If user clicks the Expression Boxplots button make the boxplots tab the active tab
    observeEvent(input$switchTab01, {
        updateTabItems(session = session,
                       inputId = "menu1",
                       selected = "boxplots")
    })
    
    # If user clicks the Expression Boxplots button make the boxplots tab the active tab
    observeEvent(input$switchTab04, {
        updateTabItems(session = session,
                       inputId = "menu1",
                       selected = "boxplots")
    })
    
    # If user clicks Heatmap button make the expression_heatmap tab the active tab
    observeEvent(input$switchTab02, {
        updateTabItems(session = session,
                       inputId = "menu1",
                       selected = "expression_heatmap")
    })
    
    # If user clicks Heatmap button make the expression_heatmap tab the active tab
    observeEvent(input$switchTab05, {
        updateTabItems(session = session,
                       inputId = "menu1",
                       selected = "expression_heatmap")
    })
    
    # If user clicks GSEA button make the GSEA tab the active tab
    observeEvent(input$switchTab03, {
        updateTabItems(session = session,
                       inputId = "menu1",
                       selected = "gsea")
    })
    
    # If user clicks GSEA button make the GSEA tab the active tab
    observeEvent(input$switchTab06, {
        updateTabItems(session = session,
                       inputId = "menu1",
                       selected = "gsea")
    })
    
    
    #### BOXPLOT TAB ####
    
    # Convert the user entered text (gene symbol) to upper case in case the user
    # entered it in lower case 
    boxplot_chosen_gene <- eventReactive(input$choose_gene, {
        # Remove any horizontal (spaces, tabs) before or after the gene symbol
        chosen_gene <- trimws(input$chosen_gene, 
                              which = "both", # leading and trailing whitespace
                              whitespace = "[ \t]" # space and tab
                              )
        # Convert to gene symbol to upper case in case the user entered it
        # in lower case or mixed case
        chosen_gene <- toupper(chosen_gene)
        # Return the uppercase gene symbol with the whitespace removed 
        chosen_gene
    })
    
    # When user clicks choose_gene button check if the gene symbol
    # input by the user is present in the GSE39396 dataset
    GSE39396_chosen_gene_present <- eventReactive(input$choose_gene, {
        if (boxplot_chosen_gene() %in% rownames(GSE39396_exp)){
            TRUE
        }
        else if (boxplot_chosen_gene() %in% rownames(GSE81838_exp)) {
            TRUE
        }
        else if (boxplot_chosen_gene() %in% rownames(GSE31279_exp)) {
            TRUE
        }
        else if (boxplot_chosen_gene() %in% rownames(GSE35602_exp)) {
            TRUE
        }
        else if (boxplot_chosen_gene() %in% rownames(GSE164665_normalised_counts)) {
            TRUE
        }
        else if (boxplot_chosen_gene() %in% rownames(GSE14548_exp)) {
            TRUE
        }
        else if (boxplot_chosen_gene() %in% rownames(GSE9899_exp)) {
            TRUE
        }
        else if (boxplot_chosen_gene() %in% rownames(GSE97284_exp)) {
            TRUE
        }
        else {
            FALSE
        }
    })
    
    # When user clicks the choose_gene button check if anything has been
    # entered into the chosen_gene input box and if nothing has been
    # entered then issue a message to user to enter a gene symbol
    observeEvent(input$choose_gene, {

        if (boxplot_chosen_gene() != "" & GSE39396_chosen_gene_present() == FALSE){
            feedbackDanger("chosen_gene",
                           boxplot_chosen_gene() == "",
                           "Please enter a gene symbol")
            feedbackWarning("chosen_gene",
                            boxplot_chosen_gene() != "" & GSE39396_chosen_gene_present() == FALSE,
                            "Chosen gene not found in any dataset. Please ensure you have typed the gene symbol correctly. Alternatively you may try an alias gene symbol.")
        }
        else if (boxplot_chosen_gene() == ""){
            feedbackWarning("chosen_gene",
                            boxplot_chosen_gene() != "" & GSE39396_chosen_gene_present() == FALSE,
                            "Chosen gene not found in any dataset. Please ensure you have typed the gene symbol correctly. Alternatively you may try an alias gene symbol.")
            feedbackDanger("chosen_gene",
                           boxplot_chosen_gene() == "",
                           "Please enter a gene symbol")
        }
        else {
            hideFeedback("chosen_gene")
        }
    })
    
    
    
    #### GSE39396 box plot ####
    
    # If the gene symbol input by the user is found then select
    # the gene expression values for this gene
    GSE39396_chosen_gene_exp <- eventReactive(input$choose_gene, {
        req(GSE39396_chosen_gene_present() == TRUE)
        as.numeric(GSE39396_exp[boxplot_chosen_gene(), ])
    })
    
    # Make a dataframe containing a column with the expression values
    # for the chosen gene and column indicating the cell type for each
    # sample in the GSE39396 dataset
    GSE39396_gg_df <- reactive({
        data.frame(gene_exp = GSE39396_chosen_gene_exp(),
                   cell_type = GSE39396_cell_types)
    })
    
    # Create a theme to use for boxplots
    boxplot_theme <- 
        theme_classic() +
        theme(
            legend.position = "none",
            plot.title = element_text(
                size = 20, family = "Helvetica", face = "bold",
                colour="black", hjust = 0.5, vjust = 1),
            axis.title.y = element_text(
                size = 18, family = "Helvetica", colour="black",
                margin = margin(r = 10)),
            axis.text.y = element_text(size=16, colour="black"),
            axis.title.x = element_text(
                size = 18, family = "Helvetica", colour="black",
                margin = margin(t = -10)), 
            axis.text.x = element_text(size=16,colour="black"),
            legend.title = element_text(size=14, face = "bold"), 
            legend.text = element_text(size=14)
        )
    
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    GSE39396_boxplot_gg_obj <- reactive({
        validate(
            need(GSE39396_chosen_gene_exp(), "Chosen gene not present in this dataset")
        )
        ggplot(GSE39396_gg_df(),
               aes(x = cell_type,
                   y = gene_exp,
                   col = cell_type)) +
            geom_boxplot(outlier.shape = NA,
                         lwd = 1) +
            geom_beeswarm(size = 2.5) +
            stat_compare_means(comparisons = list(c("Epithelial","Leukocytes" ),
                                                  c("Epithelial", "Endothelial"),
                                                  c("Epithelial", "Fibroblasts")),
                               method = "wilcox.test",
                               size = 5) +
            scale_color_manual(values = c("Epithelial" = "#99C945",
                                          "Leukocytes" = "#24796C",
                                          "Endothelial" = "#52BCA3",
                                          "Fibroblasts" = "#CC61B0")) +
            labs(x = NULL,
                 y = "Expression") +
            boxplot_theme
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    output$GSE39396_boxplot <- renderPlot({
        GSE39396_boxplot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save boxplot
    output$GSE39396_boxplot_download_button <- renderUI({
        req(GSE39396_chosen_gene_exp())
        req(GSE39396_boxplot_gg_obj())
        downloadButton("GSE39396_boxplot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
                       ) 
    })
    
    # Download handling for when download/save boxplot button clicked
    output$GSE39396_boxplot_download <- downloadHandler(
        filename = function() {
            paste(boxplot_chosen_gene(), "_GSE39396_CRC_boxplot" , ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, GSE39396_boxplot_gg_obj(), device = "png")
        }
    )
    
    
    #### GSE81838 boxplot ####
    
    # If the gene symbol input by the user is found then select
    # the gene expression values for this gene
    GSE81838_chosen_gene_exp <- eventReactive(input$choose_gene, {
        req(GSE39396_chosen_gene_present() == TRUE)
        as.numeric(GSE81838_exp[boxplot_chosen_gene(), ])
    })
    
    # Make a dataframe containing a column with the expression values
    # for the chosen gene and column indicating the cell type for each
    # sample in the GSE81838 dataset
    GSE81838_gg_df <- reactive({
        data.frame(gene_exp = GSE81838_chosen_gene_exp(),
                   cell_type = GSE81838_cell_types)
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    GSE81838_boxplot_gg_obj <- reactive({
        validate(
            need(GSE81838_chosen_gene_exp(), "Chosen gene not present in this dataset")
        )
        ggplot(GSE81838_gg_df(),
               aes(x = cell_type,
                   y = gene_exp,
                   col = cell_type)) +
            geom_boxplot(outlier.shape = NA,
                         lwd = 1) +
            geom_beeswarm(size = 2.5) +
            stat_compare_means(comparisons = list(c("Epithelium", "Stroma")),
                               method = "wilcox.test",
                               size = 5) +
            scale_color_manual(values = c("Epithelium" = "#5D69B1",
                                          "Stroma" = "#E58606")) +
            labs(x = NULL,
                 y = "Expression") +
            boxplot_theme
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    output$GSE81838_boxplot <- renderPlot({
        GSE81838_boxplot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save boxplot
    output$GSE81838_boxplot_download_button <- renderUI({
        req(GSE81838_chosen_gene_exp())
        req(GSE81838_boxplot_gg_obj())
        downloadButton("GSE81838_boxplot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save boxplot button clicked
    output$GSE81838_boxplot_download <- downloadHandler(
        filename = function() {
            paste(boxplot_chosen_gene(), "_GSE81838_TNBC_boxplot" , ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, GSE81838_boxplot_gg_obj(), device = "png")
        }
    )
    
    
    #### GSE31279 boxplot ####
    
    # If the gene symbol input by the user is found then select
    # the gene expression values for this gene
    GSE31279_chosen_gene_exp <- eventReactive(input$choose_gene, {
        req(GSE39396_chosen_gene_present() == TRUE)
        as.numeric(GSE31279_exp[boxplot_chosen_gene(), ])
    })
    
    # Make a dataframe containing a column with the expression values
    # for the chosen gene and column indicating the cell type for each
    # sample in the GSE31279 dataset
    GSE31279_gg_df <- reactive({
        data.frame(gene_exp = GSE31279_chosen_gene_exp(),
                   cell_type = GSE31279_cell_types)
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    GSE31279_boxplot_gg_obj <- reactive({
        validate(
            need(GSE31279_chosen_gene_exp(), "Chosen gene not present in this dataset")
        )
        ggplot(GSE31279_gg_df(),
               aes(x = cell_type,
                   y = gene_exp,
                   col = cell_type)) +
            geom_boxplot(outlier.shape = NA,
                         lwd = 1) +
            geom_beeswarm(size = 2.5) +
            stat_compare_means(comparisons = list(c("Epithelium", "Stroma")),
                               method = "wilcox.test",
                               size = 5) +
            scale_color_manual(values = c("Epithelium" = "#5D69B1",
                                          "Stroma" = "#E58606")) +
            labs(x = NULL,
                 y = "Expression") +
            boxplot_theme
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    output$GSE31279_boxplot <- renderPlot({
        GSE31279_boxplot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save boxplot
    output$GSE31279_boxplot_download_button <- renderUI({
        req(GSE31279_chosen_gene_exp())
        req(GSE31279_boxplot_gg_obj())
        downloadButton("GSE31279_boxplot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save boxplot button clicked
    output$GSE31279_boxplot_download <- downloadHandler(
        filename = function() {
            paste(boxplot_chosen_gene(), "_GSE31279_CRC_boxplot" , ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, GSE31279_boxplot_gg_obj(), device = "png")
        }
    )
    
    
    #### GSE35602 boxplot ####
    
    # If the gene symbol input by the user is found then select
    # the gene expression values for this gene
    GSE35602_chosen_gene_exp <- eventReactive(input$choose_gene, {
        req(GSE39396_chosen_gene_present() == TRUE)
        as.numeric(GSE35602_exp[boxplot_chosen_gene(), ])
    })
    
    # Make a dataframe containing a column with the expression values
    # for the chosen gene and column indicating the cell type for each
    # sample in the GSE35602 dataset
    GSE35602_gg_df <- reactive({
        data.frame(gene_exp = GSE35602_chosen_gene_exp(),
                   cell_type = GSE35602_cell_types)
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    GSE35602_boxplot_gg_obj <- reactive({
        validate(
            need(GSE35602_chosen_gene_exp(), "Chosen gene not present in this dataset")
        )
        ggplot(GSE35602_gg_df(),
               aes(x = cell_type,
                   y = gene_exp,
                   col = cell_type)) +
            geom_boxplot(outlier.shape = NA,
                         lwd = 1) +
            geom_beeswarm(size = 2.5) +
            stat_compare_means(comparisons = list(c("Epithelium", "Stroma")),
                               method = "wilcox.test",
                               size = 5) +
            scale_color_manual(values = c("Epithelium" = "#5D69B1",
                                          "Stroma" = "#E58606")) +
            labs(x = NULL,
                 y = "Expression") +
            boxplot_theme
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    output$GSE35602_boxplot <- renderPlot({
        GSE35602_boxplot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save boxplot
    output$GSE35602_boxplot_download_button <- renderUI({
        req(GSE35602_chosen_gene_exp())
        req(GSE35602_boxplot_gg_obj())
        downloadButton("GSE35602_boxplot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save boxplot button clicked
    output$GSE35602_boxplot_download <- downloadHandler(
        filename = function() {
            paste(boxplot_chosen_gene(), "_GSE35602_CRC_boxplot" , ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, GSE35602_boxplot_gg_obj(), device = "png")
        }
    )
    
    
    #### GSE164665 boxplot ####
    
    # If the gene symbol input by the user is found then select
    # the gene expression values for this gene
    GSE164665_chosen_gene_exp <- eventReactive(input$choose_gene, {
        req(GSE39396_chosen_gene_present() == TRUE)
        as.numeric(GSE164665_normalised_counts[boxplot_chosen_gene(), ])
    })
    
    # Make a dataframe containing a column with the expression values
    # for the chosen gene and column indicating the cell type for each
    # sample in the GSE164665 dataset
    GSE164665_gg_df <- reactive({
        data.frame(gene_exp = GSE164665_chosen_gene_exp(),
                   cell_type = GSE164665_cell_types)
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    GSE164665_boxplot_gg_obj <- reactive({
        validate(
            need(GSE164665_chosen_gene_exp(), "Chosen gene not present in this dataset")
        )
        ggplot(GSE164665_gg_df(),
               aes(x = cell_type,
                   y = gene_exp,
                   col = cell_type)) +
            geom_boxplot(outlier.shape = NA,
                         lwd = 1) +
            geom_beeswarm(size = 2.5) +
            stat_compare_means(comparisons = list(c("Epithelium", "Stroma")),
                               method = "wilcox.test",
                               size = 5) +
            scale_color_manual(values = c("Epithelium" = "#5D69B1",
                                          "Stroma" = "#E58606")) +
            labs(x = NULL,
                 y = "Expression (normalized counts)") +
            boxplot_theme
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    output$GSE164665_boxplot <- renderPlot({
        GSE164665_boxplot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save boxplot
    output$GSE164665_boxplot_download_button <- renderUI({
        req(GSE164665_chosen_gene_exp())
        req(GSE164665_boxplot_gg_obj())
        downloadButton("GSE164665_boxplot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save boxplot button clicked
    output$GSE164665_boxplot_download <- downloadHandler(
        filename = function() {
            paste(boxplot_chosen_gene(), "_GSE164665_PDAC_boxplot" , ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, GSE164665_boxplot_gg_obj(), device = "png")
        }
    )
    
    
    
    #### GSE14548 boxplot ####
    
    # If the gene symbol input by the user is found then select
    # the gene expression values for this gene
    GSE14548_chosen_gene_exp <- eventReactive(input$choose_gene, {
        req(GSE39396_chosen_gene_present() == TRUE)
        as.numeric(GSE14548_exp[boxplot_chosen_gene(), ])
    })
    
    # Make a dataframe containing a column with the expression values
    # for the chosen gene and column indicating the cell type for each
    # sample in the GSE14548 dataset
    GSE14548_gg_df <- reactive({
        data.frame(gene_exp = GSE14548_chosen_gene_exp(),
                   cell_type = GSE14548_cell_types)
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    GSE14548_boxplot_gg_obj <- reactive({
        validate(
            need(GSE14548_chosen_gene_exp(), "Chosen gene not present in this dataset")
        )
        ggplot(GSE14548_gg_df(),
               aes(x = cell_type,
                   y = gene_exp,
                   col = cell_type)) +
            geom_boxplot(outlier.shape = NA,
                         lwd = 1) +
            geom_beeswarm(size = 2.5) +
            stat_compare_means(comparisons = list(c("Epithelium", "Stroma")),
                               method = "wilcox.test",
                               size = 5) +
            scale_color_manual(values = c("Epithelium" = "#5D69B1",
                                          "Stroma" = "#E58606")) +
            labs(x = NULL,
                 y = "Expression") +
            boxplot_theme
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    output$GSE14548_boxplot <- renderPlot({
        GSE14548_boxplot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save boxplot
    output$GSE14548_boxplot_download_button <- renderUI({
        req(GSE14548_chosen_gene_exp())
        req(GSE14548_boxplot_gg_obj())
        downloadButton("GSE14548_boxplot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save boxplot button clicked
    output$GSE14548_boxplot_download <- downloadHandler(
        filename = function() {
            paste(boxplot_chosen_gene(), "_GSE14548_Breast_Cancer_boxplot" , ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, GSE14548_boxplot_gg_obj(), device = "png")
        }
    )
    
    
    #### GSE9899 boxplot ####
    
    # If the gene symbol input by the user is found then select
    # the gene expression values for this gene
    GSE9899_chosen_gene_exp <- eventReactive(input$choose_gene, {
        req(GSE39396_chosen_gene_present() == TRUE)
        as.numeric(GSE9899_exp[boxplot_chosen_gene(), ])
    })
    
    # Make a dataframe containing a column with the expression values
    # for the chosen gene and column indicating the cell type for each
    # sample in the GSE9899 dataset
    GSE9899_gg_df <- reactive({
        data.frame(gene_exp = GSE9899_chosen_gene_exp(),
                   cell_type = GSE9899_cell_types)
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    GSE9899_boxplot_gg_obj <- reactive({
        validate(
            need(GSE9899_chosen_gene_exp(), "Chosen gene not present in this dataset")
        )
        ggplot(GSE9899_gg_df(),
               aes(x = cell_type,
                   y = gene_exp,
                   col = cell_type)) +
            geom_boxplot(outlier.shape = NA,
                         lwd = 1) +
            geom_beeswarm(size = 2.5) +
            stat_compare_means(comparisons = list(c("Epithelium", "Stroma")),
                               method = "wilcox.test",
                               size = 5) +
            scale_color_manual(values = c("Epithelium" = "#5D69B1",
                                          "Stroma" = "#E58606")) +
            labs(x = NULL,
                 y = "Expression") +
            boxplot_theme
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    output$GSE9899_boxplot <- renderPlot({
        GSE9899_boxplot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save boxplot
    output$GSE9899_boxplot_download_button <- renderUI({
        req(GSE9899_chosen_gene_exp())
        req(GSE9899_boxplot_gg_obj())
        downloadButton("GSE9899_boxplot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save boxplot button clicked
    output$GSE9899_boxplot_download <- downloadHandler(
        filename = function() {
            paste(boxplot_chosen_gene(), "_GSE9899_Ovarian_Cancer_boxplot" , ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, GSE9899_boxplot_gg_obj(), device = "png")
        }
    )
    
    
    #### GSE97284 boxplot ####
    
    # If the gene symbol input by the user is found then select
    # the gene expression values for this gene
    GSE97284_chosen_gene_exp <- eventReactive(input$choose_gene, {
        req(GSE39396_chosen_gene_present() == TRUE)
        as.numeric(GSE97284_exp[boxplot_chosen_gene(), ])
    })
    
    # Make a dataframe containing a column with the expression values
    # for the chosen gene and column indicating the cell type for each
    # sample in the GSE97284 dataset
    GSE97284_gg_df <- reactive({
        data.frame(gene_exp = GSE97284_chosen_gene_exp(),
                   cell_type = GSE97284_cell_types)
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    GSE97284_boxplot_gg_obj <- reactive({
        validate(
            need(GSE97284_chosen_gene_exp(), "Chosen gene not present in this dataset")
        )
        ggplot(GSE97284_gg_df(),
               aes(x = cell_type,
                   y = gene_exp,
                   col = cell_type)) +
            geom_boxplot(outlier.shape = NA,
                         lwd = 1) +
            geom_beeswarm(size = 2.5) +
            stat_compare_means(comparisons = list(c("Epithelium", "Stroma")),
                               method = "wilcox.test",
                               size = 5) +
            scale_color_manual(values = c("Epithelium" = "#5D69B1",
                                          "Stroma" = "#E58606")) +
            labs(x = NULL,
                 y = "Expression") +
            boxplot_theme
    })
    
    # Output the plot with boxplots of the expression of the chosen
    # gene in each of the cell types
    output$GSE97284_boxplot <- renderPlot({
        GSE97284_boxplot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save boxplot
    output$GSE97284_boxplot_download_button <- renderUI({
        req(GSE97284_chosen_gene_exp())
        req(GSE97284_boxplot_gg_obj())
        downloadButton("GSE97284_boxplot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save boxplot button clicked
    output$GSE97284_boxplot_download <- downloadHandler(
        filename = function() {
            paste(boxplot_chosen_gene(), "_GSE97284_Prostate_Cancer_boxplot" , ".png", sep = "")
        },
        content = function(file) {
            ggsave(file, GSE97284_boxplot_gg_obj(), device = "png")
        }
    )
    
    
    
    #### HEATMAP TAB ####
    
    # When user clicks the gene_list_submit button on Heatmap tab
    # split the user input at each new line character to obtain
    # a vector of gene symbols
    gene_list <- eventReactive(input$gene_list_submit, {
        genes_vector <- strsplit(input$gene_list, "\n")[[1]]
        if (all(genes_vector == "")){
            ""
        }
        else{
            # Keep only unique gene symbols
            genes_vector <- unique(genes_vector)
            # Remove any empty strings (which have come from blank lines in
            # input box)
            genes_vector <- genes_vector[genes_vector != ""]
            # Remove any horizontal whitespace (space, tab) at start or end
            # i.e. before or after gene symbol
            genes_vector <- trimws(genes_vector, 
                                   which = "both", # leading and trailing whitespace
                                   whitespace = "[ \t]" # space and tab
                                   )
            # Convert gene symbols to upper case in case the user entered them
            # partially or completely in lower case
            genes_vector <- toupper(genes_vector)
        }
    })
    
    
    # Create logical vectors indicate if genes from list are present/missing in each dataset
    GSE39396_gene_list_presence <- eventReactive(input$gene_list_submit, {
        gene_list() %in% rownames(GSE39396_exp)
    })
    
    GSE81838_gene_list_presence <- eventReactive(input$gene_list_submit, {
        gene_list() %in% rownames(GSE81838_exp)
    })
    
    GSE31279_gene_list_presence <- eventReactive(input$gene_list_submit, {
        gene_list() %in% rownames(GSE31279_exp)
    })
    
    GSE35602_gene_list_presence <- eventReactive(input$gene_list_submit, {
        gene_list() %in% rownames(GSE35602_exp)
    })
    
    GSE164665_gene_list_presence <- eventReactive(input$gene_list_submit, {
        gene_list() %in% rownames(GSE164665_vst_counts)
    })
    
    GSE14548_gene_list_presence <- eventReactive(input$gene_list_submit, {
        gene_list() %in% rownames(GSE14548_exp)
    })
    
    GSE9899_gene_list_presence <- eventReactive(input$gene_list_submit, {
        gene_list() %in% rownames(GSE9899_exp)
    })
    
    GSE97284_gene_list_presence <- eventReactive(input$gene_list_submit, {
        gene_list() %in% rownames(GSE97284_exp)
    })
    
    
    # Create vectors with genes from list missing in each dataset
    GSE39396_gene_list_missing_genes <- eventReactive(input$gene_list_submit, {
        gene_list()[!GSE39396_gene_list_presence()]
    })
    
    GSE81838_gene_list_missing_genes <- eventReactive(input$gene_list_submit, {
        gene_list()[!GSE81838_gene_list_presence()]
    })
    
    GSE31279_gene_list_missing_genes <- eventReactive(input$gene_list_submit, {
        gene_list()[!GSE31279_gene_list_presence()]
    })
    
    GSE35602_gene_list_missing_genes <- eventReactive(input$gene_list_submit, {
        gene_list()[!GSE35602_gene_list_presence()]
    })
    
    GSE164665_gene_list_missing_genes <- eventReactive(input$gene_list_submit, {
        gene_list()[!GSE164665_gene_list_presence()]
    })
    
    GSE14548_gene_list_missing_genes <- eventReactive(input$gene_list_submit, {
        gene_list()[!GSE14548_gene_list_presence()]
    })
    
    GSE9899_gene_list_missing_genes <- eventReactive(input$gene_list_submit, {
        gene_list()[!GSE9899_gene_list_presence()]
    })
    
    GSE97284_gene_list_missing_genes <- eventReactive(input$gene_list_submit, {
        gene_list()[!GSE97284_gene_list_presence()]
    })
    
    
    # Create vector of genes from user's list of genes that are not present in every dataset
    gene_list_missing_genes_all_datasets <- eventReactive(input$gene_list_submit, {
        Reduce(intersect,
               list(GSE39396_gene_list_missing_genes(),
                     GSE81838_gene_list_missing_genes(),
                     GSE31279_gene_list_missing_genes(),
                     GSE35602_gene_list_missing_genes(),
                     GSE164665_gene_list_missing_genes(),
                     GSE14548_gene_list_missing_genes(),
                     GSE9899_gene_list_missing_genes(),
                     GSE97284_gene_list_missing_genes()
               ))
    })
    
    # Logical variable to indicate if there genes in the user's list of heatmap genes
    # that are missing in every dataset - controls conditionalPanel in UI
    output$heatmap_genes_missing_all_datasets <- eventReactive(input$gene_list_submit, {
        if (gene_list_present() == FALSE){
            return(FALSE)
        }
        else if (is.null(gene_list_missing_genes_all_datasets())){
            return(FALSE)
        }
        else if (!is.null(gene_list_missing_genes_all_datasets()) & length(gene_list_missing_genes_all_datasets()) > 0){
            return(TRUE)
        }
        else {
            return(FALSE)
        }
    })
    
    # output$heatmap_genes_missing_all_datasets will be used in the UI code to control
    # the appearance of a message to users using a conditional panel
    # therefore its state must be passed to the UI even though its value 
    # is never displayed by the UI
    outputOptions(output, "heatmap_genes_missing_all_datasets", suspendWhenHidden = FALSE)
    
    # List of genes to output in modal in UI as missing from all datasets
    output$heatmap_genes_missing_all_datasets_modal_content <- renderUI({
        missing_all <- gene_list_missing_genes_all_datasets()
        HTML(paste(missing_all, collapse = "<br/>"))
    })
    
    # Display modal with missing genes if user clicks the actionLink in UI
    observeEvent(input$heatmap_genes_missing_all_datasets, {
        showModal(
            modalDialog(
                htmlOutput("heatmap_genes_missing_all_datasets_modal_content"),
                footer = modalButton("Close", icon = icon("remove", lib = "glyphicon")),
                easyClose = TRUE,
                size = "s"
            )
        )
    })
    
    
    # Create vector of genes from user's list of genes that are not present in at least one dataset
    gene_list_missing_genes_any_dataset <- eventReactive(input$gene_list_submit, {
        unique(c(GSE39396_gene_list_missing_genes(),
                 GSE81838_gene_list_missing_genes(),
                 GSE31279_gene_list_missing_genes(),
                 GSE35602_gene_list_missing_genes(),
                 GSE164665_gene_list_missing_genes(),
                 GSE14548_gene_list_missing_genes(),
                 GSE9899_gene_list_missing_genes(),
                 GSE97284_gene_list_missing_genes()
        ))
    })
    
    # Logical variable to indicate if there genes in the user's list of heatmap genes
    # that are missing in at least one dataset - controls conditionalPanel in UI
    output$heatmap_genes_missing_any_dataset <- eventReactive(input$gene_list_submit, {
        if (gene_list_present() == FALSE){
            return(FALSE)
        }
        else if (is.null(gene_list_missing_genes_any_dataset())){
            return(FALSE)
        }
        else if (!is.null(gene_list_missing_genes_any_dataset()) & !is.null(gene_list_missing_genes_all_datasets())){
            missing_any_not_all <- setdiff(gene_list_missing_genes_any_dataset(), gene_list_missing_genes_all_datasets())
            if (!is.null(missing_any_not_all) & length(missing_any_not_all) > 0){
                return(TRUE)
            }
            else {
                return(FALSE)
            }
        }
        else {
            return(FALSE)
        }
    })
    
    # output$heatmap_genes_missing_any_dataset will be used in the UI code to control
    # the appearance of a message to users using a conditional panel
    # therefore its state must be passed to the UI even though its value 
    # is never displayed by the UI
    outputOptions(output, "heatmap_genes_missing_any_dataset", suspendWhenHidden = FALSE)
    
    # List of genes to output in modal in UI as missing from at least one dataset
    output$heatmap_genes_missing_any_dataset_modal_content <- renderUI({
        missing_any <- setdiff(gene_list_missing_genes_any_dataset(), gene_list_missing_genes_all_datasets())
        HTML(paste(missing_any, collapse = "<br/>"))
    })
    
    # Display modal with genes missing from at least one dataset (but not all)
    # if user clicks the actionLink in UI
    observeEvent(input$heatmap_genes_missing_any_dataset, {
        showModal(
            modalDialog(
                htmlOutput("heatmap_genes_missing_any_dataset_modal_content"),
                footer = modalButton("Close", icon = icon("remove", lib = "glyphicon")),
                easyClose = TRUE,
                size = "s"
            )
        )
    })
    
    # When user clicks the gene_list_submit button on Heatmap tab
    # check if at least one of the genes in the list of genes
    # entered by the user is present in any dataset
    GSE39396_gene_list_present <- eventReactive(input$gene_list_submit, {
        any(gene_list() %in% rownames(GSE39396_exp))
        })
    
    GSE81838_gene_list_present <- eventReactive(input$gene_list_submit, {
        any(gene_list() %in% rownames(GSE81838_exp))
    })
    
    GSE31279_gene_list_present <- eventReactive(input$gene_list_submit, {
        any(gene_list() %in% rownames(GSE31279_exp))
    })
    
    GSE35602_gene_list_present <- eventReactive(input$gene_list_submit, {
        any(gene_list() %in% rownames(GSE35602_exp))
    })
    
    GSE164665_gene_list_present <- eventReactive(input$gene_list_submit, {
        any(gene_list() %in% rownames(GSE164665_vst_counts))
    })
    
    GSE14548_gene_list_present <- eventReactive(input$gene_list_submit, {
        any(gene_list() %in% rownames(GSE14548_exp))
    })
    
    GSE9899_gene_list_present <- eventReactive(input$gene_list_submit, {
        any(gene_list() %in% rownames(GSE9899_exp))
    })
    
    GSE97284_gene_list_present <- eventReactive(input$gene_list_submit, {
        any(gene_list() %in% rownames(GSE97284_exp))
    })
    
    gene_list_present <- eventReactive(input$gene_list_submit, {
        any(GSE39396_gene_list_present(),
            GSE81838_gene_list_present(),
            GSE31279_gene_list_present(),
            GSE35602_gene_list_present(),
            GSE164665_gene_list_present(),
            GSE14548_gene_list_present(),
            GSE9899_gene_list_present(),
            GSE97284_gene_list_present())
    })
    
    
    
    # When user clicks the choose_gene button check if that user has input some
    # text and also check that the text contains at least one gene symbol
    # which is present in at least one of the datasets
    # If these conditions aren't met issue messages to user to alert user to
    # these problems
    observeEvent(input$gene_list_submit, {
        hideFeedback("gene_list")
        if (length(gene_list()) == 1){
            # If nothing entered into input box or only new line characters
            if (gene_list() == ""){
                # This warning will not be triggered but is included to
                # clear/hide this warning if it is already present
                # feedbackWarning("gene_list",
                #                 input$gene_list != "" & gene_list_present() == FALSE,
                #                 "None of the chosen genes are present in any dataset. Please ensure you have typed the gene symbol(s) correctly.")
                # This warning will be issued as the user has not entered
                # any text/gene symbol (only possible text is a sting of new
                # line characters)
                feedbackDanger("gene_list",
                               gene_list() == "",
                               "Please enter a gene symbol(s)")
            }
            # If the user has input text/gene symbols but none of the text/gene
            # symbols match a gene symbol in any of the datasets
            else if (gene_list_present() == FALSE){
                # This warning will be issued as the user 
                feedbackWarning("gene_list",
                                input$gene_list != "" & gene_list_present() == FALSE,
                                "None of the chosen genes are present in any dataset. Please ensure you have typed the gene symbol(s) correctly. Alternatively you may try an alias for the gene symbol(s).")
            }
            else {
                hideFeedback("gene_list")
            }
        }
        else {
            if (gene_list_present() == FALSE){
                feedbackWarning("gene_list",
                                input$gene_list != "" & gene_list_present() == FALSE,
                                "None of the chosen genes are present in any dataset. Please ensure you have typed the gene symbol(s) correctly. Alternatively you may try an alias for the gene symbol(s).")
            }
            else {
                hideFeedback("gene_list")
            }
        }
    })

    
    # Check if any gene(s) in the user input list (for heatmap module)
    # is/are missing in any dataset
    output$heatmap_genes_missing <- eventReactive(input$gene_list_submit, {
        req(gene_list())
        req(gene_list_present())
        if (GSE39396_gene_list_present() & !all(gene_list() %in% rownames(GSE39396_exp))){
            return(TRUE)
        }
        else if (GSE81838_gene_list_present() & !all(gene_list() %in% rownames(GSE81838_exp))){
            return(TRUE)
        }
        else if (GSE31279_gene_list_present() & !all(gene_list() %in% rownames(GSE31279_exp))){
            return(TRUE)
        }
        else if (GSE35602_gene_list_present() & !all(gene_list() %in% rownames(GSE35602_exp))){
            return(TRUE)
        }
        else if (GSE164665_gene_list_present() & !all(gene_list() %in% rownames(GSE164665_vst_counts))){
            return(TRUE)
        }
        else if (GSE14548_gene_list_present() & !all(gene_list() %in% rownames(GSE14548_exp))){
            return(TRUE)
        }
        else if (GSE9899_gene_list_present() & !all(gene_list() %in% rownames(GSE9899_exp))){
            return(TRUE)
        }
        else if (GSE9899_gene_list_present() & !all(gene_list() %in% rownames(GSE97284_exp))){
            return(TRUE)
        }
        else{
            return(FALSE)
        }
    })

    #output$heatmap_genes_missing_test_text <- renderText({heatmap_genes_missing()})
    #output$heatmap_genes_missing_test <- reactive({heatmap_genes_missing()})
    
    # output$heatmap_genes_missing will be used in the UI code to control
    # appearance of messages to users using conditional Panels
    # therefore its state must be passed to the UI even though its value 
    # is never displayed by the UI
    outputOptions(output, "heatmap_genes_missing", suspendWhenHidden = FALSE)
    
    #### GSE39396 heatmap ####
    
    # If at least one gene in the user's input list is present
    # then select the rows of the dataframe which correspond to genes
    # present in the user's input gene list
    gene_list_df <- reactive({
        req(GSE39396_gene_list_present() == TRUE)
        found_genes <- gene_list()[gene_list() %in% rownames(GSE39396_exp)]
        GSE39396_exp[found_genes,]
    })
    
    # Z-score scale the expression of the genes present in the
    # user's input gene list
    gene_list_df_scaled <- reactive({
        req(GSE39396_gene_list_present() == TRUE)
        t(scale(t(gene_list_df())))
    })
    
    #output$genes <- renderPrint({gene_list_df()})
    
    # Crreate a heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    GSE39396_heatmap <- reactive({
        req(gene_list())
        validate(
            need(GSE39396_gene_list_present(),"Chosen gene(s) not found in this dataset")
        )
        top_annot <- HeatmapAnnotation(
            "Cell Type" = anno_block(
                gp = gpar(fill = c(
                    "Epithelial" = "#99C945",
                    "Leukocytes" = "#24796C",
                    "Endothelial" = "#52BCA3",
                    "Fibroblasts" = "#CC61B0"
                )),
                labels = c("Epithelial",
                           "Leukocytes",
                           "Endothelial",
                           "Fibroblasts"),
                labels_gp = gpar(col = "white",
                                 fontface = "bold")
            )
        )
        
        # Determine the font size to use to label genes (rows) when
        # plotting the heatmap
        row_label_size <- 400 / (nrow(gene_list_df()) * 1.5)
        # Make max possible font size 10
        if (row_label_size > 10){
            row_label_size <- 10
        }
        # Make minimum possible font size 1
        else if (row_label_size < 1){
            row_label_size <- 1
        }
        
        # Create heatmap of Z-score scale gene expression
        hm <- Heatmap(gene_list_df_scaled(),
                      name = "Z-score",
                      top_annotation = top_annot,
                      column_split = GSE39396_cell_types,
                      column_title = NULL,
                      show_column_names = FALSE,
                      row_names_gp = gpar(fontsize = row_label_size),
                      cluster_rows = FALSE,
                      clustering_method_rows = "ward.D2",
                      cluster_columns = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_dend = FALSE,
                      heatmap_legend_param = list(
                          title_gp = gpar(fontsize = 16,
                                          fontface = "bold"),
                          labels_gp = gpar(fontsize = 14)
                      )
        )
        
        hm
    })
    
    # Output the heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    output$GSE39396_gene_list_heatmap <- renderPlot({
        req(GSE39396_heatmap())
        draw(GSE39396_heatmap(),
             merge_legend = TRUE) # Put Z-score scale and cell type legend in same column
    })
    
    # If the heatmap is successfully made display button to download/save heatmap
    output$GSE39396_heatmap_download_button <- renderUI({
        req(gene_list())
        req(gene_list_df_scaled())
        req(GSE39396_heatmap())
        downloadButton("GSE39396_heatmap_download",
                       "Save Plot",
                       #class = "rightAlignButton"
                       ) 
    })
    
    # Download handling for when download/save heatmap button clicked
    output$GSE39396_heatmap_download <- downloadHandler(
        filename = function() {
            paste("GSE39396_CRC_heatmap", ".png", sep = "")
        },
        content = function(file) {
            png(file, width = 7, height = 7, units = "in", res = 300)
            draw(GSE39396_heatmap(),
                 merge_legend = TRUE)
            dev.off()
        }
    )
    
    
    #### GSE81838 heatmap ####
    
    # If at least one gene in the user's input list is present
    # then select the rows of the dataframe which correspond to genes
    # present in the user's input gene list
    GSE81838_gene_list_df <- reactive({
        req(gene_list())
        req(GSE81838_gene_list_present() == TRUE)
        found_genes <- gene_list()[gene_list() %in% rownames(GSE81838_exp)]
        GSE81838_exp[found_genes,]
    })
    
    # Z-score scale the expression of the genes present in the
    # user's input gene list
    GSE81838_gene_list_df_scaled <- reactive({
        req(GSE81838_gene_list_present() == TRUE)
        t(scale(t(GSE81838_gene_list_df())))
    })
    
    
    # Create a heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    GSE81838_heatmap <- reactive({
        req(gene_list())
        validate(
            need(GSE81838_gene_list_present(), "Chosen gene(s) not found in this dataset")
        )
        top_annot <- HeatmapAnnotation(
            "Cell Type" = anno_block(
                gp = gpar(fill = c(
                    "Epithelium" = "#5D69B1",
                    "Stroma" = "#E58606"
                )),
                labels = c("Epithelium", "Stroma"),
                labels_gp = gpar(col = "white",
                                 fontface = "bold")
            )
        )
        
        # Determine the font size to use to label genes (rows) when
        # plotting the heatmap
        row_label_size <- 400 / (nrow(GSE81838_gene_list_df()) * 1.5)
        # Make max possible font size 10
        if (row_label_size > 10){
            row_label_size <- 10
        }
        # Make minimum possible font size 1
        else if (row_label_size < 1){
            row_label_size <- 1
        }
        
        # Create heatmap of Z-score scale gene expression
        hm <- Heatmap(GSE81838_gene_list_df_scaled(),
                      name = "Z-score",
                      top_annotation = top_annot,
                      column_split = GSE81838_cell_types,
                      column_title = NULL,
                      show_column_names = FALSE,
                      row_names_gp = gpar(fontsize = row_label_size),
                      cluster_rows = FALSE,
                      clustering_method_rows = "ward.D2",
                      cluster_columns = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_dend = FALSE,
                      heatmap_legend_param = list(
                          title_gp = gpar(fontsize = 16,
                                          fontface = "bold"),
                          labels_gp = gpar(fontsize = 14)
                      )
        )
        
        hm
    })
    
    # Output the heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    output$GSE81838_gene_list_heatmap <- renderPlot({
        req(GSE81838_heatmap())
        draw(GSE81838_heatmap(),
             merge_legend = TRUE) # Put Z-score scale and cell type legend in same column
    })
    
    # If the heatmap is successfully made display button to download/save heatmap
    output$GSE81838_heatmap_download_button <- renderUI({
        req(gene_list())
        req(GSE81838_gene_list_df_scaled())
        req(GSE81838_heatmap())
        downloadButton("GSE81838_heatmap_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save heatmap button clicked
    output$GSE81838_heatmap_download <- downloadHandler(
        filename = function() {
            paste("GSE81838_TNBC_heatmap", ".png", sep = "")
        },
        content = function(file) {
            png(file, width = 7, height = 7, units = "in", res = 300)
            draw(GSE81838_heatmap(),
                 merge_legend = TRUE)
            dev.off()
        }
    )
    
    ##### GSE31279 heatmap #####
    
    # If at least one gene in the user's input list is present
    # then select the rows of the dataframe which correspond to genes
    # present in the user's input gene list
    GSE31279_gene_list_df <- reactive({
        req(GSE31279_gene_list_present() == TRUE)
        found_genes <- gene_list()[gene_list() %in% rownames(GSE31279_exp)]
        GSE31279_exp[found_genes,]
    })
    
    # Z-score scale the expression of the genes present in the
    # user's input gene list
    GSE31279_gene_list_df_scaled <- reactive({
        req(GSE31279_gene_list_present() == TRUE)
        t(scale(t(GSE31279_gene_list_df())))
    })
    
    # Create a heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    GSE31279_heatmap <- reactive({
        req(gene_list())
        validate(
            need(GSE31279_gene_list_present(), "Chosen gene(s) not found in this dataset")
        )
        top_annot <- HeatmapAnnotation(
            "Cell Type" = anno_block(
                gp = gpar(fill = c(
                    "Epithelium" = "#5D69B1",
                    "Stroma" = "#E58606"
                )),
                labels = c("Epithelium", "Stroma"),
                labels_gp = gpar(col = "white",
                                 fontface = "bold")
            )
        )
        
        # Determine the font size to use to label genes (rows) when
        # plotting the heatmap
        row_label_size <- 400 / (nrow(GSE31279_gene_list_df()) * 1.5)
        # Make max possible font size 10
        if (row_label_size > 10){
            row_label_size <- 10
        }
        # Make minimum possible font size 1
        else if (row_label_size < 1){
            row_label_size <- 1
        }
        
        # Create heatmap of Z-score scale gene expression
        hm <- Heatmap(GSE31279_gene_list_df_scaled(),
                      name = "Z-score",
                      top_annotation = top_annot,
                      column_split = GSE31279_cell_types,
                      column_title = NULL,
                      show_column_names = FALSE,
                      row_names_gp = gpar(fontsize = row_label_size),
                      cluster_rows = FALSE,
                      clustering_method_rows = "ward.D2",
                      cluster_columns = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_dend = FALSE,
                      heatmap_legend_param = list(
                          title_gp = gpar(fontsize = 16,
                                          fontface = "bold"),
                          labels_gp = gpar(fontsize = 14)
                      )
        )
        
        hm
    })
    
    # Output the heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    output$GSE31279_gene_list_heatmap <- renderPlot({
        req(GSE31279_heatmap())
        draw(GSE31279_heatmap(),
             merge_legend = TRUE) # Put Z-score scale and cell type legend in same column
    })
    
    # If the heatmap is successfully made display button to download/save heatmap
    output$GSE31279_heatmap_download_button <- renderUI({
        req(gene_list())
        req(GSE31279_gene_list_df_scaled())
        req(GSE31279_heatmap())
        downloadButton("GSE31279_heatmap_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save heatmap button clicked
    output$GSE31279_heatmap_download <- downloadHandler(
        filename = function() {
            paste("GSE31279_CRC_heatmap", ".png", sep = "")
        },
        content = function(file) {
            png(file, width = 7, height = 7, units = "in", res = 300)
            draw(GSE31279_heatmap(),
                 merge_legend = TRUE)
            dev.off()
        }
    )
    
    
    ##### GSE35602 heatmap #####
    
    # If at least one gene in the user's input list is present
    # then select the rows of the dataframe which correspond to genes
    # present in the user's input gene list
    GSE35602_gene_list_df <- reactive({
        req(GSE35602_gene_list_present() == TRUE)
        found_genes <- gene_list()[gene_list() %in% rownames(GSE35602_exp)]
        GSE35602_exp[found_genes,]
    })
    
    # Z-score scale the expression of the genes present in the
    # user's input gene list
    GSE35602_gene_list_df_scaled <- reactive({
        req(GSE35602_gene_list_present() == TRUE)
        t(scale(t(GSE35602_gene_list_df())))
    })
    
    # Create a heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    GSE35602_heatmap <- reactive({
        req(gene_list())
        validate(
            need(GSE35602_gene_list_present(), "Chosen gene(s) not found in this dataset")
        )
        top_annot <- HeatmapAnnotation(
            "Cell Type" = anno_block(
                gp = gpar(fill = c(
                    "Epithelium" = "#5D69B1",
                    "Stroma" = "#E58606"
                )),
                labels = c("Epithelium", "Stroma"),
                labels_gp = gpar(col = "white",
                                 fontface = "bold")
            )
        )
        
        # Determine the font size to use to label genes (rows) when
        # plotting the heatmap
        row_label_size <- 400 / (nrow(GSE35602_gene_list_df()) * 1.5)
        # Make max possible font size 10
        if (row_label_size > 10){
            row_label_size <- 10
        }
        # Make minimum possible font size 1
        else if (row_label_size < 1){
            row_label_size <- 1
        }
        
        # Create heatmap of Z-score scale gene expression
        hm <- Heatmap(GSE35602_gene_list_df_scaled(),
                      name = "Z-score",
                      top_annotation = top_annot,
                      column_split = GSE35602_cell_types,
                      column_title = NULL,
                      show_column_names = FALSE,
                      row_names_gp = gpar(fontsize = row_label_size),
                      cluster_rows = FALSE,
                      clustering_method_rows = "ward.D2",
                      cluster_columns = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_dend = FALSE,
                      heatmap_legend_param = list(
                          title_gp = gpar(fontsize = 16,
                                          fontface = "bold"),
                          labels_gp = gpar(fontsize = 14)
                      )
        )
        
        hm
    })
    
    # Output the heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    output$GSE35602_gene_list_heatmap <- renderPlot({
        req(GSE35602_heatmap())
        draw(GSE35602_heatmap(),
             merge_legend = TRUE) # Put Z-score scale and cell type legend in same column
    })
    
    # If the heatmap is successfully made display button to download/save heatmap
    output$GSE35602_heatmap_download_button <- renderUI({
        req(gene_list())
        req(GSE35602_gene_list_df_scaled())
        req(GSE35602_heatmap())
        downloadButton("GSE35602_heatmap_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save heatmap button clicked
    output$GSE35602_heatmap_download <- downloadHandler(
        filename = function() {
            paste("GSE35602_CRC_heatmap", ".png", sep = "")
        },
        content = function(file) {
            png(file, width = 7, height = 7, units = "in", res = 300)
            draw(GSE35602_heatmap(),
                 merge_legend = TRUE)
            dev.off()
        }
    )
    
    
    
    ##### GSE164665 heatmap #####
    
    # If at least one gene in the user's input list is present
    # then select the rows of the dataframe which correspond to genes
    # present in the user's input gene list
    GSE164665_gene_list_df <- reactive({
        req(GSE164665_gene_list_present() == TRUE)
        found_genes <- gene_list()[gene_list() %in% rownames(GSE164665_vst_counts)]
        GSE164665_vst_counts[found_genes,]
    })
    
    # Z-score scale the expression of the genes present in the
    # user's input gene list
    GSE164665_gene_list_df_scaled <- reactive({
        req(GSE164665_gene_list_present() == TRUE)
        t(scale(t(GSE164665_gene_list_df())))
    })
    
    # Create a heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    GSE164665_heatmap <- reactive({
        req(gene_list())
        validate(
            need(GSE164665_gene_list_present(), "Chosen gene(s) not found in this dataset")
        )
        top_annot <- HeatmapAnnotation(
            "Cell Type" = anno_block(
                gp = gpar(fill = c(
                    "Epithelium" = "#5D69B1",
                    "Stroma" = "#E58606"
                )),
                labels = c("Epithelium", "Stroma"),
                labels_gp = gpar(col = "white",
                                 fontface = "bold")
            )
        )
        
        # Determine the font size to use to label genes (rows) when
        # plotting the heatmap
        row_label_size <- 400 / (nrow(GSE164665_gene_list_df()) * 1.5)
        # Make max possible font size 10
        if (row_label_size > 10){
            row_label_size <- 10
        }
        # Make minimum possible font size 1
        else if (row_label_size < 1){
            row_label_size <- 1
        }
        
        # Create heatmap of Z-score scale gene expression
        hm <- Heatmap(GSE164665_gene_list_df_scaled(),
                      name = "Z-score",
                      top_annotation = top_annot,
                      column_split = GSE164665_cell_types,
                      column_title = NULL,
                      show_column_names = FALSE,
                      row_names_gp = gpar(fontsize = row_label_size),
                      cluster_rows = FALSE,
                      clustering_method_rows = "ward.D2",
                      cluster_columns = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_dend = FALSE,
                      heatmap_legend_param = list(
                          title_gp = gpar(fontsize = 16,
                                          fontface = "bold"),
                          labels_gp = gpar(fontsize = 14)
                      )
        )
        
        hm
    })
    
    # Output the heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    output$GSE164665_gene_list_heatmap <- renderPlot({
        req(GSE164665_heatmap())
        draw(GSE164665_heatmap(),
             merge_legend = TRUE) # Put Z-score scale and cell type legend in same column
    })
    
    # If the heatmap is successfully made display button to download/save heatmap
    output$GSE164665_heatmap_download_button <- renderUI({
        req(gene_list())
        req(GSE164665_gene_list_df_scaled())
        req(GSE164665_heatmap())
        downloadButton("GSE164665_heatmap_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save heatmap button clicked
    output$GSE164665_heatmap_download <- downloadHandler(
        filename = function() {
            paste("GSE164665_PDAC_heatmap", ".png", sep = "")
        },
        content = function(file) {
            png(file, width = 7, height = 7, units = "in", res = 300)
            draw(GSE164665_heatmap(),
                 merge_legend = TRUE)
            dev.off()
        }
    )
    
    
    #### GSE14548 heatmap ####
    
    # If at least one gene in the user's input list is present
    # then select the rows of the dataframe which correspond to genes
    # present in the user's input gene list
    GSE14548_gene_list_df <- reactive({
        req(gene_list())
        req(GSE14548_gene_list_present() == TRUE)
        found_genes <- gene_list()[gene_list() %in% rownames(GSE14548_exp)]
        GSE14548_exp[found_genes,]
    })
    
    # Z-score scale the expression of the genes present in the
    # user's input gene list
    GSE14548_gene_list_df_scaled <- reactive({
        req(GSE14548_gene_list_present() == TRUE)
        t(scale(t(GSE14548_gene_list_df())))
    })
    
    
    # Create a heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    GSE14548_heatmap <- reactive({
        req(gene_list())
        validate(
            need(GSE14548_gene_list_present(), "Chosen gene(s) not found in this dataset")
        )
        top_annot <- HeatmapAnnotation(
            "Cell Type" = anno_block(
                gp = gpar(fill = c(
                    "Epithelium" = "#5D69B1",
                    "Stroma" = "#E58606"
                )),
                labels = c("Epithelium", "Stroma"),
                labels_gp = gpar(col = "white",
                                 fontface = "bold")
            )
        )
        
        # Determine the font size to use to label genes (rows) when
        # plotting the heatmap
        row_label_size <- 400 / (nrow(GSE14548_gene_list_df()) * 1.5)
        # Make max possible font size 10
        if (row_label_size > 10){
            row_label_size <- 10
        }
        # Make minimum possible font size 1
        else if (row_label_size < 1){
            row_label_size <- 1
        }
        
        # Create heatmap of Z-score scale gene expression
        hm <- Heatmap(GSE14548_gene_list_df_scaled(),
                      name = "Z-score",
                      top_annotation = top_annot,
                      column_split = GSE14548_cell_types,
                      column_title = NULL,
                      show_column_names = FALSE,
                      row_names_gp = gpar(fontsize = row_label_size),
                      cluster_rows = FALSE,
                      clustering_method_rows = "ward.D2",
                      cluster_columns = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_dend = FALSE,
                      heatmap_legend_param = list(
                          title_gp = gpar(fontsize = 16,
                                          fontface = "bold"),
                          labels_gp = gpar(fontsize = 14)
                      )
        )
        
        hm
    })
    
    # Output the heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    output$GSE14548_gene_list_heatmap <- renderPlot({
        req(GSE14548_heatmap())
        draw(GSE14548_heatmap(),
             merge_legend = TRUE) # Put Z-score scale and cell type legend in same column
    })
    
    # If the heatmap is successfully made display button to download/save heatmap
    output$GSE14548_heatmap_download_button <- renderUI({
        req(gene_list())
        req(GSE14548_gene_list_df_scaled())
        req(GSE14548_heatmap())
        downloadButton("GSE14548_heatmap_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save heatmap button clicked
    output$GSE14548_heatmap_download <- downloadHandler(
        filename = function() {
            paste("GSE14548_Breast_Cancer_heatmap", ".png", sep = "")
        },
        content = function(file) {
            png(file, width = 7, height = 7, units = "in", res = 300)
            draw(GSE14548_heatmap(),
                 merge_legend = TRUE)
            dev.off()
        }
    )
    
    
    #### GSE9899 heatmap ####
    
    # If at least one gene in the user's input list is present
    # then select the rows of the dataframe which correspond to genes
    # present in the user's input gene list
    GSE9899_gene_list_df <- reactive({
        req(gene_list())
        req(GSE9899_gene_list_present() == TRUE)
        found_genes <- gene_list()[gene_list() %in% rownames(GSE9899_exp)]
        GSE9899_exp[found_genes,]
    })
    
    # Z-score scale the expression of the genes present in the
    # user's input gene list
    GSE9899_gene_list_df_scaled <- reactive({
        req(GSE9899_gene_list_present() == TRUE)
        t(scale(t(GSE9899_gene_list_df())))
    })
    
    
    # Create a heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    GSE9899_heatmap <- reactive({
        req(gene_list())
        validate(
            need(GSE9899_gene_list_present(), "Chosen gene(s) not found in this dataset")
        )
        top_annot <- HeatmapAnnotation(
            "Cell Type" = anno_block(
                gp = gpar(fill = c(
                    "Epithelium" = "#5D69B1",
                    "Stroma" = "#E58606"
                )),
                labels = c("Epithelium", "Stroma"),
                labels_gp = gpar(col = "white",
                                 fontface = "bold")
            )
        )
        
        # Determine the font size to use to label genes (rows) when
        # plotting the heatmap
        row_label_size <- 400 / (nrow(GSE9899_gene_list_df()) * 1.5)
        # Make max possible font size 10
        if (row_label_size > 10){
            row_label_size <- 10
        }
        # Make minimum possible font size 1
        else if (row_label_size < 1){
            row_label_size <- 1
        }
        
        # Create heatmap of Z-score scale gene expression
        hm <- Heatmap(GSE9899_gene_list_df_scaled(),
                      name = "Z-score",
                      top_annotation = top_annot,
                      column_split = GSE9899_cell_types,
                      column_title = NULL,
                      show_column_names = FALSE,
                      row_names_gp = gpar(fontsize = row_label_size),
                      cluster_rows = FALSE,
                      clustering_method_rows = "ward.D2",
                      cluster_columns = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_dend = FALSE,
                      heatmap_legend_param = list(
                          title_gp = gpar(fontsize = 16,
                                          fontface = "bold"),
                          labels_gp = gpar(fontsize = 14)
                      )
        )
        
        hm
    })
    
    # Output the heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    output$GSE9899_gene_list_heatmap <- renderPlot({
        req(GSE9899_heatmap())
        draw(GSE9899_heatmap(),
             merge_legend = TRUE) # Put Z-score scale and cell type legend in same column
    })
    
    # If the heatmap is successfully made display button to download/save heatmap
    output$GSE9899_heatmap_download_button <- renderUI({
        req(gene_list())
        req(GSE9899_gene_list_df_scaled())
        req(GSE9899_heatmap())
        downloadButton("GSE9899_heatmap_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save heatmap button clicked
    output$GSE9899_heatmap_download <- downloadHandler(
        filename = function() {
            paste("GSE9899_Ovarian_Cancer_heatmap", ".png", sep = "")
        },
        content = function(file) {
            png(file, width = 7, height = 7, units = "in", res = 300)
            draw(GSE9899_heatmap(),
                 merge_legend = TRUE)
            dev.off()
        }
    )
    
    
    #### GSE97284 heatmap ####
    
    # If at least one gene in the user's input list is present
    # then select the rows of the dataframe which correspond to genes
    # present in the user's input gene list
    GSE97284_gene_list_df <- reactive({
        req(gene_list())
        req(GSE97284_gene_list_present() == TRUE)
        found_genes <- gene_list()[gene_list() %in% rownames(GSE97284_exp)]
        GSE97284_exp[found_genes,]
    })
    
    # Z-score scale the expression of the genes present in the
    # user's input gene list
    GSE97284_gene_list_df_scaled <- reactive({
        req(GSE97284_gene_list_present() == TRUE)
        t(scale(t(GSE97284_gene_list_df())))
    })
    
    
    # Create a heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    GSE97284_heatmap <- reactive({
        req(gene_list())
        validate(
            need(GSE97284_gene_list_present(), "Chosen gene(s) not found in this dataset")
        )
        top_annot <- HeatmapAnnotation(
            "Cell Type" = anno_block(
                gp = gpar(fill = c(
                    "Epithelium" = "#5D69B1",
                    "Stroma" = "#E58606"
                )),
                labels = c("Epithelium", "Stroma"),
                labels_gp = gpar(col = "white",
                                 fontface = "bold")
            )
        )
        
        # Determine the font size to use to label genes (rows) when
        # plotting the heatmap
        row_label_size <- 400 / (nrow(GSE97284_gene_list_df()) * 1.5)
        # Make max possible font size 10
        if (row_label_size > 10){
            row_label_size <- 10
        }
        # Make minimum possible font size 1
        else if (row_label_size < 1){
            row_label_size <- 1
        }
        
        # Create heatmap of Z-score scale gene expression
        hm <- Heatmap(GSE97284_gene_list_df_scaled(),
                      name = "Z-score",
                      top_annotation = top_annot,
                      column_split = GSE97284_cell_types,
                      column_title = NULL,
                      show_column_names = FALSE,
                      row_names_gp = gpar(fontsize = row_label_size),
                      cluster_rows = FALSE,
                      clustering_method_rows = "ward.D2",
                      cluster_columns = FALSE,
                      cluster_column_slices = FALSE,
                      show_column_dend = FALSE,
                      heatmap_legend_param = list(
                          title_gp = gpar(fontsize = 16,
                                          fontface = "bold"),
                          labels_gp = gpar(fontsize = 14)
                      )
        )
        
        hm
    })
    
    # Output the heatmap of the z-score scaled expression of the genes
    # present in the user's input gene list
    output$GSE97284_gene_list_heatmap <- renderPlot({
        req(GSE97284_heatmap())
        draw(GSE97284_heatmap(),
             merge_legend = TRUE) # Put Z-score scale and cell type legend in same column
    })
    
    # If the heatmap is successfully made display button to download/save heatmap
    output$GSE97284_heatmap_download_button <- renderUI({
        req(gene_list())
        req(GSE97284_gene_list_df_scaled())
        req(GSE97284_heatmap())
        downloadButton("GSE97284_heatmap_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save heatmap button clicked
    output$GSE97284_heatmap_download <- downloadHandler(
        filename = function() {
            paste("GSE97284_Prostate_Cancer_heatmap", ".png", sep = "")
        },
        content = function(file) {
            png(file, width = 7, height = 7, units = "in", res = 300)
            draw(GSE97284_heatmap(),
                 merge_legend = TRUE)
            dev.off()
        }
    )
    
    #### GSEA TAB ####
    
    observeEvent(input$gsea_gene_list_submit, {
        req(gsea_gene_list_present())
        # Initiate progress bar for GSEA
        show_modal_progress_line(0.1, "Performing GSEA (10%)", easing = "easeOut", duration = 3000, color = "red")
        
        # Update progress bar after GSEA completed for GSE39396
        try(isTruthy(GSE39396_gsea_plot_gg_obj()))
        update_modal_progress(0.2, "Performing GSEA (20%)")
        
        # Update progress bar after GSEA completed for GSE81838
        try(isTruthy(GSE81838_gsea_plot_gg_obj()))
        update_modal_progress(0.35, "Performing GSEA (35%)")
        
        # Update progress bar after GSEA completed for GSE31279
        try(isTruthy(GSE31279_gsea_plot_gg_obj()))
        update_modal_progress(0.5, "Performing GSEA (50%)")
        
        # Update progress bar after GSEA completed for GSE35602
        try(isTruthy(GSE35602_gsea_plot_gg_obj()))
        update_modal_progress(0.6, "Performing GSEA (60%)")
        
        # Update progress bar after GSEA completed for GSE164665
        try(isTruthy(GSE164665_gsea_plot_gg_obj()))
        update_modal_progress(0.75, "Performing GSEA (75%)")
        
        # Update progress bar after GSEA completed for GSE164665
        try(isTruthy(GSE14548_gsea_plot_gg_obj()))
        update_modal_progress(0.85, "Performing GSEA (85%)")
        
        # Update progress bar after GSEA completed for GSE164665
        try(isTruthy(GSE9899_gsea_plot_gg_obj()))
        update_modal_progress(0.95, "Performing GSEA (95%)")
        
        # Update progress bar after GSEA completed for GSE164665
        try(isTruthy(GSE97284_gsea_plot_gg_obj()))
        update_modal_progress(0.1, "Performing GSEA (100%)")

        # Remove progress bar after GSEA completed for all datasets 
        remove_modal_progress()
    })
    
    
    # When user selects an ontology set the corresponding string
    # which appears at the start of that ontology's genesets
    ontology_identifier <- eventReactive(input$ontology,{
        if (input$ontology == "Hallmark"){
            "HALLMARK"
        }
        else if (input$ontology == "BioCarta"){
            "BIOCARTA"
        }
        else if (input$ontology == "KEGG"){
            "KEGG"
        }
        else if (input$ontology == "Reactome"){
            "REACTOME"
        }
        else if (input$ontology == "PID"){
            "PID"
        }
        else if (input$ontology == "WikiPathways"){
            "WP"
        }
    })
        
    
    # When user selects an ontology create a vector with the name of all
    # the gene sets in that ontology
    ontology_gene_sets <- eventReactive(input$ontology,{
        ontology_gene_set_pattern <- paste0("^", ontology_identifier(), "_")
        ontology_inds <- grep(ontology_gene_set_pattern, names(gene_sets_list))
        ontology_gene_sets <- names(gene_sets_list)[ontology_inds]
        ontology_gene_sets
    })
    
    # Update the gene set selectInput to reflect the available gene sets
    # based on the chosen ontology
    observeEvent(input$ontology, {
        ontology_gene_set_pattern <- paste0("^", ontology_identifier(), "_")
        ontology_gene_sets_available <- gsub(ontology_gene_set_pattern, "", ontology_gene_sets())
        ontology_gene_sets_available <- gsub("_", " ", ontology_gene_sets_available)
        updateSelectizeInput(session = session,
                             "gene_set",
                             label = "Select gene set:",
                             choices = ontology_gene_sets_available,
                             options = list(maxOptions = length(ontology_gene_sets_available)),
                             server = TRUE)
    })
    
    
    # When user clicks the gene_list_submit button on Heatmap tab
    # split the user input at each new line character to obtain
    # a vector of gene symbols
    gsea_gene_list <- eventReactive(input$gsea_gene_list_submit, {
        if (input$gene_list_type == "existing"){
            if (input$gene_set != ""){
                full_geneset_name <- gsub(" ", "_", input$gene_set)
                full_geneset_name <- paste0(ontology_identifier(), "_", full_geneset_name)
                chosen_gene_set <- gene_sets_list[[full_geneset_name]]
                chosen_gene_set
            }
            else {
                ""
            }
        }
        else if (input$gene_list_type == "custom"){
            # Split gene list input on new line character
            genes <- strsplit(input$gsea_gene_list, "\n")[[1]]
            # Keep only unique gene symbols
            genes <- unique(genes)
            # Remove any empty strings
            genes <- genes[genes != ""]
            # Remove any horizontal whitespace (space, tab) at start or end
            # i.e. before or after gene symbol
            genes <- trimws(genes, 
                            which = "both", # leading and trailing whitespace
                            whitespace = "[ \t]" # space and tab
            )
            # Convert gene symbols to upper case in case the user entered them
            # partially or completely in lower case 
            genes <- toupper(genes)
            
            # If no input then make genes an empty string
            if (length(genes) == 0){
                #print("gsea genes nothing")
                ""
            }
            # Else return the unique gene symbols (with any empty strings
            # removed and all in upper case
            else {
                #print(genes)
                genes
            }
        }
    })
    
    
    # When user clicks the gene_list_submit button on GSEA tab
    # check if at least one of the genes in the list of genes
    # entered by the user is present in the dataset
    gsea_gene_list_present <- eventReactive(input$gsea_gene_list_submit, {
        if (sum(gsea_gene_list() %in% rownames(GSE39396_exp)) > 0){
            TRUE
        }
        else if (sum(gsea_gene_list() %in% rownames(GSE81838_exp)) > 0){
            TRUE
        }
        else if (sum(gsea_gene_list() %in% rownames(GSE35602_exp)) > 0){
            TRUE
        }
        else if (sum(gsea_gene_list() %in% rownames(GSE31279_exp)) > 0){
            TRUE
        }
        else if (sum(gsea_gene_list() %in% rownames(GSE14548_exp)) > 0){
            TRUE
        }
        else if (sum(gsea_gene_list() %in% rownames(GSE164665_DGEA_result)) > 0){
            TRUE
        }
        else if (sum(gsea_gene_list() %in% rownames(GSE9899_exp)) > 0){
            TRUE
        }
        else if (sum(gsea_gene_list() %in% rownames(GSE97284_exp)) > 0){
            TRUE
        }
        else {
            FALSE
        }
    })
    
    
    # Create logical vectors indicating if genes from GSEA gene list are
    # present/missing in each dataset
    GSE39396_gsea_gene_list_presence <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list() %in% rownames(GSE39396_exp)
    })
    
    GSE81838_gsea_gene_list_presence <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list() %in% rownames(GSE81838_exp)
    })
    
    GSE31279_gsea_gene_list_presence <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list() %in% rownames(GSE31279_exp)
    })
    
    GSE35602_gsea_gene_list_presence <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list() %in% rownames(GSE35602_exp)
    })
    
    GSE164665_gsea_gene_list_presence <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list() %in% rownames(GSE164665_DGEA_result)
    })
    
    GSE14548_gsea_gene_list_presence <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list() %in% rownames(GSE14548_exp)
    })
    
    GSE9899_gsea_gene_list_presence <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list() %in% rownames(GSE9899_exp)
    })
    
    GSE97284_gsea_gene_list_presence <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list() %in% rownames(GSE97284_exp)
    })
    
    
    
    # Create vectors with genes from GSEA gene list missing in each dataset
    GSE39396_gsea_gene_list_missing_genes <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list()[!GSE39396_gsea_gene_list_presence()]
    })
    
    GSE81838_gsea_gene_list_missing_genes <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list()[!GSE81838_gsea_gene_list_presence()]
    })
    
    GSE31279_gsea_gene_list_missing_genes <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list()[!GSE31279_gsea_gene_list_presence()]
    })
    
    GSE35602_gsea_gene_list_missing_genes <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list()[!GSE35602_gsea_gene_list_presence()]
    })
    
    GSE164665_gsea_gene_list_missing_genes <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list()[!GSE164665_gsea_gene_list_presence()]
    })
    
    GSE14548_gsea_gene_list_missing_genes <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list()[!GSE14548_gsea_gene_list_presence()]
    })
    
    GSE9899_gsea_gene_list_missing_genes <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list()[!GSE9899_gsea_gene_list_presence()]
    })
    
    GSE97284_gsea_gene_list_missing_genes <- eventReactive(input$gsea_gene_list_submit, {
        gsea_gene_list()[!GSE97284_gsea_gene_list_presence()]
    })
    
    # Create vector of genes from user's list of genes that are not present
    # in every dataset
    gsea_gene_list_missing_genes_all_datasets <- 
        eventReactive(input$gsea_gene_list_submit, {
            Reduce(intersect,
                   list(GSE39396_gsea_gene_list_missing_genes(),
                        GSE81838_gsea_gene_list_missing_genes(),
                        GSE31279_gsea_gene_list_missing_genes(),
                        GSE35602_gsea_gene_list_missing_genes(),
                        GSE164665_gsea_gene_list_missing_genes(),
                        GSE14548_gsea_gene_list_missing_genes(),
                        GSE9899_gsea_gene_list_missing_genes(),
                        GSE97284_gsea_gene_list_missing_genes()
                   ))
        })
    
    # Logical variable to indicate if there are genes in the user selected list
    # of GSEA genes that are missing in every dataset - controls a
    # conditionalPanel in UI
    output$gsea_genes_missing_all_datasets <- eventReactive(input$gsea_gene_list_submit, {
        if (gsea_gene_list_present() == FALSE){
            return(FALSE)
        }
        else if (is.null(gsea_gene_list_missing_genes_all_datasets())){
            return(FALSE)
        }
        else if (!is.null(gsea_gene_list_missing_genes_all_datasets()) & length(gsea_gene_list_missing_genes_all_datasets()) > 0){
            return(TRUE)
        }
        else {
            return(FALSE)
        }
    })
    
    # output$gsea_genes_missing_all_datasets will be used in the UI code to control
    # the appearance of a message to users using a conditional panel
    # therefore its state must be passed to the UI even though its value 
    # is never displayed by the UI
    outputOptions(output, "gsea_genes_missing_all_datasets", suspendWhenHidden = FALSE)
    
    # List of genes in GSEA gene list to output in modal in UI as missing from
    # all datasets
    output$gsea_genes_missing_all_datasets_modal_content <- renderUI({
        missing_all <- gsea_gene_list_missing_genes_all_datasets()
        HTML(paste(missing_all, collapse = "<br/>"))
    })
    
    
    # Display modal with missing genes (missing in all datasets) if user clicks
    # the actionLink in UI
    observeEvent(input$gsea_genes_missing_all_datasets, {
        showModal(
            modalDialog(
                htmlOutput("gsea_genes_missing_all_datasets_modal_content"),
                footer = modalButton("Close", icon = icon("remove", lib = "glyphicon")),
                easyClose = TRUE,
                size = "s"
            )
        )
    })
    
    
    # Create vector of genes from user's list of genes that are not present in at least one dataset
    gsea_gene_list_missing_genes_any_dataset <- eventReactive(input$gsea_gene_list_submit, {
        unique(c(GSE39396_gsea_gene_list_missing_genes(),
                 GSE81838_gsea_gene_list_missing_genes(),
                 GSE31279_gsea_gene_list_missing_genes(),
                 GSE35602_gsea_gene_list_missing_genes(),
                 GSE164665_gsea_gene_list_missing_genes(),
                 GSE14548_gsea_gene_list_missing_genes(),
                 GSE9899_gsea_gene_list_missing_genes(),
                 GSE97284_gsea_gene_list_missing_genes()
        ))
    })
    
    # Logical variable to indicate if there genes in the user's list of GSEA
    # genes that are missing in at least one dataset - controls a
    # conditionalPanel in UI
    output$gsea_genes_missing_any_dataset <- eventReactive(input$gsea_gene_list_submit, {
        if (gsea_gene_list_present() == FALSE){
            return(FALSE)
        }
        else if (is.null(gsea_gene_list_missing_genes_any_dataset())){
            return(FALSE)
        }
        else if (!is.null(gsea_gene_list_missing_genes_any_dataset()) & !is.null(gsea_gene_list_missing_genes_all_datasets())){
            missing_any_not_all <- setdiff(gsea_gene_list_missing_genes_any_dataset(), gsea_gene_list_missing_genes_all_datasets())
            if (!is.null(missing_any_not_all) & length(missing_any_not_all) > 0){
                return(TRUE)
            }
            else {
                return(FALSE)
            }
        }
        else {
            return(FALSE)
        }
    })
    
    
    # List of genes (from GSEA gene list) to output in modal in UI as missing
    # from at least one dataset
    output$gsea_genes_missing_any_dataset_modal_content <- renderUI({
        missing_any <- setdiff(gsea_gene_list_missing_genes_any_dataset(), gsea_gene_list_missing_genes_all_datasets())
        HTML(paste(missing_any, collapse = "<br/>"))
    })
    
    
    # Display modal with genes missing from at least one dataset (but not all)
    # if user clicks the actionLink in UI
    observeEvent(input$gsea_genes_missing_any_dataset, {
        showModal(
            modalDialog(
                htmlOutput("gsea_genes_missing_any_dataset_modal_content"),
                footer = modalButton("Close", icon = icon("remove", lib = "glyphicon")),
                easyClose = TRUE,
                size = "s"
            )
        )
    })
    
    # output$gsea_genes_missing_any_dataset will be used in the UI code to control
    # the appearance of a message to users using a conditional panel
    # therefore its state must be passed to the UI even though its value 
    # is never displayed by the UI
    outputOptions(output, "gsea_genes_missing_any_dataset", suspendWhenHidden = FALSE)
    
    
    # When user clicks the choose_gene button check if that user has input some
    # text and also check that the text contains at least one gene symbol
    # which is present in at least one of the datasets
    # If these conditions aren't met issue messages to user to alert user to
    # these problems
    observeEvent(input$gsea_gene_list_submit, {
        hideFeedback("gsea_gene_list")
        hideFeedback("gene_set")
        if (input$gene_list_type == "existing"){
            if (input$gene_set == ""){
                feedbackDanger("gene_set",
                               input$gene_set == "",
                               "Please select a gene set")
            }
            else if (gsea_gene_list_present() == FALSE){
                feedbackWarning("gene_set",
                               gsea_gene_list_present() == FALSE,
                               "None of the genes in this gene set are present in any dataset so GSEA can not be performed")
            }
            else {
                hideFeedback("gene_set")
            }
        }
        else if (input$gene_list_type == "custom"){
            if (length(gsea_gene_list()) == 1){
                # If nothing entered into input box or only new line characters
                if (gsea_gene_list() == ""){
                    # This warning will be issued as the user has not entered
                    # any text/gene symbol (only possible text is a sting of new
                    # line characters)
                    feedbackDanger("gsea_gene_list",
                                   gsea_gene_list() == "",
                                   "Please enter a gene symbol(s)")
                }
                # If the user has input text/gene symbols but none of the text/gene
                # symbols match a gene symbol in any of the datasets
                else if (gsea_gene_list_present() == FALSE){
                    # This warning will be issued as the user 
                    feedbackWarning("gsea_gene_list",
                                    input$gsea_gene_list != "" & gsea_gene_list_present() == FALSE,
                                    "None of the chosen genes are present in any dataset. Please ensure you have typed the gene symbol(s) correctly. Alternatively you may try an alias for the gene symbol(s).")
                }
                else {
                    hideFeedback("gsea_gene_list")
                }
            }
            else {
                if (gsea_gene_list_present() == FALSE){
                    feedbackWarning("gsea_gene_list",
                                    input$gsea_gene_list != "" & gsea_gene_list_present() == FALSE,
                                    "None of the chosen genes are present in any dataset. Please ensure you have typed the gene symbol(s) correctly. Alternatively you may try an alias for the gene symbol(s).")
                }
                else {
                    hideFeedback("gsea_gene_list")
                }
            } 
        }
        
    })
    
    
    
    # Title to use for GSEA plots based on gene set selected/ or based on name
    # input by user for user defined gene lists
    geneset_title <- eventReactive(input$gsea_gene_list_submit, { 
        if (input$gene_list_type == "existing") {
            paste(ontology_identifier(), input$gene_set)
        }
        else if (input$gene_list_type == "custom"){
            if (!is.null(input$custom_gene_set_name)){
                input$custom_gene_set_name
            }
            else {
                ""
            }
        }
    })
    
    
    ##### GSE39396 GSEA #####
    
    # Update the cell type options the user can select for group1 based on the
    # option selected for group2 so that the user can't select the same cell
    # type in both groups
    observeEvent(input$GSE39396_gsea_group2, {

        if (!is.null(input$GSE39396_gsea_group2)){

            updated_choices <- unique(GSE39396_cell_types)[!unique(GSE39396_cell_types) %in% input$GSE39396_gsea_group2]

            updateSelectInput(session,
                              "GSE39396_gsea_group1",
                              label = "Compare:",
                              choices = updated_choices,
                              selected = input$GSE39396_gsea_group1
                              )
        }
        else{
            updateSelectInput(session,
                              "GSE39396_gsea_group1",
                              label = "Compare:",
                              choices = unique(GSE39396_cell_types),
                              selected = input$GSE39396_gsea_group1
                              )
        }
    },
    ignoreNULL = FALSE)
    
    # Update the cell type options the user can select for group2 based on the
    # option selected for group1 so that the user can't select the same cell
    # type in both groups
    observeEvent(input$GSE39396_gsea_group1, {
        if (!is.null(input$GSE39396_gsea_group1)){
            updateSelectInput(session,
                              "GSE39396_gsea_group2",
                              label = "to:",
                              choices = unique(GSE39396_cell_types)[!unique(GSE39396_cell_types) %in% input$GSE39396_gsea_group1],
                              selected = input$GSE39396_gsea_group2
            )
        }
        else if (is.null(input$GSE39396_gsea_group1)){
            updateSelectInput(session,
                              "GSE39396_gsea_group2",
                              label = "to:",
                              choices = unique(GSE39396_cell_types),
                              selected = input$GSE39396_gsea_group2
            )
        }
        
    },
    ignoreNULL = FALSE)
    
    gsea_group1_initial <- TRUE
    gsea_group2_initial <- TRUE
    
    values <- reactiveValues(gsea_group1_initial = 0,
                             gsea_group2_initial = 0,
                             GSE39396_group1_cell_types = "Fibroblasts",
                             GSE39396_group2_cell_types = "Epithelial")

    
    observeEvent(input$update_gsea, {
        values$GSE39396_group1_cell_types <- input$GSE39396_gsea_group1
        values$GSE39396_group2_cell_types <- input$GSE39396_gsea_group2
        feedbackWarning("GSE39396_gsea_group1", is.null(values$GSE39396_group1_cell_types), "Please select at least one cell type")
        feedbackWarning("GSE39396_gsea_group2", is.null(values$GSE39396_group2_cell_types), "Please select at least one cell type")
    })
    
    # Identify which genes which genes from the user-defined GSEA gene list 
    # are present in GSE39396
    GSE39396_gsea_genes <- reactive({
        req(gsea_gene_list_present() == TRUE)
        found_genes <- gsea_gene_list()[gsea_gene_list() %in% rownames(GSE39396_exp)]
        found_genes
    })
    
    # Make logical variable indicating whether to show the inputs for
    # selecting cell types for the GSE39396 dataset to the user or not
    # based on wheteher genes in list are present in the dataset and whether
    # the GSEA has been succesffuly calculated/plotted. This controls
    # the conditionalPanel within the GSE39396 GSEA box
    output$GSE39396_show_comparison_groups <- reactive({
        isTruthy(GSE39396_gsea_genes()) & isTruthy(GSE39396_gsea_plot_gg_obj())
    })
    
    # Make output$GSE39396_show_comparison_groups state be passed to UI even
    # though this output value is never displayed in the UI
    outputOptions(output, "GSE39396_show_comparison_groups", suspendWhenHidden = FALSE)
    

    GSE39396_gsea_group1_inds <- reactive({
        req(gsea_gene_list_present())
        req(GSE39396_gsea_genes())
        GSE39396_cell_types %in% values$GSE39396_group1_cell_types
    })

    GSE39396_gsea_group2_inds <- reactive({
        req(gsea_gene_list_present())
        req(GSE39396_gsea_genes())
        GSE39396_cell_types %in% values$GSE39396_group2_cell_types
    })
    
    
    GSE39396_dgea_groups <- reactive({
        req(GSE39396_gsea_group1_inds())
        req(GSE39396_gsea_group2_inds())
        groups <- rep("Other", length(GSE39396_cell_types))
        groups[GSE39396_gsea_group1_inds()] <- "Group1"
        groups[GSE39396_gsea_group2_inds()] <- "Group2"
        groups <- as.factor(groups)
        groups
    })
    
    # Perform differential analysis
    
    GSE39396_top_table <- reactive({
        req(GSE39396_dgea_groups())
        # Set up the experimental design matrix
        design <- model.matrix(~0 + GSE39396_dgea_groups())
        # Change the column names
        # If all 4 cell types are included in the comparison groups 
        # (Group1 and Group2) then there will be only two columns in
        # the design matrix
        if (ncol(design) == 2){
            colnames(design) <- c("Group1", "Group2")
        }
        # Else if one (or more) cell types is not included in either
        # Group1 or Group2 then there will be three columns in 
        # the design matrix
        else {
            colnames(design) <- c("Group1", "Group2", "Other")
        }
        
        # Fit the linear model to the expression set
        fit <- lmFit(GSE39396_exp, design)
        
        # Set up a contrast matrix with Group1 compared to Group2 
        cont_matrix <- makeContrasts(Group1 - Group2, levels = design)
        
        #  Fit to the contrast matrix
        fit2 <- contrasts.fit(fit, cont_matrix)
        fit2 <- eBayes(fit2)
        
        # Calculate differential expression and order by log-odds
        top_table_res <- topTable(fit2, adjust.method = "fdr", sort.by = "B", number = Inf)
        
        # Return result of topTable 
        top_table_res
        
    })
    
    # Order GSE39396_top_table
    GSE39396_top_table_ordered <- reactive({
        GSE39396_top_table()[order(GSE39396_top_table()$t, decreasing = TRUE),]
    })
    
    # Create ranked gene list (vector)
    GSE39396_ranked_gene_list <- reactive({
        # Get the t statistic for each gene
        ranked_gene_stats <- GSE39396_top_table_ordered()$t
        # Make the names of vector the gene symbols
        names(ranked_gene_stats) <- rownames(GSE39396_top_table_ordered())
        
        ranked_gene_stats
    })
    
    # Put the gene list into correct format for running GSEA with 
    #clusterProfiler i.e. a term2gene dataframe
    term2gene <- reactive({
        req(gsea_gene_list_present())
        data.frame(term = "geneset",
                   gene = gsea_gene_list())
    })
    
    # Perform GSEA for GSE39396
    GSE39396_gsea_res <- reactive({
        validate(
            need(GSE39396_gsea_genes(), "None of the genes are present in this dataset")
        )
        set.seed(123)
        gsea_result <- 
            GSEA(
            geneList = GSE39396_ranked_gene_list(),
            exponent = 1,
            minGSSize = 1,
            maxGSSize = Inf,
            eps = 0,
            pvalueCutoff = 1,
            pAdjustMethod = "BH",
            TERM2GENE = term2gene(),
            TERM2NAME = NA,
            verbose = TRUE,
            seed = 123,
            by = "fgsea",
            nPerm = 10000
        )
        
        # Return GSEA result
        gsea_result
    })
    
    # Make GSEA plot for GSE39396
    GSE39396_gsea_plot_gg_obj <- reactive ({
        plot_gsea(GSE39396_gsea_res(),
                  geneSetID = "geneset",
                  title = geneset_title(),
                  left_group = values$GSE39396_group1_cell_types,
                  right_group = values$GSE39396_group2_cell_types,
                  left_group_label_x_coord = 1500,
                  right_group_label_x_coord = 20500)
    })
    
    # Output GSEA plot for GSE39396
    output$GSE39396_gsea_plot <- renderPlot({
        GSE39396_gsea_plot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save gsea_plit
    output$GSE39396_gsea_plot_download_button <- renderUI({
        req(GSE39396_gsea_genes())
        req(GSE39396_gsea_plot_gg_obj())
        downloadButton("GSE39396_gsea_plot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save gsea_plot button clicked
    output$GSE39396_gsea_plot_download <- downloadHandler(
        filename = function() {
            # If no gene set name entered for user defined signature
            # add "user_defined_gene_set_" prefix to save plot file name
            if (is.null(geneset_title()) | geneset_title() == ""){
                paste("user_defined_gene_set_", "GSE39396_CRC_gsea_plot", ".png", sep = "")
            }
            # Else use gene set name for the save plot filename
            else {
                paste(gsub(" ", "_", geneset_title()), "_GSE39396_CRC_gsea_plot", ".png", sep = "")
            }
        },
        content = function(file) {
            ggsave(file, GSE39396_gsea_plot_gg_obj(), device = "png")
        }
    )
    
    
    ##### GSE81838 GSEA #####
    
    # Identify which genes which genes from the user-defined GSEA gene list 
    # are present in GSE81838
    GSE81838_gsea_genes <- reactive({
        req(gsea_gene_list_present() == TRUE)
        found_genes <- gsea_gene_list()[gsea_gene_list() %in% rownames(GSE81838_exp)]
        found_genes
    })
    
    # Perform differential analysis
    
    GSE81838_top_table <- reactive({
        
        req(GSE81838_gsea_genes())
        
        # Set up the experimental design matrix
        design <- model.matrix(~0 + GSE81838_cell_types)
        
        # Change the column names
        if (ncol(design) == 2){
            colnames(design) <- levels(GSE81838_cell_types)
        }
        
        # Fit the linear model to the expression set
        fit <- lmFit(GSE81838_exp, design)
        
        # Set up a contrast matrix with Group1 compared to Group2 
        cont_matrix <- makeContrasts(Stroma - Epithelium, levels = design)
        
        #  Fit to the contrast matrix
        fit2 <- contrasts.fit(fit, cont_matrix)
        fit2 <- eBayes(fit2)
        
        # Calculate differential expression and order by log-odds
        top_table_res <- topTable(fit2, adjust.method = "fdr", sort.by = "B", number = Inf)
        
        #print(head(top_table_res))
        # Return result of topTable 
        top_table_res
        
    })
    
    # Order GSE81838_top_table
    GSE81838_top_table_ordered <- reactive({
        GSE81838_top_table()[order(GSE81838_top_table()$t, decreasing = TRUE),]
    })
    
    # Create ranked gene list (vector)
    GSE81838_ranked_gene_list <- reactive({
        # Get the t statistic for each gene
        ranked_gene_stats <- GSE81838_top_table_ordered()$t
        # Make the names of vector the gene symbols
        names(ranked_gene_stats) <- rownames(GSE81838_top_table_ordered())
        
        ranked_gene_stats
    })
    
    # Perform GSEA for GSE81838 for the selected gene list/signature
    GSE81838_gsea_res <- reactive({
        validate(
            need(GSE81838_gsea_genes(), "None of the genes are present in this dataset")
        )
        set.seed(123)
        gsea_result <- 
            GSEA(
                geneList = GSE81838_ranked_gene_list(),
                exponent = 1,
                minGSSize = 1,
                maxGSSize = Inf,
                eps = 0,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                TERM2GENE = term2gene(),
                TERM2NAME = NA,
                verbose = TRUE,
                seed = 123,
                by = "fgsea",
                nPerm = 10000
            )
        
        # Return GSEA result
        gsea_result
    })
    
    
    # Create GSEA plot for GSE39396
    GSE81838_gsea_plot_gg_obj <- reactive ({
        plot_gsea(GSE81838_gsea_res(),
                  geneSetID = "geneset",
                  title = geneset_title(),
                  left_group = levels(GSE81838_cell_types)[2],
                  right_group = levels(GSE81838_cell_types)[1],
                  left_group_label_x_coord = 1000,
                  right_group_label_x_coord = nrow(GSE81838_exp) - 500)
    })
    
    # Output GSEA plot for GSE39396
    output$GSE81838_gsea_plot <- renderPlot({
        GSE81838_gsea_plot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save gsea_plit
    output$GSE81838_gsea_plot_download_button <- renderUI({
        req(GSE81838_gsea_genes())
        req(GSE81838_gsea_plot_gg_obj())
        downloadButton("GSE81838_gsea_plot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save gsea_plot button clicked
    output$GSE81838_gsea_plot_download <- downloadHandler(
        filename = function() {
            # If no gene set name entered for user defined signature
            # add "user_defined_gene_set_" prefix to save plot file name
            if (is.null(geneset_title()) | geneset_title() == ""){
                paste("user_defined_gene_set_", "GSE81838_TNBC_gsea_plot", ".png", sep = "")
            }
            # Else use gene set name for the save plot filename
            else {
                paste(gsub(" ", "_", geneset_title()), "_GSE81838_TNBC_gsea_plot", ".png", sep = "")
            }
        },
        content = function(file) {
            ggsave(file, GSE81838_gsea_plot_gg_obj(), device = "png")
        }
    )
    
    
    
    
    
    
    ##### GSE31279 GSEA #####
    
    # Identify which genes which genes from the user-defined GSEA gene list 
    # are present in GSE81838
    GSE31279_gsea_genes <- reactive({
        req(gsea_gene_list_present() == TRUE)
        found_genes <- gsea_gene_list()[gsea_gene_list() %in% rownames(GSE31279_exp)]
        found_genes
    })
    
    # Perform differential analysis
    
    GSE31279_top_table <- reactive({
        
        req(GSE31279_gsea_genes())
        
        # Set up the experimental design matrix
        design <- model.matrix(~0 + GSE31279_cell_types)
        
        # Change the column names
        if (ncol(design) == 2){
            colnames(design) <- levels(GSE31279_cell_types)
        }
        
        # Fit the linear model to the expression set
        fit <- lmFit(GSE31279_exp, design)
        
        # Set up a contrast matrix with Group1 compared to Group2 
        cont_matrix <- makeContrasts(Stroma - Epithelium, levels = design)
        
        #  Fit to the contrast matrix
        fit2 <- contrasts.fit(fit, cont_matrix)
        fit2 <- eBayes(fit2)
        
        # Calculate differential expression and order by log-odds
        top_table_res <- topTable(fit2, adjust.method = "fdr", sort.by = "B", number = Inf)
        
        # Return result of topTable 
        top_table_res
        
    })
    
    # Order GSE31279_top_table
    GSE31279_top_table_ordered <- reactive({
        GSE31279_top_table()[order(GSE31279_top_table()$t, decreasing = TRUE),]
    })
    
    # Create ranked gene list (vector)
    GSE31279_ranked_gene_list <- reactive({
        # Get the t statistic for each gene
        ranked_gene_stats <- GSE31279_top_table_ordered()$t
        # Make the names of vector the gene symbols
        names(ranked_gene_stats) <- rownames(GSE31279_top_table_ordered())
        
        ranked_gene_stats
    })
    
    # Perform GSEA for GSE31279 for the selected gene list/signature
    GSE31279_gsea_res <- reactive({
        validate(
            need(GSE31279_gsea_genes(), "None of the genes are present in this dataset")
        )
        set.seed(123)
        gsea_result <- 
            GSEA(
                geneList = GSE31279_ranked_gene_list(),
                exponent = 1,
                minGSSize = 1,
                maxGSSize = Inf,
                eps = 0,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                TERM2GENE = term2gene(),
                TERM2NAME = NA,
                verbose = TRUE,
                seed = 123,
                by = "fgsea",
                nPerm = 10000
            )
        
        # Return GSEA result
        gsea_result
    })
    
    # Create GSEA plot object for GSE31279
    GSE31279_gsea_plot_gg_obj <- reactive ({
        plot_gsea(GSE31279_gsea_res(),
                  geneSetID = "geneset",
                  title = geneset_title(),
                  left_group = levels(GSE31279_cell_types)[2],
                  right_group = levels(GSE31279_cell_types)[1],
                  left_group_label_x_coord = 1000,
                  right_group_label_x_coord = nrow(GSE31279_exp) - 500)
    })
    
    # Output GSEA plot for GSE31279
    output$GSE31279_gsea_plot <- renderPlot({
        GSE31279_gsea_plot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save gsea_plot
    output$GSE31279_gsea_plot_download_button <- renderUI({
        req(GSE31279_gsea_genes())
        req(GSE31279_gsea_plot_gg_obj())
        downloadButton("GSE31279_gsea_plot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save gsea_plot button clicked
    output$GSE31279_gsea_plot_download <- downloadHandler(
        filename = function() {
            # If no gene set name entered for user defined signature
            # add "user_defined_gene_set_" prefix to save plot file name
            if (is.null(geneset_title()) | geneset_title() == ""){
                paste("user_defined_gene_set_", "GSE31279_CRC_gsea_plot", ".png", sep = "")
            }
            # Else use gene set name for the save plot filename
            else {
                paste(gsub(" ", "_", geneset_title()), "_GSE31279_CRC_gsea_plot", ".png", sep = "")
            }
        },
        content = function(file) {
            ggsave(file, GSE31279_gsea_plot_gg_obj(), device = "png")
        }
    )
    
    
    
    
    ##### GSE35602 GSEA #####
    
    # Identify which genes which genes from the user-defined GSEA gene list 
    # are present in GSE35602
    GSE35602_gsea_genes <- reactive({
        req(gsea_gene_list_present() == TRUE)
        found_genes <- gsea_gene_list()[gsea_gene_list() %in% rownames(GSE35602_exp)]
        found_genes
    })
    
    # Perform differential analysis
    
    GSE35602_top_table <- reactive({
        
        req(GSE35602_gsea_genes())
        
        # Set up the experimental design matrix
        design <- model.matrix(~0 + GSE35602_cell_types)
        
        # Change the column names
        if (ncol(design) == 2){
            colnames(design) <- levels(GSE35602_cell_types)
        }
        
        # Fit the linear model to the expression set
        fit <- lmFit(GSE35602_exp, design)
        
        # Set up a contrast matrix with Group1 compared to Group2 
        cont_matrix <- makeContrasts(Stroma - Epithelium, levels = design)
        
        #  Fit to the contrast matrix
        fit2 <- contrasts.fit(fit, cont_matrix)
        fit2 <- eBayes(fit2)
        
        # Calculate differential expression and order by log-odds
        top_table_res <- topTable(fit2, adjust.method = "fdr", sort.by = "B", number = Inf)
        
        # Return result of topTable 
        top_table_res
        
    })
    
    # Order GSE35602_top_table
    GSE35602_top_table_ordered <- reactive({
        GSE35602_top_table()[order(GSE35602_top_table()$t, decreasing = TRUE),]
    })
    
    # Create ranked gene list (vector)
    GSE35602_ranked_gene_list <- reactive({
        # Get the t statistic for each gene
        ranked_gene_stats <- GSE35602_top_table_ordered()$t
        # Make the names of vector the gene symbols
        names(ranked_gene_stats) <- rownames(GSE35602_top_table_ordered())
        
        ranked_gene_stats
    })
    
    # Perform GSEA for GSE35602 for the selected gene list/signature
    GSE35602_gsea_res <- reactive({
        validate(
            need(GSE35602_gsea_genes(), "None of the genes are present in this dataset")
        )
        set.seed(123)
        gsea_result <- 
            GSEA(
                geneList = GSE35602_ranked_gene_list(),
                exponent = 1,
                minGSSize = 1,
                maxGSSize = Inf,
                eps = 0,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                TERM2GENE = term2gene(),
                TERM2NAME = NA,
                verbose = TRUE,
                seed = 123,
                by = "fgsea",
                nPerm = 10000
            )
        
        # Return GSEA result
        gsea_result
    })
    
    
    # Create GSEA plot for GSE35602
    GSE35602_gsea_plot_gg_obj <- reactive ({
        plot_gsea(GSE35602_gsea_res(),
                  geneSetID = "geneset",
                  title = geneset_title(),
                  left_group = levels(GSE35602_cell_types)[2],
                  right_group = levels(GSE35602_cell_types)[1],
                  left_group_label_x_coord = 1000,
                  right_group_label_x_coord = nrow(GSE35602_exp) - 500)
    })
    
    # Output GSEA plot for GSE35602
    output$GSE35602_gsea_plot <- renderPlot({
        GSE35602_gsea_plot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save gsea_plit
    output$GSE35602_gsea_plot_download_button <- renderUI({
        req(GSE35602_gsea_genes())
        req(GSE35602_gsea_plot_gg_obj())
        downloadButton("GSE35602_gsea_plot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save gsea_plot button clicked
    output$GSE35602_gsea_plot_download <- downloadHandler(
        filename = function() {
            # If no gene set name entered for user defined signature
            # add "user_defined_gene_set_" prefix to save plot file name
            if (is.null(geneset_title()) | geneset_title() == ""){
                paste("user_defined_gene_set_", "GSE35602_CRC_gsea_plot", ".png", sep = "")
            }
            # Else use gene set name for the save plot filename
            else {
                paste(gsub(" ", "_", geneset_title()), "_GSE35602_CRC_gsea_plot", ".png", sep = "")
            }
        },
        content = function(file) {
            ggsave(file, GSE35602_gsea_plot_gg_obj(), device = "png")
        }
    )
    
    
    
    
    ##### GSE164665 GSEA #####
    
    # Identify which genes which genes from the user-defined GSEA gene list 
    # are present in GSE81838
    GSE164665_gsea_genes <- reactive({
        req(gsea_gene_list_present() == TRUE)
        found_genes <- gsea_gene_list()[gsea_gene_list() %in% rownames(GSE164665_DGEA_result)]
        found_genes
    })
    
    # Use precomputed differential analysis as top_table
    GSE164665_top_table <- reactive({
        
        req(GSE164665_gsea_genes())
        
        GSE164665_DGEA_result
    })
    
    # Order GSE164665_top_table
    GSE164665_top_table_ordered <- reactive({
        GSE164665_top_table()[order(GSE164665_top_table()$stat, decreasing = TRUE),]
    })
    
    # Create ranked gene list (vector)
    GSE164665_ranked_gene_list <- reactive({
        # Get the t statistic for each gene
        ranked_gene_stats <- GSE164665_top_table_ordered()$stat
        # Make the names of vector the gene symbols
        names(ranked_gene_stats) <- rownames(GSE164665_top_table_ordered())
        
        ranked_gene_stats
    })
    
    # Perform GSEA for GSE164665 for the selected gene list/signature
    GSE164665_gsea_res <- reactive({
        validate(
            need(GSE164665_gsea_genes(), "None of the genes are present in this dataset")
        )
        set.seed(123)
        gsea_result <- 
            GSEA(
                geneList = GSE164665_ranked_gene_list(),
                exponent = 1,
                minGSSize = 1,
                maxGSSize = Inf,
                eps = 0,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                TERM2GENE = term2gene(),
                TERM2NAME = NA,
                verbose = TRUE,
                seed = 123,
                by = "fgsea",
                nPerm = 10000
            )
        
        # Return GSEA result
        gsea_result
    })
    
    
    # Create GSEA plot for GSE164665
    GSE164665_gsea_plot_gg_obj <- reactive ({
        plot_gsea(GSE164665_gsea_res(),
                  geneSetID = "geneset",
                  title = geneset_title(),
                  left_group = levels(GSE164665_cell_types)[2],
                  right_group = levels(GSE164665_cell_types)[1],
                  left_group_label_x_coord = 1000,
                  right_group_label_x_coord = nrow(GSE164665_DGEA_result) - 500)
    })
    
    # Output GSEA plot for GSE164665
    output$GSE164665_gsea_plot <- renderPlot({
        GSE164665_gsea_plot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save gsea_plit
    output$GSE164665_gsea_plot_download_button <- renderUI({
        req(GSE164665_gsea_genes())
        req(GSE164665_gsea_plot_gg_obj())
        downloadButton("GSE164665_gsea_plot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save gsea_plot button clicked
    output$GSE164665_gsea_plot_download <- downloadHandler(
        filename = function() {
            # If no gene set name entered for user defined signature
            # add "user_defined_gene_set_" prefix to save plot file name
            if (is.null(geneset_title()) | geneset_title() == ""){
                paste("user_defined_gene_set_", "GSE164665_PDAC_gsea_plot", ".png", sep = "")
            }
            # Else use gene set name for the save plot filename
            else {
                paste(gsub(" ", "_", geneset_title()), "_GSE164665_PDAC_gsea_plot", ".png", sep = "")
            }
        },
        content = function(file) {
            ggsave(file, GSE164665_gsea_plot_gg_obj(), device = "png")
        }
    )
    
    
    
    ##### GSE14548 GSEA #####
    
    # Identify which genes which genes from the user-defined GSEA gene list 
    # are present in GSE14548
    GSE14548_gsea_genes <- reactive({
        req(gsea_gene_list_present() == TRUE)
        found_genes <- gsea_gene_list()[gsea_gene_list() %in% rownames(GSE14548_exp)]
        found_genes
    })
    
    # Make top table equal to the read in differential analysis result (from global.R)
    GSE14548_top_table <- reactive({
        req(GSE14548_gsea_genes())
        GSE14548_DGEA_result
    })
    
    # Order GSE14548_top_table
    GSE14548_top_table_ordered <- reactive({
        GSE14548_top_table()[order(GSE14548_top_table()$t, decreasing = TRUE),]
    })
    
    # Create ranked gene list (vector)
    GSE14548_ranked_gene_list <- reactive({
        # Get the t statistic for each gene
        ranked_gene_stats <- GSE14548_top_table_ordered()$t
        # Make the names of vector the gene symbols
        names(ranked_gene_stats) <- rownames(GSE14548_top_table_ordered())
        
        ranked_gene_stats
    })
    
    # Perform GSEA for GSE14548 for the selected gene list/signature
    GSE14548_gsea_res <- reactive({
        validate(
            need(GSE14548_gsea_genes(), "None of the genes are present in this dataset")
        )
        set.seed(123)
        gsea_result <- 
            GSEA(
                geneList = GSE14548_ranked_gene_list(),
                exponent = 1,
                minGSSize = 1,
                maxGSSize = Inf,
                eps = 0,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                TERM2GENE = term2gene(),
                TERM2NAME = NA,
                verbose = TRUE,
                seed = 123,
                by = "fgsea",
                nPerm = 10000
            )
        
        # Return GSEA result
        gsea_result
    })
    
    
    # Create GSEA plot for GSE39396
    GSE14548_gsea_plot_gg_obj <- reactive ({
        plot_gsea(GSE14548_gsea_res(),
                  geneSetID = "geneset",
                  title = geneset_title(),
                  left_group = levels(GSE14548_cell_types)[2],
                  right_group = levels(GSE14548_cell_types)[1],
                  left_group_label_x_coord = 1000,
                  right_group_label_x_coord = nrow(GSE14548_exp) - 500)
    })
    
    # Output GSEA plot for GSE39396
    output$GSE14548_gsea_plot <- renderPlot({
        GSE14548_gsea_plot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save gsea_plit
    output$GSE14548_gsea_plot_download_button <- renderUI({
        req(GSE14548_gsea_genes())
        req(GSE14548_gsea_plot_gg_obj())
        downloadButton("GSE14548_gsea_plot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save gsea_plot button clicked
    output$GSE14548_gsea_plot_download <- downloadHandler(
        filename = function() {
            # If no gene set name entered for user defined signature
            # add "user_defined_gene_set_" prefix to save plot file name
            if (is.null(geneset_title()) | geneset_title() == ""){
                paste("user_defined_gene_set_", "GSE14548_Breast_Cancer_gsea_plot", ".png", sep = "")
            }
            # Else use gene set name for the save plot filename
            else {
                paste(gsub(" ", "_", geneset_title()), "_GSE14548_Breast_Cancer_gsea_plot", ".png", sep = "")
            }
        },
        content = function(file) {
            ggsave(file, GSE14548_gsea_plot_gg_obj(), device = "png")
        }
    )
    
    
    
    ##### GSE9899 GSEA #####
    
    # Identify which genes which genes from the user-defined GSEA gene list 
    # are present in GSE9899
    GSE9899_gsea_genes <- reactive({
        req(gsea_gene_list_present() == TRUE)
        found_genes <- gsea_gene_list()[gsea_gene_list() %in% rownames(GSE9899_exp)]
        found_genes
    })
    
    # Make top table equal to the read in differential analysis result (from global.R)
    GSE9899_top_table <- reactive({
        req(GSE9899_gsea_genes())
        GSE9899_DGEA_result
    })
    
    # Order GSE9899_top_table
    GSE9899_top_table_ordered <- reactive({
        GSE9899_top_table()[order(GSE9899_top_table()$t, decreasing = TRUE),]
    })
    
    # Create ranked gene list (vector)
    GSE9899_ranked_gene_list <- reactive({
        # Get the t statistic for each gene
        ranked_gene_stats <- GSE9899_top_table_ordered()$t
        # Make the names of vector the gene symbols
        names(ranked_gene_stats) <- rownames(GSE9899_top_table_ordered())
        
        ranked_gene_stats
    })
    
    # Perform GSEA for GSE9899 for the selected gene list/signature
    GSE9899_gsea_res <- reactive({
        validate(
            need(GSE9899_gsea_genes(), "None of the genes are present in this dataset")
        )
        set.seed(123)
        gsea_result <- 
            GSEA(
                geneList = GSE9899_ranked_gene_list(),
                exponent = 1,
                minGSSize = 1,
                maxGSSize = Inf,
                eps = 0,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                TERM2GENE = term2gene(),
                TERM2NAME = NA,
                verbose = TRUE,
                seed = 123,
                by = "fgsea",
                nPerm = 10000
            )
        
        # Return GSEA result
        gsea_result
    })
    
    
    # Create GSEA plot for GSE39396
    GSE9899_gsea_plot_gg_obj <- reactive ({
        plot_gsea(GSE9899_gsea_res(),
                  geneSetID = "geneset",
                  title = geneset_title(),
                  left_group = levels(GSE9899_cell_types)[2],
                  right_group = levels(GSE9899_cell_types)[1],
                  left_group_label_x_coord = 1000,
                  right_group_label_x_coord = nrow(GSE9899_exp) - 500)
    })
    
    # Output GSEA plot for GSE39396
    output$GSE9899_gsea_plot <- renderPlot({
        GSE9899_gsea_plot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save gsea_plit
    output$GSE9899_gsea_plot_download_button <- renderUI({
        req(GSE9899_gsea_genes())
        req(GSE9899_gsea_plot_gg_obj())
        downloadButton("GSE9899_gsea_plot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save gsea_plot button clicked
    output$GSE9899_gsea_plot_download <- downloadHandler(
        filename = function() {
            # If no gene set name entered for user defined signature
            # add "user_defined_gene_set_" prefix to save plot file name
            if (is.null(geneset_title()) | geneset_title() == ""){
                paste("user_defined_gene_set_", "GSE9899_Ovarian_Cancer_gsea_plot", ".png", sep = "")
            }
            # Else use gene set name for the save plot filename
            else {
                paste(gsub(" ", "_", geneset_title()), "_GSE9899_Ovarian_Cancer_gsea_plot", ".png", sep = "")
            }
        },
        content = function(file) {
            ggsave(file, GSE9899_gsea_plot_gg_obj(), device = "png")
        }
    )
    
    
    
    ##### GSE97284 GSEA #####
    
    # Identify which genes which genes from the user-defined GSEA gene list 
    # are present in GSE97284
    GSE97284_gsea_genes <- reactive({
        req(gsea_gene_list_present() == TRUE)
        found_genes <- gsea_gene_list()[gsea_gene_list() %in% rownames(GSE97284_exp)]
        found_genes
    })
    
    # Make top table equal to the read in differential analysis result (from global.R)
    GSE97284_top_table <- reactive({
        req(GSE97284_gsea_genes())
        GSE97284_DGEA_result
    })
    
    # Order GSE97284_top_table
    GSE97284_top_table_ordered <- reactive({
        GSE97284_top_table()[order(GSE97284_top_table()$t, decreasing = TRUE),]
    })
    
    # Create ranked gene list (vector)
    GSE97284_ranked_gene_list <- reactive({
        # Get the t statistic for each gene
        ranked_gene_stats <- GSE97284_top_table_ordered()$t
        # Make the names of vector the gene symbols
        names(ranked_gene_stats) <- rownames(GSE97284_top_table_ordered())
        
        ranked_gene_stats
    })
    
    # Perform GSEA for GSE97284 for the selected gene list/signature
    GSE97284_gsea_res <- reactive({
        validate(
            need(GSE97284_gsea_genes(), "None of the genes are present in this dataset")
        )
        set.seed(123)
        gsea_result <- 
            GSEA(
                geneList = GSE97284_ranked_gene_list(),
                exponent = 1,
                minGSSize = 1,
                maxGSSize = Inf,
                eps = 0,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                TERM2GENE = term2gene(),
                TERM2NAME = NA,
                verbose = TRUE,
                seed = 123,
                by = "fgsea",
                nPerm = 10000
            )
        
        # Return GSEA result
        gsea_result
        
    })
    
    
    # Create GSEA plot for GSE39396
    GSE97284_gsea_plot_gg_obj <- reactive ({
        plot_gsea(GSE97284_gsea_res(),
                  geneSetID = "geneset",
                  title = geneset_title(),
                  left_group = levels(GSE97284_cell_types)[2],
                  right_group = levels(GSE97284_cell_types)[1],
                  left_group_label_x_coord = 1000,
                  right_group_label_x_coord = nrow(GSE97284_exp) - 500)
    })
    
    # Output GSEA plot for GSE39396
    output$GSE97284_gsea_plot <- renderPlot({
        GSE97284_gsea_plot_gg_obj()
    })
    
    # If the boxplot is successfully made then display button to
    # download/save gsea_plit
    output$GSE97284_gsea_plot_download_button <- renderUI({
        req(GSE97284_gsea_genes())
        req(GSE97284_gsea_plot_gg_obj())
        downloadButton("GSE97284_gsea_plot_download",
                       "Save Plot",
                       #class = "rightAlignButton"
        ) 
    })
    
    # Download handling for when download/save gsea_plot button clicked
    output$GSE97284_gsea_plot_download <- downloadHandler(
        filename = function() {
            # If no gene set name entered for user defined signature
            # add "user_defined_gene_set_" prefix to save plot file name
            if (is.null(geneset_title()) | geneset_title() == ""){
                paste("user_defined_gene_set_", "GSE97284_Prostate_Cancer_gsea_plot", ".png", sep = "")
            }
            # Else use gene set name for the save plot filename
            else {
                paste(gsub(" ", "_", geneset_title()), "_GSE97284_Prostate_Cancer_gsea_plot", ".png", sep = "")
            }
        },
        content = function(file) {
            ggsave(file, GSE97284_gsea_plot_gg_obj(), device = "png")
        }
    )
    
    
}) # close shinyServer
