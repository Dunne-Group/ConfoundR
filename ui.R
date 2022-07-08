#### Load required packages ####
library(shiny)
library(shinydashboard)
library(dashboardthemes)
library(shinyFeedback)
library(shinycssloaders)

#### DASHBOARD ####
dashboardPage(
    
  #### HEADER ####
  dashboardHeader(
    title = "ConfoundR"
  ),
  
  #### SIDEBAR ####
  dashboardSidebar(
    
    # Side bar menu items
    sidebarMenu(id="menu1",
                menuItem("Introduction",
                         icon = icon("info-circle"),
                         tabName = "intro"),
                menuItem("Expression Boxplots",
                         icon = icon("chart-bar"),
                         tabName="boxplots"),
                menuItem("Expression Heatmap",
                         icon = icon("th"),
                         tabName = "expression_heatmap"),
                menuItem("GSEA",
                         icon = icon("chart-line"),
                         tabName="gsea")
                
                
    ), # Close sidebarMenu
    
    # Text to display at bottom of side bar
    helpText("Developed by",
             br(),
             "Ryan Byrne",
             br(),
             a("Molecular Pathology Research Group", href="https://dunne-lab.com/", target="_blank"),
             br(),
             "Queen's University Belfast",
             style = "padding-left:1em; padding-right:1em; position:absolute; bottom:1em; text-align:center;"
    ) # Close helpText
  ), # Close dashboardSidebar
  
    
  #### BODY #####
  dashboardBody(
    
    # CSS
    tags$head(
      tags$style(HTML("
        .modal-dialog{
            overflow-y: initial !important
        }
        .modal-body{
            max-height: 75vh;
            overflow-y: auto;
        }
        ")
      )
    ),
    
    
    # JavaScript to load new tab at the top of page
    # when switch tab buttons are clicked
    tags$script(HTML(" 
        $(document).ready(function () {
        
            $('#menu1').on('change', function (e) {
               window.scrollTo(0,0);
               });
        
            $('#switchTab01').on('click', function (e) {
               window.scrollTo(0,0);
               });
        
            $('#switchTab01').on('click', function (e) {
               window.scrollTo(0,0);
               });
            
            $('#switchTab02').on('click', function (e) {
               window.scrollTo(0,0);
               });
               
            $('#switchTab03').on('click', function (e) {
               window.scrollTo(0,0);
               });
            
            $('#switchTab04').on('click', function (e) {
               window.scrollTo(0,0);
               });
            
            $('#switchTab05').on('click', function (e) {
               window.scrollTo(0,0);
               });
               
            $('#switchTab06').on('click', function (e) {
               window.scrollTo(0,0);
               });
               
               });"
    )),
    
    # Enable shinyFeedback to allow helpful feedback messages to be
    # given to user
    useShinyFeedback(),
    
    # UI for each tab
    tabItems(
      
      #### INTRODUCTION TAB ####
      tabItem(tabName="intro",
              
              h2("Introduction", style="padding-left: 1em"),
              
              img(src='ConfoundR_Schematic.png',
                  alt = "ConfoundR Overview",
                  style = "display: block; margin-left: auto;margin-right: auto; width: 70%; margin-top: 2em; margin-bottom: 2em;"
              ),
              
              #### Aims ####
              h3("Aims", style="padding-left: 1em"),
              # Text explaining purpose/capabilities of ConfoundR
              p("ConfoundR is an interactive web application developed in R with Shiny. 
                The goals of this app are to enable users to:",
                style="padding-left: 2em"),
              tags$ol(
                tags$li("Compare the expression of an individual gene in tumour epithelium with the expression of the same gene in stromal/non-epithelial cells", style = "margin-top: 1em; margin-bottom: 1em"),
                tags$li("Compare the expression of multiple genes in tumour epithelium with the expression of the same genes in stromal/non-epithelial cells", style = "margin-bottom: 1em"),
                tags$li("Compare the expression/enrichment of a gene set/gene signature in tumour epithelium with the expression/enrichment of the same gene set/gene signature in stromal/non-epithelial cells"),
                style="padding-left: 5em"
              ),
              
              br(),
              
              p("By enabling the above comparisons the app enables users to identify individual genes, or gene sets/gene signatures which show differential expression/enrichment across tumour cell populations 
                and therefore may be confounded when examined in bulk transcriptmomic tumour samples as a consequence of different proportions of tumour cell populations."
                , style="padding-left: 2em"),
              
              br(),
              
              p("If you are ready to carry out this analysis, you can skip the below introductions, and examine the expression of individual genes using the ", 
                actionButton('switchTab01', label = "Expression Boxplots", style = "padding: 3px 10px;"), 
                " page, or examine the expression of multiple genes using either a heatmap on the ",
                actionButton('switchTab02', label = "Expression Heatmap", style = "padding: 3px 10px;"), 
                " page or the Gene Set Enrichment Analysis (GSEA) method on the ",
                actionButton('switchTab03', label = "GSEA", style = "padding: 3px 10px;"),
                " page. For more information on the datasets used in the ConfoundR app or for more details on the analysis modules available in the app see below."
                , style="padding-left: 2em"),
              
              br(),
              
              ##### Datasets Intro #####
              h3("Datasets", style="padding-left: 1em"),
              
              # Brief intro about all datasets
              p("Eight datasets, available from Gene Expression Omnibus (GEO), are used in the ConfoundR app:",
                style="padding-left: 2em"),
              tags$ul(
                tags$li("three from colorectal cancer (CRC) (GSE39396, GSE35602 & GSE31279)"),
                tags$li("two from breast cancer (GSE8183 & GSE14548) one of which is from triple negative breast cancer (TNBC) (GSE81838)"),
                tags$li("one from pancreatic ductal adenocarcinoma (PDAC) (GSE164665)"),
                tags$li("one from ovarian cancer (GSE9899)"),
                tags$li("one from prostate cancer (GSE97284)"),
                style="padding-left: 5em"
              ),
              p("Each of these datasets consists of matched tumour epithelium and stromal/non-epithelial samples that have been obtained from bulk tumour specimens using either laser capture microdissection (LCM) or fluorescence-activated cell sorting (FACS). 
                To view more detailed descriptions of each of the datasets and the pre-processing steps applied ",
                actionLink("dataset_info_toggle", "click here"),
                ".",
                style="padding-left: 2em"),
              
              
              ###### Detailed Dataset Info ######
              # Conditional panel with detailed info on datasets and pre-processing steps
              conditionalPanel(
                # show panel if dataset_info_toggle actionLink is clicked an even number of times
                condition = "input.dataset_info_toggle % 2 != 0",
                
                # GSE39396 info
                h4("GSE39396", style="padding-left: 3em"),            
                p("Cohort of tumours from 6 CRC patients. 
                  Each tumour was separated by FACS into epithelial cells (EPCAM+), leukocytes (CD45+), endothelial cells (CD31+) and fibroblasts (FAP+). 
                  This cohort was profiled on the Affymetrix HT HG-U133+ PM Array Plate and background corrected, quantile normalised, summarised and log2 transformed using the robust multi-array average (RMA) method. 
                  The probe IDs were aligned to gene symbols using the microarray annotation file. 
                  For genes represented by multiple probes a single measurement for the gene was selected by choosing the probe with highest mean value across the samples."
                  , style="padding-left: 4em"),
                
                # GSE35602 info
                h4("GSE35602", style="padding-left: 3em"),            
                p("Cohort of 13 matched tumour epithelium and stroma LCM from human CRC tissue (4 normal samples were removed). 
                  This cohort was profiled using Agilent-014850 Whole Human Genome Microarray 4x44K G4112F and lowess normalised. The series matrix was downloaded from GEO and the 4 normal samples removed. 
                  The probe IDs were aligned to gene symbols using the microarray annotation file. 
                  For genes represented by multiple probes a single measurement for the gene was selected  by choosing the probe with highest mean value across the samples."
                  , style="padding-left: 4em"),
                
                # GSE31279 info
                h4("GSE31279", style="padding-left: 3em"),
                p("Cohort of 8 matched tumour epithelium and stroma LCM from human CRC tissue and 2 cases with LCM tumour stroma only: GSM775275 and GSM775279. 
                  Profiled by Illumina humanRef-8 v2.0 expression beadchip and quantile normalised. 
                  The series matrix, which included matched whole tumour and normal for 35 patients, with 10 of those patients having LCM tumour stroma +/- epithelium, was downloaded from GEO. 
                  The 18 LCM samples of interest (8 matched tumour epithelium and stroma samples plus the 2 stroma samples with no matched tumour epithelium) were extracted from the matrix. 
                  Genes with no expression across all samples were removed. The probe IDs were aligned to gene symbols using the microarray annotation file. 
                  For genes represented by multiple probes a single measurement for the gene was selected by choosing the probe with highest mean value across the samples."
                  , style="padding-left: 4em"),
                        
                # GSE81838 info
                h4("GSE81838", style="padding-left: 3em"),
                p("Cohort of 10 matched tumour epithelium and stroma LCM from human triple negative breast cancer profiled using the Affymetrix Human Gene 1.0 ST Array. RMA normalisation and log2 transformation performed. Probe IDs were aligned to gene symbol using the microarray annotation file. 
                  For genes represented by multiple probes a single measurement for the gene was selected by choosing the probe with highest mean value across the samples."   
                  , style="padding-left: 4em"),
                    
                # GSE14548 info
                h4("GSE14548", style="padding-left: 3em"),
                p("Cohort of 14 fresh frozen primary breast cancer biopsies separated into the epithelial and stroma compartments using LCM and profiled using the Affymetrix Human X3P Array. In the epithelial compartment, normal and malignant (ductal carcinoma in situ (DCIS) or invasive ductal carcinoma (IDC)) epithelium were captured. 
                  In the stroma compartment normal stroma away from the malignant lesion, the DCIS-associated stroma and/or IDC-associated stroma whenever possible. The series matrix was downloaded from GEO and only the 9 matched IDC epithelium and IDC-associated stroma samples were selected. Probe IDs were aligned to gene symbol using the microarray annotation file. 
                  For genes represented by multiple probes a single measurement for the gene was selected by choosing the probe with highest mean value across the samples."
                  , style="padding-left: 4em"),
                    
                # GSE164665 info
                h4("GSE164665", style="padding-left: 3em"),
                p("Cohort of 19 matched tumour epithelium and stroma LCM samples from human pancreatic ductal adenocarcinoma (PDAC) sequenced using Illumina NextSeq 500. 
                  The matrices with the raw gene read counts for the tumour epithelium and stroma samples were downloaded from GEO and merged. 
                  Genes with low expression (less than 10 counts across all samples or counts in less than three samples were removed). 
                  Normalised counts were calculated using size factors estimated by the DESeq2 function estimateSizeFactors. 
                  Variance stabilised transformed, normalised counts were generated using the vst function from DESeq2 with the blind argument set to FALSE."
                  , style="padding-left: 4em"),
                        
                # GSE9899 info
                h4("GSE9899", style="padding-left: 3em"),
                p("Cohort of 295 human ovarian cancer samples including five matched tumour epithelium and stroma LCM profiled by Affymetrix Human Genome U133 Plus 2.0 Array. 
                  The series matrix was downloaded from GEO and the five matched tumour epithelium and stroma LCM samples were selected. Probe IDs were aligned to gene symbol using the microarray annotation file. 
                  For genes represented by multiple probes a single measurement for the gene was selected by choosing the probe with highest mean value across the samples."
                  , style="padding-left: 4em"),
                
                # GSE97284 info
                h4("GSE97284", style="padding-left: 3em"),
                p("Cohort of LCM tissue specimens from 12 low grade (Gleason 3+3) and 13 high grade (Gleason 8 and higher) radical prostatectomy and 5 cystoprostatectomy cases profiled using the Affymetrix Human Gene 1.0 ST Array. 
                  For each radical prostatectomy case tumour epithelium, prostatic intraepithelial neoplasia (PIN) epithelium and benign epithelium were captured along with adjacent stroma (tumour-associated stroma, PIN-associated stroma and benign-associated stroma). 
                  For cystoprostatectomy cases benign epithelium and adjacent stroma were captured. The series matrix was downloaded from GEO and the 25 matched tumour epithelium and tumour-associated stroma samples were selected. 
                  Probe IDs were aligned to gene symbol using the microarray annotation file. 
                  For genes represented by multiple probes a single measurement for the gene was selected by choosing the probe with highest mean value across the samples."
                  , style="padding-left: 4em")
                
              ), # close conditionalPanel - detailed dataset info
                  
              br(), # line break
              
              
              ##### Expression Boxplots Description #####
              h3("Expression Boxplots", style="padding-left: 1em"),
              
              p("The ", actionButton('switchTab04', label = "Expression Boxplots", style = "padding: 3px 10px;"), " module allows users to compare the expression of a single gene between epithelial and stromal/non-epithelial cells in each of the datasets. 
                This module enables the user to enter a gene symbol and the analysis module will draw boxplots for the expression of the chosen gene in epithelial and non-epithelial cells allowing a visual comparison of gene expresssion (for each dataset in which the chosen gene is found). 
                In addition, a Mann-Whitney U test is performed comparing the expression of the chosen gene in epithelial and stromal/non-epithelial cells and the p-value is displayed on the plots. 
                Boxplots are made using the", a("ggplot2", href="https://ggplot2.tidyverse.org", target="_blank"), "package and the Mann-Whitney U tests are performed via the ", 
                a("ggpubr", href="https://rpkgs.datanovia.com/ggpubr", target="_blank")," package.",
                style="padding-left: 2em"),
              
              br(), # line break
                  
              
              ##### Expression Heatmap Description #####
              h3("Expression Heatmap", style="padding-left: 1em"),
              
              p("The ", actionButton('switchTab05', label = "Expression Heatmap", style = "padding: 3px 10px;"), " module allows users to visually compare the expression of multiple genes between epithelial and stromal/non-epithelial cells in each of the datasets. 
                This module enables the user to enter a list of gene symbols (each gene symbol should be on a new line) and the analysis module will draw a heatmap of the expression of the chosen genes (which are present in each of the datasets). 
                For heatmaps the expression values for each gene are converted to Z-scores for plotting, enabling visual comparison of the expression of each gene across samples. 
                For the RNA-seq dataset (GSE164665) the transformed, normalised counts produced by the DESeq2 function, vst, are used as the expression values and are converted to Z-scores for plotting of the heatmap. 
                In the heatmaps, the samples are split into epithelial and stromal/non-epithelial samples to aid comparison between these populations. 
                Heatmaps are made using the ", a("ComplexHeatmap", href="https://doi.org/10.1093/bioinformatics/btw313", target="_blank"), "R package.",
                style="padding-left: 2em"),
              
              br(), # line break
                  
              ##### GSEA Description #####
              h3("GSEA", style="padding-left: 1em"),
              # Text explaining GSEA module
              p("The ", actionButton('switchTab06', label = "GSEA", style = "padding: 3px 10px;"), " module enables users to perform gene set enrichment analysis (GSEA). 
                Users can select existing gene sets from the Hallmark, KEGG (Kyoto Encyclopedia of Genes and Genomes), BioCarta, Reactome and PID (Pathway Interactions Database) gene set collections, using the dropdown menus provided, 
                or use their own custom gene set by entering a list of gene symbols (each gene symbol should be on a line).",
                style="padding-left: 2em; margin-bottom: 1.5em"),
                  
              ###### Differential expression analysis methods description ######
              h4("DE analysis methods", style="padding-left: 3em"),
              # Text explaining methods used to perform differential analysis prior to pre-ranked GSEA
              p("To perform GSEA it is necessary to first perform differential expression (DE) analysis. 
                This app implements two different methods to conduct DE analysis, depending on the profiling technology used for the transcriptomic profiling of each dataset. 
                For datasets that were profiled using microarray technology DE analysis was carried out using the",
                a("limma", href="https://doi.org/10.1093/nar/gkv007", target="_blank"),
                " R package while for datasets that were profiled using RNA-Seq technology DE analysis was carried out using the",
                a("DESeq2", href="https://doi.org/10.1186/s13059-014-0550-8", target="_blank"),
                " package.",
                "For the LCM experiments DE analysis is conducted comparing stroma to tumour epithelium. 
                For the FACS sorted dataset (GSE39396), by default the DE analysis is performed comparing fibroblasts to tumour epithelial cells. 
                However, the user can change which cell populations they wish to compare, for this FACS sorted dataset, using the “Compare” and “to” selection menus provided within the GSE39396 panel on the GSEA page. 
                Multiple cell populations can be selected as one group and compared to another individual cell population or another group of multiple cell populations (e.g. Fibroblasts & Leukocytes compared to Epithelial cells)."
                , style="padding-left: 4em; margin-bottom: 2em"),
              
              ###### GSEA analysis methods description ######
              h4("GSEA analysis methods", style="padding-left: 3em"),
              # Text explaining methods used to perform gene set enrichment analysis
              p("To perform pre-ranked GSEA first DE analysis is performed in each dataset as described above comparing gene expression in stromal/non-epithelial cells with expression in epithelial cells. 
                Following DE analysis for each dataset genes are ranked based on the t-statistic (limma) or Wald statistic (DESeq2) calculated by the DE analysis package. 
                This ranked gene list is used, in combination with the gene set selected/uploaded by the user, to perform pre-ranked GESA using the GSEA function from the ",
                a("clusterProfiler", href="https://doi.org/10.1089/omi.2011.0118", target="_blank"),
                " R package using the ", a("fgsea", href="https://doi.org/10.1101/060012", target="_blank"),
                "method (specifically the fgseaSimple method) with 10000 permutations and a random seed of 123.",
                style="padding-left: 4em; margin-bottom: 2em"),
              
              ###### GSEA statistics description ######
              h4("GSEA statistics", style="padding-left: 3em"),
              p("GSEA provides three statistics: Enrichment Score (ES), Normalised Enrichment Score (NES) and p-value (p). The meaning of each of these is explained below.", style="padding-left: 4em"),
              tags$dl(
                # Enrichment score (ES) description
                tags$dt("Enrichment Score (ES)"),
                tags$dd("For the LCM datasets, a positive ES indicates that the gene set is more enriched in the stroma compared to the tumour epithelium. For the FACS dataset, a positive ES indicates that the gene set is more enriched in the cell populations selected in the “Compare” box relative to the cell populations selected in the “to” box (by default this is fibroblasts relative to epithelial cells but this can be changed by the user). On the other hand, a negative ES in the LCM datasets indicates that the gene set is less enriched in the stroma compared to the tumour epithelium (or in other words the gene set is more enriched in the tumour epithelium compared to the stroma). For the FACS dataset, a negative ES indicates that the gene set is less enriched in the cell populations selected in the “Compare” box relative to the cell populations selected in the “to” box (or in other words more enriched in the cell populations selected in the “to” box relative to the cell populations selected in the “Compare” box).",
                        style = "margin-bottom: 1em"),
                # Normalised enrichment score (NES) description
                tags$dt("Normalised Enrichment Score (NES)"),
                tags$dd("To calculate the NES, enrichment scores are calculated for 10000 random gene sets, with the same number of genes as the user selected gene set. The ES obtained by the user selected gene set is then divided by the mean of the enrichment scores for the random gene sets to obtain the NES.",
                        style = "margin-bottom: 1em"),
                # p-value description
                tags$dt("p-value (p)"),
                tags$dd("Using the enrichment scores of the 10000 random gene sets, the probability of obtaining an ES as extreme or more extreme than the ES obtained for the user selected gene set is calculated and this is the p-value (p) for the gene set. The p-value can then be used in conjunction with the NES to determine if the gene set is significantly positively/negatively enriched."),
                style="padding-left: 4em"
              )
      ), # Close intro tab 
            
      #### BOXPLOTS TAB ####
      # This tab enables the user to enter a gene symbol and boxplots
      # will be plotted for each dataset comparing the expression
      # of the user selected gene in the stroma compared to the epithelium
      tabItem(tabName="boxplots",
              h2("Expression Boxplots"),
              
              ##### Boxplot user inputs #####
              
              # Input for user to enter a gene symbol they want 
              #boxplots for 
              textInput("chosen_gene",
                        "Enter a gene symbol",
                        value = "",
                        width = NULL,
                        placeholder = "TP53"),
              # Button for user to click after they enter the gene symbol
              actionButton("choose_gene", "Go"),
              
              # Line breaks for spacing
              br(),
              br(),
              br(),
              
              ##### Boxplots for datasets #####
              
              # First row of boxplots
              fluidRow(
                
                # GSE39396 boxplot
                box(
                  title = "GSE39396 - CRC",
                  solidHeader = TRUE,
                  status = "danger", # changes box colour danger is red
                  collapsible = TRUE, # allow box to be collapsed
                  collapsed = FALSE, # box uncollapsed on initialisation
                  width = 6,
                  uiOutput("GSE39396_boxplot_download_button"), # button to download plot 
                  br(), # line break for space between button and plot
                  
                  # Output the boxplot but wrap in withSpinner
                  # which will display a spinner to show that boxplot
                  # is recalculating
                  withSpinner(
                    plotOutput("GSE39396_boxplot"),
                    type = 1,
                    color = "red"
                  )
                ),
                
                
                # GSE35602 boxplot
                box(title = "GSE35602 - CRC",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE35602_boxplot_download_button"), # button to download plot 
                    br(), # line break for space between button and plot
                    
                    # Output the boxplot but wrap in withSpinner
                    # which will display a spinner to show that boxplot
                    # is recalculating
                    withSpinner(
                      plotOutput("GSE35602_boxplot"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                )
              ),
              
              # Second Row
              fluidRow(
                
                # GSE31279 boxplot
                box(title = "GSE31279 - CRC",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE31279_boxplot_download_button"), # button to download plot 
                    br(), # line break for space between button and plot
                    
                    # Output the boxplot but wrap in withSpinner
                    # which will display a spinner to show that boxplot
                    # is recalculating
                    withSpinner(
                      plotOutput("GSE31279_boxplot"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                ),
                
                # GSE81838 boxplot
                box(title = "GSE81838 - TNBC",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE81838_boxplot_download_button"), # button to download plot
                    br(), # line break for space between button and plot
                    
                    # Output the boxplot but wrap in withSpinner
                    # which will display a spinner to show that boxplot
                    # is recalculating
                    withSpinner(
                      plotOutput("GSE81838_boxplot"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                )
                
              ),
              
              # Third row
              fluidRow(
                
                # GSE14548 boxplot
                box(title = "GSE14548 - Breast Cancer",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE14548_boxplot_download_button"), # button to download plot
                    br(), # line break for space between button and plot
                    
                    # Output the boxplot but wrap in withSpinner
                    # which will display a spinner to show that boxplot
                    # is recalculating
                    withSpinner(
                      plotOutput("GSE14548_boxplot"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                ),
                
                
                # GSE164665 boxplot
                box(title = "GSE164665 - PDAC",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE164665_boxplot_download_button"), # button to download plot
                    br(), # line break for space between button and plot
                    
                    # Output the boxplot but wrap in withSpinner
                    # which will display a spinner to show that boxplot
                    # is recalculating
                    withSpinner(
                      plotOutput("GSE164665_boxplot"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                )
              ),
              
              # Fourth row
              fluidRow(
                
                # GSE9899 boxplot
                box(title = "GSE9899 - Ovarian Cancer",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE9899_boxplot_download_button"), # button to download plot
                    br(), # line break for space between button and plot
                    
                    # Output the boxplot but wrap in withSpinner
                    # which will display a spinner to show that boxplot
                    # is recalculating
                    withSpinner(
                      plotOutput("GSE9899_boxplot"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                ),
                
                # GSE97284 boxplot
                box(title = "GSE97284 - Prostate Cancer",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE97284_boxplot_download_button"), # button to download plot
                    br(), # line break for space between button and plot
                    
                    # Output the boxplot but wrap in withSpinner
                    # which will display a spinner to show that boxplot
                    # is recalculating
                    withSpinner(
                      plotOutput("GSE97284_boxplot"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                )
                
              )
      ),
      
      #### HEATMAP TAB ####
      tabItem(tabName = "expression_heatmap",
              h2("Expression Heatmap"),
              
              ##### Heatmap user inputs #####
              
              # Input for user to enter a list of gene symbols
              # Each gene symbol should be on a new line
              textAreaInput(inputId = "gene_list",
                            label = "Enter a list of gene symbols (as shown)",
                            value = "", # value on initialisation
                            width = NULL,
                            placeholder = "TP53\nCCND2\nCDKN2A\nMDM2", # displays example gene list to user
                            rows = 10 # number of rows of input to display at one time
              ),
              
              # Button for user to click to plot heatmaps after they
              # have entered their list of gene symbols
              actionButton("gene_list_submit", "Make Heatmap"),
              
              
              # Line breaks to add space between button above and content below
              br(),
              br(),
              
              # Show message to user if some genes are missing from ALL datasets
              conditionalPanel(
                condition = "output.heatmap_genes_missing_all_datasets == true",
                p("*Some genes in your list are not present in any of the datasets (",
                  actionLink("heatmap_genes_missing_all_datasets", "View genes"), # Link to display missing genes
                  ")")
              ),
              
              # Show message to user if some genes are missing from AT LEAST ONE dataset
              conditionalPanel(
                condition = "output.heatmap_genes_missing_any_dataset == true",
                p("*Some genes in your list are not present in every dataset (",
                  actionLink("heatmap_genes_missing_any_dataset", "View genes"), # Link to display missing genes
                  ")")
              ),
              
              br(),
              
              
              ##### Heatmap plots for datasets #####
              
              # First row
              fluidRow(
                
                # GSE39396 heatmap
                box(title = "GSE39396 - CRC",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE39396_heatmap_download_button"), # button to download plot
                    
                    # Output the heatmap but wrap in withSpinner
                    # which will display a spinner to show that heatmap
                    # being recalculated/replotted
                    withSpinner(
                      plotOutput("GSE39396_gene_list_heatmap"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                ),
                
                # GSE35602 heatmap
                box(title = "GSE35602 - CRC",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE35602_heatmap_download_button"), # button to download plot
                    
                    # Output the heatmap but wrap in withSpinner
                    # which will display a spinner to show that heatmap
                    # being recalculated/replotted
                    withSpinner(
                      plotOutput("GSE35602_gene_list_heatmap"),  # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                )
              ),
              
              # Second row
              fluidRow(
                
                # GSE31279 heatmap
                box(title = "GSE31279 - CRC",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE31279_heatmap_download_button"), # button to download plot
                    
                    # Output the heatmap but wrap in withSpinner
                    # which will display a spinner to show that heatmap
                    # being recalculated/replotted
                    withSpinner(
                      plotOutput("GSE31279_gene_list_heatmap"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                ),
                
                # GSE81838 heatmap
                box(title = "GSE81838 - TNBC",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE81838_heatmap_download_button"), # button to download plot
                    
                    # Output the heatmap but wrap in withSpinner
                    # which will display a spinner to show that heatmap
                    # being recalculated/replotted
                    withSpinner(
                      plotOutput("GSE81838_gene_list_heatmap"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                )
              ),
              
              # Third row
              fluidRow(
                
                # GSE14548 heatmap
                box(title = "GSE14548 - Breast Cancer",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE14548_heatmap_download_button"), # button to download plot
                    
                    # Output the heatmap but wrap in withSpinner
                    # which will display a spinner to show that heatmap
                    # being recalculated/replotted
                    withSpinner(
                      plotOutput("GSE14548_gene_list_heatmap"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                ),
                
                # GSE164665 heatmap
                box(title = "GSE164665 - PDAC",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE164665_heatmap_download_button"), # button to download plot
                    
                    # Output the heatmap but wrap in withSpinner
                    # which will display a spinner to show that heatmap
                    # being recalculated/replotted
                    withSpinner(
                      plotOutput("GSE164665_gene_list_heatmap"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                )
              ),
              
              # Fourth row
              fluidRow(
                
                # GSE9899 heatmap
                box(title = "GSE9899 - Ovarian Cancer",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE9899_heatmap_download_button"), # button to download plot
                    
                    # Output the heatmap but wrap in withSpinner
                    # which will display a spinner to show that heatmap
                    # being recalculated/replotted
                    withSpinner(
                      plotOutput("GSE9899_gene_list_heatmap"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                ),
                
                # GSE97284 heatmap
                box(title = "GSE97284 - Prostate Cancer",
                    solidHeader = TRUE,
                    status = "danger", # changes box colour danger is red
                    collapsible = TRUE, # allow box to be collapsed
                    collapsed = FALSE, # box uncollapsed on initialisation
                    width = 6,
                    uiOutput("GSE97284_heatmap_download_button"), # button to download plot
                    
                    # Output the heatmap but wrap in withSpinner
                    # which will display a spinner to show that heatmap
                    # being recalculated/replotted
                    withSpinner(
                      plotOutput("GSE97284_gene_list_heatmap"), # output boxplot
                      type = 1, # spinner shape
                      color = "red" # spinner colour
                    )
                )
              )
              
      ),
      
      #### GSEA TAB ####
      tabItem(tabName="gsea",
              
              h2("GSEA"),
              
              ##### GSEA user inputs #####
              
              # Radio buttons to let user select whether they wish to use
              # an existing signature/gene list or their own gene list
              radioButtons(
                inputId = "gene_list_type",
                label = "Type of gene list:",
                choiceNames = c("Existing gene set/signature",
                                "User-defined gene list"),
                choiceValues = c("existing",
                                 "custom"),
              ),
              
              # Panel that displays if user wishes to use an existing
              # signature/gene list
              conditionalPanel(
                condition = "input.gene_list_type == 'existing'", # condition so panel only displays if user wishes to use an existing signature/gene list
                # Input to select which ontology/gene set collection
                selectInput("ontology", 
                            label = "Select gene set collection:",
                            choices = c("Hallmark",
                                        "KEGG",
                                        "BioCarta",
                                        "Reactome",
                                        "PID"),
                            selected = "Hallmark" # Hallmark selected on initialisation
                ),
                # Input to select the gene set form a drop down menu
                # Options in drop down menu change based on the
                # ontology/gene set collection select in selectInput above
                selectizeInput("gene_set", 
                               label = "Select gene set:",
                               choices = NULL)
              ),
              
              # Panel that displays if user wishes to use an existing
              # signature/gene list
              conditionalPanel(
                condition = "input.gene_list_type == 'custom'",  # condition so panel only displays if user wishes to use their own gene list
                
                # Input to enter a list of gene symbols to use for GSEA
                # Each gene symbol should be on a new line
                textAreaInput("gsea_gene_list",
                              "Enter a list of gene symbols (as shown)",
                              value = "", # value on initialisation
                              width = NULL,
                              placeholder = "TP53\nCCND2\nCDKN2A\nMDM2", # Example gene list which user will see
                              rows = 10 # Number of rows to display at once. If more genes (rows) entered a scroll bar appears
                ),
                # Input to enter a name for the gene list entered
                # This is optional and if entered the name is used as
                # title of GSEA plots
                textInput(inputId = "custom_gene_set_name",
                          label = "Enter a name for your gene set (optional):",
                          placeholder = "Enter a name")
              ),
              
              # Button user clicks to perform GSEA after selecting an
              # existing gene list/signature or entering their own gene
              #list
              actionButton("gsea_gene_list_submit", "Perform GSEA"),
              
              # Line breaks for space between button and plots
              br(),
              br(),
              
              # Show message to user if some genes are missing from ALL datasets
              conditionalPanel(
                condition = "output.gsea_genes_missing_all_datasets == true",
                p("*Some genes in the list are not present in any of the datasets (",
                  actionLink("gsea_genes_missing_all_datasets", "View genes"),  # Link to display missing genes
                  ")")
              ),
              # Show message to user if some genes are missing from AT LEAST ONE dataset
              conditionalPanel(
                condition = "output.gsea_genes_missing_any_dataset == true",
                p("*Some genes in the list are not present in every dataset (",
                  actionLink("gsea_genes_missing_any_dataset", "View genes"),  # Link to display missing genes
                  ")")
              ),
              
              br(),
              
              
              ##### GSEA plots for datasets ####
              
              # First row of GSEA plots
              fluidRow(
                box(
                  # GSE39396 GSEA plot
                  title = "GSE39396 - CRC",
                  solidHeader = TRUE,
                  status = "danger", # changes box colour danger is red
                  collapsible = TRUE, # allow box to be collapsed
                  collapsed = FALSE, # box uncollapsed on initialisation
                  width = 6,
                  
                  # Panel that will appear to allow user to select which
                  # cell types to compare for this dataset. Panel
                  # will only appear if genes are present in this dataset
                  # and GSEA calculation is successful
                  conditionalPanel(
                    "output.GSE39396_show_comparison_groups == true",
                    column(
                      width = 6,
                      # Input to select one or more cell types
                      # which will be compared to the cell types
                      # selected by user in below selectInput
                      # (GSE39396_gsea_group2) 
                      selectInput("GSE39396_gsea_group1", 
                                  label = "Compare:",
                                  choices = c("Endothelial",
                                              "Epithelial",
                                              "Fibroblasts",
                                              "Leukocytes"),
                                  selected = "Fibroblasts",
                                  multiple = TRUE),
                    ),
                    column(
                      width = 6,
                      # Input to select one or more cell types
                      # which the cell types selected in above 
                      # selectInput (GSE39396_gsea_group1) will
                      # compared against
                      selectInput("GSE39396_gsea_group2", 
                                  label = "to:",
                                  choices = c("Endothelial",
                                              "Epithelial",
                                              "Fibroblasts",
                                              "Leukocytes"),
                                  selected = "Epithelial",
                                  multiple = TRUE)
                    ),
                    
                    # Button to recalculate/replot the GSEA
                    # after user changes the cell types the choose
                    # to compare
                    actionButton(inputId = "update_gsea",
                                 label = "Update Plot",
                                 icon = icon("sync")),
                    
                    # Line breaks for space between buttons
                    br(),
                    br(),
                    
                    # Button to download GSEA plot
                    uiOutput("GSE39396_gsea_plot_download_button"),
                    
                  ), # close conditionalPanel
                  
                  
                  br(), # Line break for space between button and plot
                  
                  # Output the GSEA plot but wrap in withSpinner
                  # which will display a spinner to show that GSEA
                  # plot is being recalculated/replotted
                  withSpinner(
                    plotOutput("GSE39396_gsea_plot"), # output GSEA plot
                    type = 1,
                    color = "red"
                  )
                ),
                box(
                  # GSE35602 GSEA plot
                  title = "GSE35602 - CRC",
                  solidHeader = TRUE,
                  status = "danger", # changes box colour danger is red
                  collapsible = TRUE, # allow box to be collapsed
                  collapsed = FALSE, # box uncollapsed on initialisation
                  width = 6,
                  uiOutput("GSE35602_gsea_plot_download_button"), # Button to download GSEA plot
                  br(), # Line break for space between button and plot
                  
                  # Output the GSEA plot but wrap in withSpinner
                  # which will display a spinner to show that GSEA
                  # plot is being recalculated/replotted
                  withSpinner(
                    plotOutput("GSE35602_gsea_plot"), # output GSEA plot
                    type = 1, # spinner type
                    color = "red" # spinner colour
                  )
                )
              ),
              
              # Second row of GSEA plots
              fluidRow(
                
                # GSE31279 GSEA plot
                box(
                  title = "GSE31279 - CRC",
                  solidHeader = TRUE,
                  status = "danger", # changes box colour danger is red
                  collapsible = TRUE, # allow box to be collapsed
                  collapsed = FALSE, # box uncollapsed on initialisation
                  width = 6,
                  uiOutput("GSE31279_gsea_plot_download_button"), # Button to download GSEA plot
                  br(), # Line break for space between button and plot
                  
                  # Output the GSEA plot but wrap in withSpinner
                  # which will display a spinner to show that GSEA
                  # plot is being recalculated/replotted
                  withSpinner(
                    plotOutput("GSE31279_gsea_plot"), # output GSEA plot
                    type = 1, # spinner type
                    color = "red" # spinner colour
                  ) 
                ), # close box GSE31279
                
                # GSE81838 GSEA plot
                box(
                  title = "GSE81838 - TNBC",
                  solidHeader = TRUE,
                  status = "danger", # changes box colour danger is red
                  collapsible = TRUE, # allow box to be collapsed
                  collapsed = FALSE, # box uncollapsed on initialisation
                  width = 6,
                  uiOutput("GSE81838_gsea_plot_download_button"), # Button to download GSEA plot
                  br(), # Line break for space between button and plot
                  
                  # Output the GSEA plot but wrap in withSpinner
                  # which will display a spinner to show that GSEA
                  # plot is being recalculated/replotted
                  withSpinner(
                    plotOutput("GSE81838_gsea_plot"), # output GSEA plot
                    type = 1, # spinner type
                    color = "red" # spinner colour
                  )
                ) # close box GSE81838
              ), # close fluidRow - second row GSEA plots
              
              # Third row of GSEA plots
              fluidRow(
                
                # GSE14548 GSEA plot
                box(
                  title = "GSE14548 - Breast Cancer",
                  solidHeader = TRUE,
                  status = "danger", # changes box colour danger is red
                  collapsible = TRUE, # allow box to be collapsed
                  collapsed = FALSE, # box uncollapsed on initialisation
                  width = 6,
                  uiOutput("GSE14548_gsea_plot_download_button"), # Button to download GSEA plot
                  br(), # Line break for space between button and plot
                  
                  # Output the GSEA plot but wrap in withSpinner
                  # which will display a spinner to show that GSEA
                  # plot is being recalculated/replotted
                  withSpinner(
                    plotOutput("GSE14548_gsea_plot"), # output GSEA plot
                    type = 1, # spinner type
                    color = "red" # spinner colour
                  )
                ), # close box GSE14548
                
                # GSE164665 GSEA plot
                box(
                  title = "GSE164665 - PDAC",
                  solidHeader = TRUE,
                  status = "danger", # changes box colour danger is red
                  collapsible = TRUE, # allow box to be collapsed
                  collapsed = FALSE, # box uncollapsed on initialisation
                  width = 6,
                  uiOutput("GSE164665_gsea_plot_download_button"), # Button to download GSEA plot
                  br(), # Line break for space between button and plot
                  
                  # Output the GSEA plot but wrap in withSpinner
                  # which will display a spinner to show that GSEA
                  # plot is being recalculated/replotted
                  withSpinner(
                    plotOutput("GSE164665_gsea_plot"), # output GSEA plot
                    type = 1, # spinner type
                    color = "red" # spinner colour
                  )
                ) # close box GSE164665
              ), # close fluidRow - third row GSEA plots
              
              fluidRow(
                
                # GSE9899 GSEA plot
                box(
                  title = "GSE9899 - Ovarian Cancer",
                  solidHeader = TRUE,
                  status = "danger", # changes box colour danger is red
                  collapsible = TRUE, # allow box to be collapsed
                  collapsed = FALSE, # box uncollapsed on initialisation
                  width = 6,
                  uiOutput("GSE9899_gsea_plot_download_button"), # Button to download GSEA plot
                  br(), # Line break for space between button and plot
                  
                  # Output the GSEA plot but wrap in withSpinner
                  # which will display a spinner to show that GSEA
                  # plot is being recalculated/replotted
                  withSpinner(
                    plotOutput("GSE9899_gsea_plot"), # output GSEA plot
                    type = 1, # spinner type
                    color = "red" # spinner colour
                  )
                ), # close box GSE9899
                
                
                # GSE97284 GSEA plot
                box(
                  title = "GSE97284 - Prostate Cancer",
                  solidHeader = TRUE,
                  status = "danger", # changes box colour danger is red
                  collapsible = TRUE, # allow box to be collapsed
                  collapsed = FALSE, # box uncollapsed on initialisation
                  width = 6,
                  uiOutput("GSE97284_gsea_plot_download_button"), # Button to download GSEA plot
                  br(), # Line break for space between button and plot
                  
                  # Output the GSEA plot but wrap in withSpinner
                  # which will display a spinner to show that GSEA
                  # plot is being recalculated/replotted
                  withSpinner(
                    plotOutput("GSE97284_gsea_plot"), # output GSEA plot
                    type = 1, # spinner type
                    color = "red" # spinner colour
                  )
                ) # close box GSE97284
                
              ) # close fluidRow
      ) # close gsea tabItem
    ) # close tabItems
  ), # close dashboardBody
  
  skin = "red" # colour to use for skin
  
) # close dashboardPage
