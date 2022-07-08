#### READ IN AND PREPARE DATASETS ####

##### GSE39396 #####

###### Read in GSE39396 expression matrix ######
GSE39396_exp <- read.delim("./data/GSE39396/GSE39396_collapsed.txt",
                           row.names = 1,
                           header = TRUE,
                           sep = ",")

# Remove the group (extra gene symbol column) and
# selectedRowID (probe ID) columns
GSE39396_exp$group <- NULL
GSE39396_exp$selectedRowID <- NULL

###### Read in GSE39396 sample metadata file ######
GSE39396_meta <- read.delim("./data/GSE39396/GSE39396_sample_metadata.txt",
                            row.names = 1,
                            header = TRUE,
                            sep = "\t")

# Make a factor vector indicating cell type of each sample in GSE39396
GSE39396_cell_types <- factor(GSE39396_meta$Cell_Type,
                              levels = c("Epithelial",
                                         "Leukocytes",
                                         "Endothelial",
                                         "Fibroblasts"))



##### GSE81838 #####

###### Read in GSE81838 expression matrix ######
GSE81838_exp <- read.delim("./data/GSE81838/GSE81838_Expression_Matrix_RMA_Normalised_Log_Transformed_MaxMean_Collapsed.txt",
                           row.names = 1,
                           header = TRUE,
                           sep = "\t")

###### Read in GSE81838 sample metadata file ######
GSE81838_meta <- read.delim("./data/GSE81838/GSE81838_sample_metadata.txt",
                            row.names = 1,
                            header = TRUE,
                            sep = "\t")

# Make a factor vector indicating cell type of each sample in GSE81838
# Get sample type from metadata
GSE81838_cell_types <- as.character(GSE81838_meta$Tissue_type)
# Change "Tumour" to "Epithelium"
GSE81838_cell_types[GSE81838_cell_types == "Tumour"] <- "Epithelium"
# Make vector into a factor
GSE81838_cell_types <- as.factor(GSE81838_cell_types)



##### GSE31279 #####

###### Read in GSE31279 expression matrix ######
GSE31279_exp <- read.delim("./data/GSE31279/GSE31279_Exp_matrix.txt",
                           row.names = 1,
                           header = TRUE,
                           sep = ",")

###### Read in GSE31279 sample metadata file ######
GSE31279_meta <- read.delim("./data/GSE31279/GSE31279_Annotation_file.txt",
                            row.names = 1,
                            header = TRUE,
                            sep = "\t")

# Make a factor vector indicating cell type of each sample in GSE31279_Annotation_file
# Get sample type from metadata
GSE31279_cell_types <- as.factor(GSE31279_meta$Location)



##### GSE35602 #####

###### Read in GSE35602 expression matrix ######
GSE35602_exp <- read.delim("./data/GSE35602/GSE35602_Exp_matrix.txt",
                           row.names = 1,
                           header = TRUE,
                           sep = ",")

###### Read in GSE35602 sample metadata file ######
GSE35602_meta <- read.delim("./data/GSE35602/GSE35602_Annotation_file.txt",
                            row.names = 1,
                            header = TRUE,
                            sep = "\t")

# Keep only meta data for samples contained in expression matrix
GSE35602_meta <- GSE35602_meta[colnames(GSE35602_exp),]

# Make a factor vector indicating cell type of each sample in GSE35602_Annotation_file
# Get sample type from metadata
GSE35602_cell_types <- as.factor(GSE35602_meta$Location)



#### GSE164665 ####

###### Read in GSE164665 raw counts matrix ######
# GSE164665_raw_counts <- read.delim("./data/GSE164665/GSE164665_all_samples_raw_gene_counts_MARCH_and_SEPT_corrected.txt",
#                                    row.names = 1,
#                                    header = TRUE,
#                                    sep = "\t")

###### Read in GSE164665 normalised counts matrix ######
GSE164665_normalised_counts <- read.delim("./data/GSE164665/GSE164665_normalised_counts_filtered_low_expression_genes.txt",
                                          row.names = 1,
                                          header = TRUE,
                                          sep = "\t")

###### Read in GSE164665 vst transformed counts matrix ######
GSE164665_vst_counts <- read.delim("./data/GSE164665/GSE164665_vst_transformed_counts_filtered_low_expression_genes.txt",
                                   row.names = 1,
                                   header = TRUE,
                                   sep = "\t")

##### Read in precomputed GSE164665 DESeq2 DGEA result #####
GSE164665_DGEA_result <- readRDS("./data/GSE164665/GSE164665_stroma_compared_to_epithelium_DGE_result_DESeq2.rds")

###### Read in GSE164665 sample metadata file ######
GSE164665_meta <- read.delim("./data/GSE164665/GSE164665_sample_metadata.txt",
                             row.names = 1,
                             header = TRUE,
                             sep = "\t")

# Convert location column to factor
GSE164665_meta$Location <- as.factor(GSE164665_meta$Location)

# Make a factor vector indicating cell type of each sample in GSE35602_Annotation_file
# Get sample type from metadata
GSE164665_cell_types <- as.factor(GSE164665_meta$Location)


##### GSE14548 #####

###### Read in GSE14548 expression matrix ######
GSE14548_exp <- read.delim("./data/GSE14548/GSE14548_Exp_matrix.txt",
                          row.names = 1,
                          header = TRUE,
                          sep = "\t")

###### Read in GSE14548 sample metadata file ######
GSE14548_meta <- read.delim("./data/GSE14548/GSE14548_sample_metadata.txt",
                           row.names = 1,
                           header = TRUE,
                           sep = "\t")

# Make a factor vector indicating cell type of each sample in GSE14548
# Get sample type from metadata
GSE14548_cell_types <- as.character(GSE14548_meta$Location)
# Make vector into a factor
GSE14548_cell_types <- as.factor(GSE14548_cell_types)

###### Read in precomputed GSE14548 DESeq2 DGEA result ######
GSE14548_DGEA_result <- 
    readRDS("./data/GSE14548/GSE14548_stroma_compared_to_epithelium_DGE_result_limma.rds")



##### GSE9899 #####

###### Read in GSE9899 expression matrix ######
GSE9899_exp <- read.delim("./data/GSE9899/GSE9899_Exp_matrix.txt",
                           row.names = 1,
                           header = TRUE,
                           sep = "\t")

###### Read in GSE9899 sample metadata file ######
GSE9899_meta <- read.delim("./data/GSE9899/GSE9899_sample_metadata.txt",
                            row.names = 1,
                            header = TRUE,
                            sep = "\t")

# Make a factor vector indicating cell type of each sample in GSE9899
# Get sample type from metadata
GSE9899_cell_types <- as.character(GSE9899_meta$Location)
# Make vector into a factor
GSE9899_cell_types <- as.factor(GSE9899_cell_types)

###### Read in precomputed GSE9899 DESeq2 DGEA result ######
GSE9899_DGEA_result <- 
    readRDS("./data/GSE9899/GSE9899_stroma_compared_to_epithelium_DGE_result_limma.rds")



##### GSE97284 #####

###### Read in GSE97284 expression matrix ######
GSE97284_exp <- read.delim("./data/GSE97284/GSE97284_Exp_matrix.txt",
                          row.names = 1,
                          header = TRUE,
                          sep = "\t")

###### Read in GSE97284 sample metadata file ######
GSE97284_meta <- read.delim("./data/GSE97284/GSE97284_sample_metadata.txt",
                           row.names = 1,
                           header = TRUE,
                           sep = "\t")

# Make a factor vector indicating cell type of each sample in GSE97284
# Get sample type from metadata
GSE97284_cell_types <- as.character(GSE97284_meta$Location)
# Make vector into a factor
GSE97284_cell_types <- as.factor(GSE97284_cell_types)

###### Read in precomputed GSE97284 DESeq2 DGEA result ######
GSE97284_DGEA_result <- 
    readRDS("./data/GSE97284/GSE97284_stroma_compared_to_epithelium_DGE_result_limma.rds")




#### READ IN GENESETS ####

# # Split genesets into to a list with each geneset a separate element of the list
# gene_sets_list <-
#     split(gene_sets$gene_symbol, gene_sets$gs_name)

gene_sets_list <- readRDS("./data/genesets/hallmark_selected_c2_genesets_msigdbr_v7_4_1.rds")


#### DEFINE GSEA PLOT FUNCTION ####

# Create custom function to plot a GSEA result (based on gseaplot2 function
# from the clusterProfiler package)
plot_gsea <- 
    function (x, geneSetID, left_group_label_x_coord, right_group_label_x_coord,
              title = "", color = "green", base_size = 11, 
              rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
              ES_geom = "line", left_group = NULL, right_group = NULL) {
        ES_geom <- match.arg(ES_geom, c("line", "dot"))
        geneList <- position <- NULL
        if (length(geneSetID) == 1) {
            gsdata <- enrichplot:::gsInfo(x, geneSetID)
        }
        else {
            gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
        }
        
        p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
            theme(panel.grid.major = element_line(colour = "grey92"), 
                  panel.grid.minor = element_line(colour = "grey92"), 
                  panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
            scale_x_continuous(expand = c(0, 0))
        if (ES_geom == "line") {
            es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
                                  size = 1)
        }
        else {
            es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
                                   size = 1, data = subset(gsdata, position == 1))
        }
        p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8), 
                                      legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
        p.res <- p.res + ylab("Running Enrichment Score (ES)") +
            theme(axis.text.x = element_blank(), 
                  axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                  plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
                                       unit = "cm"))
        # Add line at y = 0
        p.res <- p.res + geom_hline(yintercept = 0, color = "gray70") 
        
        # Get ES, NES and p-value and format for adding to plot
        selected_gene_set_gsea_result <- as.data.frame(x)
        ES <- format(selected_gene_set_gsea_result$enrichmentScore, digits = 3, nsmall = 2)
        NES <- format(selected_gene_set_gsea_result$NES, digits = 3, nsmall = 2)
        p_val <- format(selected_gene_set_gsea_result$pvalue, digits = 3, nsmall = 2)
        
        # Combine ES, NES and p-val into string which will be plot subtitle
        plot_subtitle <- paste("ES =", ES, "  NES =", NES, "  p =", p_val)
        
        # Add a dashed red line at y = EnrichmentScore
        p.res <- p.res +
            geom_hline(yintercept = selected_gene_set_gsea_result$enrichmentScore,
                       linetype = "dashed",
                       color = "red")
        
        # Add title (gene set name) and subtitle (ES, NES, p-value) to plot
        if (is.null(title) | is.na(title) | title == ""){
            p.res <- p.res +
                labs(title = NULL,
                     subtitle = plot_subtitle) + 
                theme(plot.title = element_blank(),
                      plot.subtitle = element_text(hjust = 0.5,
                                                   color = "black",
                                                   face = "bold"))
        }
        else {
            p.res <- p.res +
                labs(title = title,
                     subtitle = plot_subtitle) + 
                theme(plot.title = element_text(hjust = 0.5,
                                                color = "black",
                                                face = "bold"),
                      plot.subtitle = element_text(hjust = 0.5,
                                                   color = "black",
                                                   face = "bold"))
        }
        
        
        i <- 0
        for (term in unique(gsdata$Description)) {
            idx <- which(gsdata$ymin != 0 & gsdata$Description == 
                             term)
            gsdata[idx, "ymin"] <- i
            gsdata[idx, "ymax"] <- i + 1
            i <- i + 1
        }
        p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
                                                                 ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
            theme_classic(base_size) + theme(legend.position = "none", 
                                             plot.margin = margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(), 
                                             axis.text = element_blank(), axis.line.x = element_blank()) + 
            scale_x_continuous(expand = c(0, 0)) + 
            scale_y_continuous(expand = c(0, 0))
        
        if (length(geneSetID) == 1) {
            v <- seq(1, sum(gsdata$position), length.out = 9)
            v <- seq(1, length(gsdata$position), length.out = 11)
            v <- v[-1]
            v <- v[-length(v)]
            
            cross_zero_ind <- sum(gsdata$geneList < 0) + 1
            rev_cross_zero_ind <- length(gsdata$geneList) - sum(gsdata$geneList > 0)
            
            v_pos <- seq(1, cross_zero_ind, length.out = 6)
            v_pos <- v_pos[-1]
            v_neg <- seq(cross_zero_ind, length(gsdata$geneList), length.out = 6)
            v_neg <- v_neg[-1]
            v_neg <- v_neg[-length(v_neg)]
            v <- c(v_pos, v_neg)
            
            inv <- findInterval(rev(cumsum(gsdata$position)), v)
            inv <- findInterval(rev(gsdata$x), v)
            
            if (min(inv) == 0) 
                inv <- inv + 1
            col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
            ymin <- min(p2$data$ymin)
            yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
            xmin <- which(!duplicated(inv))
            xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
            d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                            xmax = xmax, col = col[unique(inv)])
            p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
                                      ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
                                 alpha = 0.9, inherit.aes = FALSE)
        }
        df2 <- p$data
        df2$y <- p$data$geneList[df2$x]
        p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
                                                   y = ~y, yend = 0), color = "grey")
        p.pos <- p.pos + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") + 
            theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, 
                                       l = 0.2, unit = "cm"))
        p.pos <- p.pos + annotate(geom = "text", x=left_group_label_x_coord, y=max(df2$y) - (diff(range(df2$y))/20), label=paste(left_group, collapse = "/"), hjust = "inward", color="black")
        p.pos <- p.pos + annotate(geom = "text", x=right_group_label_x_coord, y=max(df2$y) - (diff(range(df2$y))/20), label=paste(right_group, collapse = "/"), hjust = "inward", color="black")
        
        if (!is.null(title) && !is.na(title) && title != "") 
            p.res <- p.res + ggtitle(title)
        if (length(color) == length(geneSetID)) {
            p.res <- p.res + scale_color_manual(values = color)
            if (length(color) == 1) {
                p.res <- p.res + theme(legend.position = "none")
                p2 <- p2 + scale_color_manual(values = "black")
            }
            else {
                p2 <- p2 + scale_color_manual(values = color)
            }
        }
        if (pvalue_table) {
            pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
            rownames(pd) <- pd$Description
            pd <- pd[, -1]
            pd <- round(pd, 4)
            tp <- tableGrob2(pd, p.res)
            p.res <- p.res + 
                theme(legend.position = "none") + 
                annotation_custom(tp, 
                                  xmin = quantile(p.res$data$x, 0.5), 
                                  xmax = quantile(p.res$data$x, 0.95),
                                  ymin = quantile(p.res$data$runningScore, 0.75),
                                  ymax = quantile(p.res$data$runningScore, 0.9))
        }
        plotlist <- list(p.res, p2, p.pos)[subplots]
        n <- length(plotlist)
        plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                               axis.ticks.x = element_line(), axis.text.x = element_text())
        if (length(subplots) == 1){ 
            return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
                                                              r = 0.2,
                                                              b = 0.2,
                                                              l = 0.2,
                                                              unit = "cm")))
        }
        if (length(rel_heights) > length(subplots)) {
            rel_heights <- rel_heights[subplots]
        }
        plot_grid(plotlist = plotlist,
                  ncol = 1,
                  align = "v",
                  rel_heights = rel_heights)
    }


################################ END ###########################################
