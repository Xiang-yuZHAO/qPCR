#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#https://github.com/daattali/shinyalert
#https://shinyapps.dreamrs.fr/shinyWidgets/
#https://bs4dash.rinterface.com/reference/
library(shiny)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(stringr)
library(tidyverse)
library(bs4Dash)
library(shinyWidgets)
library(readxl)
library(plotly)
library(DT)
library(digest)
library(shinyalert)
library(markdown)
library(knitr)

#### 1. 函数 ####
options(encoding = "UTF-8", warn = -1, 
        show.error.messages = F)


# df <- readxl::read_excel(path = "d:/ct1.xlsx", sheet = 2)
# df <- melt(df, id.vars = "Sample", variable.name = "Gene", value.name = "Ct")
# 
# df_standard <- read.csv(file = "D:/ct3.csv")


# * 1.1 PCR analyze 差异分析 ----

# HKG: housekeeping gene
# GOI: gene of interest
# Control: 对照组
# Treated: 实验组
# df: 数据框，第一列是“Sample”，第二列是“Gene”，第三列是“Ct"

#比较不同样本间(Control和Treated)同一个基因(gene)的倍数差异，不需要设置HKG
pcr_dct <- function(df, Control, Treated, gene) {
  dt <- df %>%
    filter(Sample %in% c(Control, Treated),
           Gene == gene) %>%
    mutate(
      Ct = as.numeric(Ct) ,
      Sample_type = if_else(Sample == Control, "Control", "Treated"),
      Control_mean = mean(Ct[Sample == Control]),
      `2deltaCt` = 2 ^ (Ct - Control_mean)
    ) %>%
    select(-Ct,-Control_mean,-Sample_type) %>%
    rename("Relative expression" = `2deltaCt`) %>%
    as.data.frame()
  return(dt)
  
}

#pcr_dct(df,Control = "brain", Treated = "kidney", gene = "c_myc")

#比较不同样本间(Control和Treated)同一个基因(GOI)的倍数差异，需要HKG
pcr_ddct <- function(df, Control, Treated, HKG, GOI) {
  dt <- df %>%
    filter(Sample %in% c(Control, Treated),
           Gene %in% c(HKG, GOI)) %>%
    mutate(
      Ct = as.numeric(Ct) ,
      Sample_type = if_else(Sample == Control, "Control", "Treated"),
      Gene_class = if_else(Gene == HKG, "HKG", "GOI"),
      HKG_mean = ifelse(Sample == Control, 
                        mean(Ct[Sample == Control & Gene == HKG]),
                        mean(Ct[Sample == Treated & Gene == HKG]))
    ) %>%
    filter(Gene_class == "GOI") %>%
    mutate(
      deltaCt = Ct - HKG_mean
    ) %>%
    as.data.frame()
  
  Calibrator <- mean(dt[dt$Sample_type == "Control",]$deltaCt) # average ∆Ct control group
  
  dt <- dt %>%
    mutate(`2^-deltadeltaCt` = 2^(-(deltaCt - Calibrator))) %>%
    select(Sample,Gene,`2^-deltadeltaCt`) %>%
    rename("Relative expression" = `2^-deltadeltaCt`) %>%
    as.data.frame()
  
  return(dt)
}

#pcr_ddct(df, Control = "lung", Treated = "brain", HKG = "GAPDH", GOI = "c_myc")


pcr_curve <- function(df, df_standard, Control, Treated, HKG, GOI) {
  res_E <- pcr_standard(df = df_standard, plot = F)
  E_GOI <- res_E[row.names(res_E) == GOI,]$E 
  E_HKG <- res_E[row.names(res_E) == HKG,]$E 
  
  dt <- df %>%
    filter(Sample %in% c(Control, Treated),
           Gene %in% c(HKG, GOI)) %>%
    mutate(
      Ct = as.numeric(Ct) ,
      Sample_type = if_else(Sample == Control, "Control", "Treated"),
      Gene_class = if_else(Gene == HKG, "HKG", "GOI"),
      mean_Ct = ifelse(Gene == HKG,
                       mean(Ct[Gene == HKG & Sample == Control]),
                       mean(Ct[Gene == GOI & Sample == Control]))
      
    ) %>%
    mutate(
      delta_Ct = mean_Ct - Ct
    ) %>%
    mutate(
      mean_delta_Ct_HKG = ifelse(Sample == Control,
                                 mean(delta_Ct[Gene == HKG & Sample == Control]),
                                 mean(delta_Ct[Gene == HKG & Sample == Treated])),
    )  %>%
    filter(Gene == GOI) %>%
    mutate(
      geneExpRatio = (E_GOI^delta_Ct)/(E_HKG^mean_delta_Ct_HKG)
    ) %>%
    select(Sample, Gene,geneExpRatio) %>%
    rename("Relative expression" = geneExpRatio)
  
  return(dt)
}

#pcr_curve(df, df_standard, Control = "brain",Treated = "kidney", HKG = "GAPDH", GOI = "c_myc")


pcr_analyze <- function(df, method = 'delta_delta_ct', Control, Treated_samples, ...) {
  pcr_rel <- switch(
    method,
    'delta_delta_ct' = pcr_ddct,
    'delta_ct' = pcr_dct,
    'relative_curve' = pcr_curve
  )
  
  
  results <- lapply(Treated_samples, function(treat_sample) {
    pcr_rel(df, Control = Control, Treated = treat_sample, ...)
  })
  
  results <- do.call(rbind, results)
  
  sample_num <- length(c(Control, Treated_samples))
  rep_num <- (nrow(results)/(sample_num-1))/2
  
  results_control <- results[results$Sample == Control,]
  results_control <- results_control[1:rep_num,]
  
  results <- results[results$Sample != Control,]
  results <- rbind(results_control,results)
  
  # results$Sample <- factor(x = results$Sample, 
  #                          levels = c(Control, Treated_samples))
  
  return(results)

}





# df_plot <- pcr_analyze(df, method = 'delta_delta_ct',
#             Treated_samples = c("kidney","brain"),
#             Control = "lung",
#             HKG = "GAPDH", GOI = "c_myc")
# 
# pcr_analyze(df, method = 'delta_ct', Treated_samples = c("kidney","lung"),
#             Control = "brain", gene = "GAPDH")
# 
# 
# pcr_analyze(df = df, method = 'relative_curve', df_standard = df_standard,
#             Control = "brain",Treated_samples = c("kidney","lung"),
#             HKG = "GAPDH", GOI = "c_myc")


# * 1.2 PCR assess 质量评估 ----
# 输入的数据格式: amount(ng), gene1, gene2, gene3, gene4...
pcr_standard <- function(df, plot = FALSE) {
  if (plot) {
    pcr_plot_assess(df)
  } else {
    # calculate trend; intercep, slop and r_squared
    genes <- colnames(df)[-which(colnames(df) == "amount")]
    lm_results <- lapply(genes, function(gene) {
      lm_model <- lm(formula = df[[gene]] ~ amount, data = df)
      lm_res <- data.frame(intercept = coef(lm_model)[1],
                           slope = coef(lm_model)[2],
                           r_squared = summary(lm_model)$r.squared)
      row.names(lm_res) <- gene
      return(lm_res)
    })
    res <- do.call(rbind, lm_results)
    
    res$E <- 10^(-1/res$slope)
    res$`% Efficiency` <- (res$E -1)*100
    return(res)
  }
}

# pcr_standard(df = df_standard, plot = T)
# pcr_standard(df = df_standard, plot = F)

# * 1.3 PCR test 统计学检验 ----

pcr_test <- function(df_plot, method = "t.test") {
  fun_test <- switch(
    method,
    't.test' = t.test,
    'wilcox.test' = wilcox.test
  )
  
  groups <- unique(as.character(df_plot$Sample))
  group_combn <- combn(groups, 2)
  
  p_values <- sapply(1:ncol(group_combn), function(i) {
    group1 <- group_combn[1, i]
    group2 <- group_combn[2, i]
    
    test_result <- fun_test(df_plot$`Relative expression`[df_plot$Sample == group1], 
                            df_plot$`Relative expression`[df_plot$Sample == group2])
    
    return(test_result$p.value)
  })
  
  res <- data.frame(sample = apply(group_combn, 2, paste, collapse = "-"),
                    p_value = p_values)
  res <- res %>%
    mutate(p_sig = case_when(
      p_value > 0.05 ~ "ns",
      p_value <= 0.05 & p_value > 0.01 ~ "*",
      p_value <= 0.01 & p_value > 0.001 ~ "**",
      p_value <= 0.001 ~ "***"
    ))
  
  
  return(res)
}

#pcr_test(df_plot,method = "t.test")

# * 1.4 其它函数 ----
getName <- function(df, object){
  as.character(unique(df[[object]])) 
}

get_P_pos <- function(df){
  max(as.numeric(df[["Relative expression"]]))+0.5
}

pcr_plot_analyze <- function(df_plot, ref_group, label_pVal_y,
                             test_method = "t.test", p_label = "p.signif") {
  
  other_group <- unique(as.character(df_plot[["Sample"]]))
  other_group <- other_group[!grepl(ref_group, other_group)]
  df_plot$Sample <- factor(x = df_plot$Sample,
                           levels = c(ref_group, other_group))
  
  ggplot(df_plot, aes(x = Sample, y = `Relative expression`, fill = Gene))+
    geom_bar(stat = "summary", fun = "mean",
             position = position_dodge(), alpha = .8) +
    geom_point(position = position_dodge(width = 0.9),size = 1, fill = "black") +
    stat_summary(fun.data = "mean_sd", geom = "errorbar", color = "black",
                 position = position_dodge(.9), width = 0.15)+
    scale_y_continuous(expand = c(0,0), 
                       limits = c(0,as.numeric(label_pVal_y)+1))+ 
    stat_compare_means(
      method = test_method, 
      label = p_label, 
      ref.group = ref_group,
      label.y = label_pVal_y,
      method.args = list(paired = F, p.adjust.method = "none"),
      size = 6, angle = 0,
      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), 
                         symbols = c("***", "**", "*", "ns")))+
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(face = "bold", size = 15, angle = 60, hjust = 1, vjust = 1),
          axis.text.y = element_text(face = "bold", size = 15),
          title = element_text(size = 17),
          legend.position = "none",
          plot.margin = margin(t = 50, l = 10)) 
}
# pcr_plot_analyze(df_plot, test_method = "t.test", 
#                  ref_group = "brain",p_label = "p.format")

# 
# df_standard <- readxl::read_excel(path = "C:/Users/ZHAO/Desktop/qPCR demo.xlsx",
#                                   sheet = 3)

pcr_plot_assess <- function(df) {
  df <- melt(df, id.vars = "amount", variable.name = "gene", value.name = "Ct")
  ggplot(df, aes(x = log10(amount), y =Ct, color = gene)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y~x, se = T) +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(face = "bold", size = 15),
          title = element_text(size = 17))
}
# pcr_plot_assess(df = df_standard)

# df_melt <- readxl::read_excel(path = "C:/Users/ZHAO/Desktop/qPCR demo.xlsx",
#                               sheet = 2)
# df_sample <- readxl::read_excel(path = "C:/Users/ZHAO/Desktop/qPCR demo.xlsx",
#                                 sheet = 4)

pcr_melt <- function(df_melt, df_sample, sample, gene){
  df_plot <- merge(x = df_melt, y = df_sample, by = "Well")
  df_plot <- df_plot[df_plot$Sample == sample & df_plot$Gene == gene,]
  
  ggplot(df_plot, aes(x = Temperature, y = Derivative, color = Well))+
    geom_point(alpha = .4)+
    geom_line() + 
    theme_bw()  +
    theme(axis.title = element_text(face = "bold", size = 15),
          axis.text = element_text(face = "bold", size = 15),
          title = element_text(size = 17),)
  
}
# pcr_melt(df_melt = df_melt,df_sample = df_sample,sample = "lung", gene = "c_myc")

format_scientific <- function(x) {
  format(x, scientific = F, digits = 3)
}




#### 2. UI ####
ui <- dashboardPage(
  dark = NULL,
  help =  NULL,

  # * 2.1 header ----
  header = dashboardHeader(title = "qPCR",
                           dashboardBrand(title = "qPCR shiny app")
                           # tags$head(tags$style(HTML(
                           #   "#myDiv { 
                           #     position: absolute;
                           #     top: 1.4rem!important;
                           #     right: -14rem!important;
                           #     
                           #   }"
                           # ))),
                           # div(
                           #   id = "myDiv",
                           #   switchInput(
                           #     inputId = "Id013",
                           #     onLabel = "中",
                           #     offLabel = "EN",
                           #     size = "mini"
                           #   )
                           #  
                           # )
                          
                           
  ),
  
  # * 2.2 sidebar ----
  sidebar = dashboardSidebar(id = "sidebar",
                             skin = "light",
                             sidebarUserPanel(name = "qPCR Master",
                                              image = "icon/qPCR20231029.png"),  
                             #sidebarHeader(title = "Sidebar"),
                             sidebarMenu(id = "updateTabItems_sidebarmenu",
                                         #sidebarHeader(title = "sidebarmenu header"),
                                         menuItem(text = "Homepage", 
                                                  tabName = "Home",
                                                  icon = shiny::icon("font-awesome")),
                                         menuItem(text = "Data Analysis", 
                                                  tabName = "qpcr_data_process",
                                                  icon = shiny::icon("bar-chart")),
                                         menuItem(text = "qPCR Principles", 
                                                  tabName = "qpcr_basic",
                                                  icon = shiny::icon("file")),
                                         menuItem(text = "Experimental Design", 
                                                  tabName = "qpcr_exp_design",
                                                  icon = shiny::icon("vials")),
                                         menuItem(text = "Analysis Instructions", 
                                                  tabName = "data_process_info",
                                                  icon = shiny::icon("file-text")),
                                         menuItem(text = "User Assistance", tabName = "help",
                                                  icon = shiny::icon("circle-question"))
                             )
  ),
  # controlbar = dashboardControlbar(id = "controlbal",skin = "light",pinned = T,
  #                                  controlbarMenu(
  #                                    controlbarItem(title = "controlbar-item")
  #                                  )),
  # * 2.3 body ----
  body = dashboardBody(
    tags$head(includeHTML(path = "www/google-analytics.html")),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
    ),
    tabItems(
    # ** 2.3.1 主页 ----
     tabItem(tabName = "Home",
             h3("Welcome to the qPCR Relative Expression Data Analysis Tool."),
             #br(),
             p("This platform is designed for researchers and learners to provide qPCR data analysis tools. Whether you are an experienced researcher or a beginner in experimental design and data analysis, you can find the tools and guidance you need on this platform to understand and process your data more intuitively and conveniently."),
             hr(),
             h4("Our website features:"),
             br(),
             div(
               id = "center-image",
               imageOutput(outputId = "img_home",
                           width = "auto", height = "auto")
             ),
             
             p("1. qPCR Experiment Principles: Provides you with fundamental knowledge of qPCR experiments to help you understand the experimental principles more deeply."),
             p("2. Experimental Design: Offers specific steps and important tips for experiment design, making it easy to create feasible experiment plans."),
             p("3. Data Analysis: Our data analysis tools are capable of efficient data processing, supporting data upload, automated analysis, and result visualization."),
             p("4. Analysis Instructions: Provides detailed data analysis explanations and examples to help you better understand your analysis results."),
             p("5. User Assistance: If you encounter any difficulties during usage, our help options can provide answers to your questions and usage tips."),
             hr(),
             h4("To get started with our application:"),
             br(),
             p("1. Upload Your Data: Go to the \"Data Analysis\" option and upload your pre-processed qPCR data."),
             p("2. Choose Analysis Way: Select the analysis way that suit your needs."),
             p("3. View Analysis Results: Review your analysis results and conclusions. If you have any confusion, you can refer to our \"Analysis Instructions\" for understanding."),
             p("4. Seek Help: If you encounter any issues during usage, you can always refer to our \"User Assistance.\""),
             actionBttn(
               inputId = "Id112",
               label = "Start data analysis", 
               style = "stretch",
               color = "primary"
             ),
             hr(),
             strong("News"),
             br(),
             HTML("▪ Release Date: October 28, 2023<br>
                  What's New in This Version:<br>
                  Enhanced Data Visualization; Advanced Statistical Analysis; User-Friendly Interface; Improved Data Import/Export <br>
                  Stay tuned for future updates as we continue to enhance our app and provide a more efficient, reliable, and user-friendly qPCR data analysis experience.")
             
             ),
    # ** 2.3.2 数据处理 ----
     tabItem(
      tabName = "qpcr_data_process",
      fluidRow(column(
        6,
        box(
          title = "Data upload",
          color = "primary",
          width = 12,
          height = 230,
          fluidRow(column(
            width = 6,
            fileInput(
              inputId = "file",
              label = "Upload file",
              placeholder = ".xls or .xlsx format",
              multiple = F,
              accept = c(".xls", ".xlsx")
            )
          ),
          
          column(
            width = 6,
            pickerInput(
              inputId = "Id088",
              label = "Relative quantification method",
              choices = c("delta delta Ct", "delta Ct", "standard curve"),
              options = list(style = "btn-secondary")
            )
          )),
         
          
          div(
            width = 12, 
            fluidRow(
              column(width = 10,p("Please upload your data in the format provided in the demo data. Do not delete any sheets, modify the column names, or change the sheet names."),),
              column(width = 2,
                     downloadBttn(outputId = "downloaddemo",
                         label = "Demo data", color = "primary",
                         style = "unite", size = "xs")
                     )
            )
            
            
          )
          
        ),
      ),
      
      column(
        width = 6,
        conditionalPanel(condition = "output.isfileup",
                         uiOutput("dynamic_box"))
      )),
      
      
      conditionalPanel(
        condition = "output.isfileup",
        tabBox(
          id = "tabcard",
          title = "Result Analysis and Visualization",
          selected = "Relative expression",
          width = 12,
          status = "primary",
          solidHeader = FALSE,
          type = "tabs",
          tabPanel(
            title = "Relative expression",
            fluidRow(
              column(width = 4, plotlyOutput("plot_relExp")),
              column(width = 4, dataTableOutput("table_rel")),
              column(width = 4, dataTableOutput("table_pVal"))
            ),
            dropdownButton(
              tags$h3("Plot parameters"),
              
              selectInput(
                inputId = 'plot_test_method',
                label = 'Statistical testing methods',
                choices = c("t.test", "wilcox.test"),
                selected = "t.test"
              ),
              
              selectInput(
                inputId = 'plot_p_label',
                label = 'P-value notation',
                choices = c("p.format", "p.signif"),
                selected = "p.signif"
              ),
              
              circle = F,
              status = "danger",
              icon = icon("gear"),
              width = "30px",
              
              tooltip = tooltipOptions(title = "Click to see inputs !")
            )
          ),
          tabPanel(title = "Melting curve",
                   fluidRow(
                     column(width = 2, uiOutput("dynamic_plot_melt")),
                     column(width = 6, plotlyOutput("plot_melt")),
                     column(width = 4, dataTableOutput("table_sample"))
                   )),
          tabPanel(title = "Standard curve",
                   fluidRow(
                     column(width = 6, plotlyOutput("plot_standard")),
                     column(width = 6, dataTableOutput("table_standard"))
                   ))
        )
      )
      ),
    # ** 2.3.3 qPCR原理 ----
     tabItem(tabName = "qpcr_basic",
             tabBox(
               id = "tabcard",
               title = "qPCR Principles",
               width = 12,
               selected = "Introduction",
               status = "primary",
               solidHeader = FALSE,
               type = "tabs",
               tabPanel(
                 title = "Introduction",
                 p("qPCR (quantitative PCR), also known as real-time PCR, is a powerful technique in the field of molecular biology. It involves the amplification of specific sequences within DNA or cDNA templates by using specific primers, heat-stable DNA polymerase, and thermal cycling, resulting in the replication of these sequences thousands to millions of times. Unlike traditional PCR, qPCR allows real-time measurement of PCR product amplification during the exponential phase of the reaction, enabling the precise determination of the initial quantity of target molecules with high accuracy. This is achieved by monitoring the increase in fluorescence signal, which is proportional to the quantity of PCR products."),
                 hr(),
                 strong("Applications of qPCR:"),
                 br(),
                 p("🔹 Gene expression analysis (most common): It determines the level of gene expression in different samples by assessing the relative abundance of transcripts."),
                 p("🔹 Viral titer determination: Quantifies the number of virus copies in a sample. Standard curves are generated using known genome equivalents or nucleic acid from titrated virus control materials and compared to the standard curve for titer determination."),
                 p("🔹 Copy number variation analysis: Analyzes genomic duplications or deletions."),
                 p("🔹 Allelic gene typing: Unlike the methods mentioned above, this requires endpoint fluorescence to determine SNP genotypes. Primer and probe design are crucial to ensure low cross-reactivity between allele-specific reactions."),
                 hr(),
                 strong("Advantages of qPCR:"),
                 br(),
                 p("🔹 It can monitor the progress of each cycle in the amplification process in real-time."),
                 p("🔹 It can accurately determine the quantity of amplification products in each cycle, allowing for precise quantification of the starting material in the sample."),
                 p("🔹 It has a wider detection dynamic range."),
                 p("🔹 Amplification and detection are completed in the same tube, eliminating the need for subsequent PCR processing steps."),
                 
               ),
               tabPanel(
                 title = "Chemical Principles",
                 p("There are two common reagents for fluorescence quantification: TaqMan and SYBR Green."),
                 h4("1. TaqMan"),
                 p("TaqMan probes are designed based on the principle of Fluorescence Resonance Energy Transfer (FRET)."),
                 strong("Components:"),
                 p("🔹  Reporter: When excited by light, it generates a fluorescent signal."),
                 p("🔹  Quencher: When in close proximity to the fluorescent substance, it quenches (eliminates) the fluorescent signal."),
                 p("🔹  Probe:A specific DNA sequence that can hybridize with the target DNA fragment."),
                 strong("Working Principle:"),
                 p("🔹  At the beginning of the qPCR experiment, the TaqMan probe hybridizes with the target DNA. At this point, the reporter and quencher are close, and the fluorescence is quenched, making it undetectable."),
                 p("🔹  As PCR proceeds, DNA polymerase synthesizes a new DNA strand from the 5' end to the 3' end. When the polymerase encounters the TaqMan probe, its 5'-3' exonuclease activity cleaves the reporter. At this point, the reporter and quencher separate, and the fluorescent substance is excited, generating a detectable signal."),
                 p("🔹  In each PCR cycle, new reporters are cleaved, producing new fluorescent signals. Therefore, the increase in the detected fluorescent signal is directly proportional to the quantity of the target DNA, enabling quantitative DNA detection."),
                 p("One common fluorescent reporter-quencher pair is FAM (emits green fluorescence) and Black Hole Quencher 1."),
                 strong("Advantages:"),
                 p("🔹  High specificity: TaqMan probes can be highly specific because they can be designed to be fully complementary to the target sequence, reducing the occurrence of false-positive results."),
                 p("🔹  Accurate quantification: TaqMan probes can provide accurate quantification results because the generation of the fluorescence signal is directly related to the extent of target DNA amplification."),
                 p("🔹  Multiplexing capability: Using TaqMan probes labeled with different fluorescent substances, multiple target sequences can be simultaneously detected in a single reaction."),
                 strong("Disadvantages:"),
                 p("🔹  High cost: The synthesis of TaqMan probes involves fluorescent substances and quenchers, making it relatively expensive."),
                 p("🔹  Complex design: TaqMan probes need to be specifically designed for each target sequence."),
                 p("🔹  High optimization requirements: Experiment conditions can significantly impact probe binding and fluorescence production, requiring careful optimization for the best results."),
                 hr(),
                 h4("2. SYBR Green"),
                 strong("These are the most commonly used qPCR reagents."),
                 strong("Working Principle:"),
                 p("SYBR Green is a fluorescent dye that can bind to double-stranded DNA and emit a fluorescence signal. In qPCR reactions, SYBR Green binds to the double-stranded DNA in the amplification product, and the intensity of the fluorescence signal is proportional to the quantity of the amplification product."),
                 strong("Advantages:"),
                 p("🔹  Easy to use, requiring the design of forward and reverse primers only, with no need for specific probe design."),
                 p("🔹  Lower cost compared to TaqMan probes."),
                 p("🔹  Flexibility: SYBR Green can be more easily used for different target genes because you only need to change the primer sequences."),
                 strong("Disadvantages:"),
                 p("🔹  Lack of specificity, as it can bind to any double-stranded DNA. This means that in qPCR reactions, if there are non-specific amplification products (such as primer dimers), SYBR Green will also bind and emit a fluorescence signal, affecting the accuracy of quantitative results."),
                 p("🔹  Difficulty in detecting multiple PCR products in one reaction due to the non-specific binding of SYBR Green."),
                 hr(),
                 h4("Others"),
                 p("For example, common qPCR reagent chemistries include Molecular Beacons, Hybridization probes, Eclipse probes, Amplifluor chemistry, Scorpions primers, LUX primers, and BD QZyme."),
                 a("You can click here to access more information.", href = "https://www.bio-rad.com/en-us/applications-technologies/qpcr-assay-development?ID=5a762865-ad72-dbfa-2c23-30243bed845a", target='_blank')
                 
               ),
               tabPanel(
                 title = "Terminology",
                 div(
                   id = "center-image",
                   imageOutput(outputId = "img1",
                               width = "auto", height = "auto")
                 ),
                 HTML("<p>1. <b>Baseline</b>: The baseline in qPCR refers to the initial cycles of the amplification reaction where the fluorescence signal changes minimally. At this stage, the reaction has not yet entered the exponential amplification phase, and the fluorescence signal primarily arises from background noise. Setting the baseline is used to reduce the impact of background noise on data analysis.</p>"),
                 HTML("<p>2. <b>Threshold</b>: The threshold is a set fluorescence signal level in qPCR experiments, typically set above the baseline and during the exponential phase of PCR, where rapid amplification occurs. The threshold is often set as ten times the standard deviation of the baseline. The threshold is established to distinguish the relevant amplification signal from the background and is usually automatically set by the qPCR program.</p>"),
                 HTML("<p>3. <b>Ct (Cycle threshold)</b>: Ct refers to the number of PCR cycles it takes for the fluorescence signal to reach the set threshold. The Ct value has a linear relationship with the logarithm of the initial template copy number. A lower Ct value corresponds to a higher initial template concentration, while a higher Ct value corresponds to a lower initial template concentration.</p>"),
                 HTML("<p>4. <b>Standard Curve</b>: A standard curve is a method used for quantitative PCR (qPCR). It is established by using known concentrations of template DNA in the PCR reaction. By creating serial dilutions of template DNA and performing real-time quantitative PCR, a standard curve is generated. Then, with the logarithm of the initial template DNA concentration on the X-axis and Ct on the Y-axis, a linear regression analysis is performed to fit a line that accurately describes the relationship between concentration and Ct values. The standard curve also includes an R^2 value to measure the repeatability of duplicate samples. Repeating the standard curve is done to assess consistency, ensuring the accuracy of sample data.</p>"),
                 HTML("<p>5. <b>Slope</b>: The slope of the standard curve reflects the efficiency of the PCR reaction. Ideally, the amplification in each cycle should be approximately 2-fold, corresponding to a slope close to -3.32. If the slope deviates from this value, there may be issues with the PCR reaction, such as low amplification efficiency or non-specific amplification. The slope of the curve can be used to determine the reaction efficiency.</p>"),
                 HTML("<p>6. <b>Efficiency</b>: The efficiency of the PCR reaction, calculated from the slope of the standard curve, is expressed as a percentage. The formula for calculating efficiency is E = 10^(-1/slope), Efficiency(%) = (E-1)*100. An efficiency of 100% indicates perfect doubling of the template in each cycle, but an acceptable range for analysis validation is 90-110%. This efficiency range corresponds to a standard curve slope of -3.58 to -3.10.</p>"),
                 HTML("<p>7. <b>Reference Fluorescent Dye</b>: ROX is commonly used as a reference fluorescent dye in qPCR. It is an inert dye used for data normalization in qPCR experiments. During qPCR experiments, variations in fluorescence signal intensity may occur between different samples due to experimental procedures, instrument settings, or other factors. By using ROX as a reference dye, the fluorescence signal of reporter dyes (such as SYBR Green) can be corrected, reducing non-specific fluorescence variations between experimental results. ROX reference fluorescent dye has the following characteristics: <br> 🔹  Its fluorescence intensity remains constant during the experiment, unaffected by the PCR reaction, making it a stable reference signal.<br> 🔹  The fluorescence emission wavelength of ROX dye is different from other fluorescent dyes commonly used in qPCR, such as SYBR Green and TaqMan probes, preventing their signals from interfering with each other.</p>"),
                 HTML("<p>8. <b>Reference Gene</b>: A reference gene is a gene whose expression remains relatively stable across all samples, aiming to control variability in experiments and normalize the expression levels of target genes. In qPCR experiments, factors such as the initial concentration of samples, RNA quality, and reverse transcription efficiency can affect the PCR results. A stably expressed reference gene is needed to eliminate these influences. By comparing the qPCR results of the target gene and the reference gene, you can obtain the relative expression level of the target gene. This allows for an accurate comparison of the expression level of the target gene in different samples, even in the presence of unavoidable variations in experimental conditions or sample handling processes.</p>")
               ),
               tabPanel(
                 title = "Quantification methods",
                 HTML("<p>1. <b>Absolute Quantification</b>: Absolute quantification involves the amplification of known quantities of samples through serial dilution to generate a standard curve. Then, unknown samples are quantified by comparing their results to this curve. Absolute quantification is commonly used in areas such as viral quantification and gene copy number analysis.</p>"),
                 HTML("<p>2. <b>Relative Quantification</b>: Relative quantification compares the relative expression levels of a target gene in different samples, rather than directly measuring the absolute quantity of the target gene. By comparing the expression of the target gene in the target sample to the expression of the same gene in a control sample, you obtain the relative proportion of the target gene in each sample. Relative quantification is often used to assess differences in gene expression levels between different samples.<br>Common methods of relative quantification include: <br>🔹 Delta delta Ct(Most commonly used) <br>🔹 Delta Ct<br>🔹 Standard curve method </p>"),
                 actionBttn(
                   inputId = "Id113",
                   label = "You can find more information in the \"Analysis Instructions\".", 
                   style = "stretch",
                   color = "primary"
                 )
                 )
              )
             ),
    
    # ** 2.3.4 qPCR实验设计 ----
     tabItem(tabName = "qpcr_exp_design",
             HTML("<h3>1. RNA Extraction</h3>RNA Extraction is the first step in the experiment, the purpose of which is to isolate RNA from biological samples (such as cells, tissues, etc.). This process usually uses specialized RNA extraction kits, such as TRIzol, etc. These kits break cells by physical and chemical methods, and then separate and purify RNA through phenol chloroform separation and ethanol precipitation. After RNA extraction, it is usually necessary to evaluate the quality and concentration of RNA to ensure it is suitable for subsequent experiments."),
             hr(),
             HTML("<h3>2. Reverse Transcription</h3>The process of transcribing RNA into cDNA. This step uses reverse transcriptase to synthesize cDNA using RNA as a template. During the reverse transcription process, reverse transcription primers (ie. random primers or oligo dT primers) and reverse transcriptase are usually added. In this reaction, with RNA as a template and the guide of the primer, the reverse transcriptase synthesizes complementary cDNA for use as a template for subsequent qPCR analysis."),
             hr(),
             HTML("<h3>3. Preparations for qPCR</h3>
                   <h4>3.1 Primer Design</h4>
                   <b>Principles:</b><br>
                   1. Primer Amplification Fragment: The length of the amplified fragment should be about 50-150 bp, because longer products cannot be efficiently amplified.<br>
                   2. Primer Length: The length of the primer should generally be 18-24 nucleotides.<br>
                   3. Primer Specificity: The primer should be specific to the target gene sequence. To ensure the specificity of the primer, a BLAST search can be conducted on the public database to ensure the primer only recognizes the target. Primer specificity can also be initially analyzed by gel electrophoresis.<br>
                   4. The Tm value of the primer: The primer should have a similar melting temperature (within 5°C), and the GC content is about 50%.<br>
                   5. Complementarity and Hybridization between Primers: Analyze the primer pair sequence to avoid complementarity and hybridization between primers, to avoid the formation of primer dimers.<br>
                   6. Introns on both sides of the primer: Primers designed across introns can distinguish cDNA and potential contamination genomic DNA amplification.<br>
                   <b>Primer Design Software/Websites:</b><br>
                   1. <a href='http://www.oligo.net' target='_blank'>Oligo 7</a><br>
                   2. <a href ='https://primer3.ut.ee/' target='_blank'>Primer 3</a><br>
                   3. <a href = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/' target='_blank'>NCBI primer blast</a><br>
                   4. <a href = 'https://pga.mgh.harvard.edu/primerbank/' target='_blank'>Primer Bank</a><br>
                   <br>
                   <h4>3.2 qPCR Reagent Selection</h4>
                   If your goal is to quickly and economically perform qPCR experiments, and you don't need to perform multiplex qPCR, then SYBR Green may be a good choice. If you need to perform multiplex qPCR or require higher specificity and accuracy, then TaqMan may be more suitable. In general, SYBR Green is used to compare differences in target gene expression between samples.<br>
                   <br>
                   <h4>3.3 Method of Quantification Selection</h4>
                   🔹  Absolute quantification can determine the actual copy number of the target, but it is also the most laborious and complex form of quantification. This method requires meticulous planning and highly accurate standard curves. Absolute quantification is often used to determine virus titers.<br>
                   🔹  Relative quantification also requires careful planning, but the generated data is relative abundance, rather than exact copy number. This method is suitable for gene expression studies, the main quantification schemes: ΔΔCt, ΔCt, and standard curve quantification.<br>
                   <br>
                   <h4>3.4 Reference Gene Selection</h4>
                   Eliminate experimental inconsistencies such as the quality and quantity of the initial template, PCR amplification efficiency, etc. it is extremely important for qPCR experimental design. Using standard primed genes (also known as reference genes) is the most thorough method to solve almost all differences in real-time fluorescence quantitative PCR. However, when using this method, the expression of this gene must be consistent in all samples.<br>
                   <b>Common Endogenous Standards Include:</b><br>
                   🔹  β-actin: cytoskeletal component<br>
                   🔹  18S Ribosomal RNA (rRNA): Ribosomal subunit<br>
                   🔹  Cyclophilin A (PPIA): Serine/threonine phosphatase inhibitor<br>
                   🔹  Glyceraldehyde 3-phosphate dehydrogenase (GAPDH): Glycolysis pathway<br>
                   🔹  β-2 Microglobulin (B2M): Major histocompatibility complex<br>
                   🔹  β-glucuronidase (GUSB): An exogenous glycosidase in lysosomes<br>
                   🔹  Hypoxanthine phosphoribosyltransferase (HPRT1): Purine salvage synthesis pathway<br>
                   🔹  TATA box binding protein (TBP): RNA transcription.<br>
                   <br>
                   💙 You can conveniently use the <a href = 'https://ngdc.cncb.ac.cn/icg/' target='_blank'> ICG database </a> to obtain primers for reference genes of various species.<br>
                   <em> ICG is a curated knowledgebase of internal control genes (or reference genes) for RT-qPCR normalization in a variety of species. Based on literature curation, ICG provides a comprehensive collection of high-quality experimentally verified internal control genes and their application scenarios for both model and non-model organisms. </em> <br>
                   <br>
                   <b>Exogenous Standards</b><br>
                   If a highly consistent endogenous standard cannot be found for a specific group of samples, then exogenous standards are a viable option. Exogenous reference genes are synthesized or in vitro transcribed RNA, whose sequences do not exist in experimental samples. Given its exogenous nature, it will not show normal biological fluctuations in cells in different states or treatments."),
             hr(),
             HTML("<h3>4. qPCR Amplification Quantification</h3>
                  In this step, cDNA is used as a template for PCR amplification, and the change in fluorescence signal is monitored in real time for quantitative analysis. During the experiment, the fluorescence signal enhances after each cycle of amplification, reflecting the increase in amplification product. The amplification process can use methods such as SYBR Green or TaqMan for fluorescence labeling and detection. By comparing the strength of the fluorescence signals of various samples, the expression level of the target gene in various samples can be obtained.")
             
             ),
    
    # ** 2.3.5 数据处理说明 ----
     tabItem(tabName = "data_process_info",
             h3("1. Absolute Quantification"),
             p("🔹 Create a standard curve by plotting the logarithm of the known template concentration obtained from serial dilutions on the x-axis and the measured Ct values on the y-axis. "),
             p("🔹 Fit the standard curve using linear regression with the formula y = kx + b, where y represents the Ct value, x is the logarithm of the concentration of the test sample, k is the slope, and b is the intercept."),
             p("🔹 Determine the Ct values of unknown samples and then use the regression equation to calculate the absolute quantification of the target gene in these unknown samples. The specific formula is: Sample = 10^(Ct - b) / k."),
             div(
               id = "center-image",
               imageOutput(outputId = "img_abs",
                           width = "auto", height = "auto")
             ),
             hr(),
             
             h3("2. Relative Quantification"),
             h4("2.1 ΔCt"),
             p("The most basic form of relative quantification, this method is relatively simple."),
             strong("Formula:"),
             br(),
             HTML("ΔCt = Ct_treated - Ct_control<br> 
                  Relative expression = 2^ΔCt"),
             p("Ct represents the cycle threshold, which is the number of cycles required to reach a specific fluorescence threshold. The difference between the Ct of the experimental group and the Ct of the control group, denoted as ΔCt, can be used to represent the relative expression of the experimental group compared to the control group."),
             HTML("In the example shown in the figure, the fold difference is calculated as 2^ΔCt: <br>
                  ΔCt = 12 - 6 = 6 <br> 
                  Fold change = 2^ΔCt = 2^6 = 64"),
             div(
               id = "center-image",
               imageOutput(outputId = "img_dCt",
                           width = "auto", height = "auto")
             ),
             hr(),
             h4("2.2 ΔΔCt"),
             p("The ΔΔCt method improves the accuracy of relative quantification by correcting for differences in amplification efficiency of different genes based on the foundation of the ΔCt method, and this correction is achieved by using reference genes."),
             strong("Formula:"),
             br(),
             HTML("∆Ct = Ct_GOI – Ct_HKG <br> 
                  ∆∆Ct = ∆Ct_Treated – ∆Ct_Control <br> 
                  Relative expression = 2^(-ΔΔCt) <br> 
                  You can refer to the table in the figure for the calculation method."),
             div(
               id = "center-image",
               imageOutput(outputId = "img_ddCt",
                           width = "auto", height = "auto")
             ),
             hr(),
             h4("2.3 Standard Curve"),
             p("The standard curve method, also known as the Pfaffl method, is a commonly used relative quantification method in qPCR experiments. Unlike the ΔΔCt method, the Pfaffl method can better handle situations where the amplification efficiency is not 100%."),
             strong("Formula:"),
             br(),
             HTML("Relative expression = E_GOI^ΔCt_GOI / E_HKG^ΔCt_HKG <br>
                  E_GOI and E_HKG represent the amplification efficiency of the target gene (GOI) and the reference gene (HKG), respectively. These efficiencies can be obtained from the slope of the standard curve.<br>
                  ΔCt_GOI and ΔCt_HKG represent the difference in Ct values between the target gene and the reference gene in the experimental group and the control group, respectively."),
             strong("Steps:"),
             br(),
             HTML("🔹 Choose a stable reference gene to normalize the RNA or cDNA quantity in the experimental group and the control group.<br>
                  🔹 Determine the amplification efficiency of the target gene and the reference gene through the standard curve method or other methods.<br>
                  🔹 Perform qPCR amplification with RNA or cDNA from the experimental group and the control group and measure the Ct values of the target gene and the reference gene.<br>
                  🔹 Calculate ΔCt_GOI and ΔCt_HKG.<br>
                  🔹 Calculate the relative expression level. <br>
                  You can refer to the table in the figure for the calculation method."),
             div(
               id = "center-image",
               imageOutput(outputId = "img_sc",
                           width = "auto", height = "auto")
             ),
             hr(),
             HTML("💙 You can click here to access more information about 
                  <a href = 'https://toptipbio.com/delta-delta-ct-pcr/', target='_blank'>ΔΔCt </a> and 
                  <a href = 'https://toptipbio.com/pfaffl-method-qpcr/', target='_blank'>Pfaffl </a> method."),
             ),
    
    # ** 2.3.6 帮助 ----
     tabItem(tabName = "help",
             HTML("💙 You can click 
             <a href = 'https://github.com/Xiang-yuZHAO/qPCR', target='_blank'> here </a> 
             to get the open source code."),
             hr(),
             includeHTML(path = "www/md/help.html"),
             hr(),
             HTML("💙 You can find more assistance in the following PDF (Chinese version), or you can click <a href = 'https://rise.articulate.com/share/CzJ69zEmF4MgocDs26iDJscAyM-xxAt-#/', target='_blank'>here </a> to access the English version of the help."),
             br(),
             tags$iframe(style = "height:900px; width:100%; scrolling=yes", src= "ThermoFisher_qPCR.pdf")
             )
    )
  ), 
  
  # * 2.4 footer ----
  footer = dashboardFooter(left = "qPCR shiny version 1.0",
                           right = "@ZXY")
  
  
)


#### 3. Server ####
server <- function(input, output, session){
  observeEvent(input$Id112,{
    updateTabItems(session = session,
                   inputId = "updateTabItems_sidebarmenu",
                   selected = "qpcr_data_process")
  })
  
  observeEvent(input$Id113,{
    updateTabItems(session = session,
                   inputId = "updateTabItems_sidebarmenu",
                   selected = "data_process_info")
  })
  
  # 3.1 输入 ----
  inputdf <- reactive({
    req(input$file)  
    inFile <- input$file
    df <- read_excel(inFile$datapath,sheet = 1) 
    df
  })
  
  inputdf_melt <- reactive({
    req(input$file)  
    inFile <- input$file
    df <- read_excel(inFile$datapath,sheet = 2) 
    df
  })
  
  
  inputdf_standard <- reactive({
    req(input$file) 
    inFile <- input$file
    df <- read_excel(inFile$datapath,sheet = 3) 
    df
  })
  
  inputdf_sample <- reactive({
    req(input$file) 
    inFile <- input$file
    df <- read_excel(inFile$datapath,sheet = 4) 
    df
  })
  
  
  
  # * 3.1.1 conditional panel ----
  output$isfileup <- reactive({
    return(!is.null(inputdf()))
  })
  outputOptions(output, 'isfileup', suspendWhenHidden = FALSE)
  
  
  rv <- reactiveValues()
  # 3.2 输入文件判断 ----
  observe({
    pre_hash <- c("相对定量" = "6e582af93e53fe8c3eb4f71633f27211", 
                  "熔解曲线" = "c643fc7e1e0d83005bd3b95144ed3f79", 
                  "标准曲线" = "67aa3d33dc98845f7ad07d8dbf9bbaa2", 
                  "样本信息" = "855e15d492dd333d101bd1654f1491fb")
    pre_hash_Eng <- c("Ct value", "Melting curve", "Standard curve", "Sample info")
    
    changed_sheets <- c()
    changed_sheets[1] <- digest::digest(object = inputdf(), algo = "md5")
    changed_sheets[2] <- digest::digest(object = inputdf_melt(), algo = "md5")
    changed_sheets[3] <- digest::digest(object = inputdf_standard(), algo = "md5")
    changed_sheets[4] <- digest::digest(object = inputdf_sample(), algo = "md5")
    
    diff <- setdiff(changed_sheets, pre_hash)
    file_upload <- names(pre_hash)[which(changed_sheets %in% diff)]
    
    isRel <- ifelse("相对定量" %in% file_upload, T, F)
    isMelt <- ifelse("熔解曲线" %in% file_upload, T, F)
    isCurve <- ifelse("标准曲线" %in% file_upload, T, F)
    isSample <- ifelse("样本信息" %in% file_upload, T, F)
    
    rv$isRel <- isRel
    rv$isCurve <- isCurve
    rv$isMelt <- isMelt
    rv$isSample <- isSample
    
    file_upload_confirm <- c(ifelse(isRel, " ✅ (User data has been uploaded)", " ❎(Demo data is being used)"),
                             ifelse(isMelt, " ✅  (User data has been uploaded)", "  ❎(Demo data is being used)"),
                             ifelse(isCurve, " ✅  (User data has been uploaded)", "  ❎(Demo data is being used)"),
                             ifelse(isSample, " ✅  (User data has been uploaded)", "  ❎(Demo data is being used)")
                             )
    
    
    shinyalert(title = "Upload successful", 
               text = paste(pre_hash_Eng, file_upload_confirm,collapse = "\n") , 
               type = "success",
               closeOnClickOutside = T,
               showCancelButton = F)
    
  })
  
  # 3.3 示例数据 ----
  output$downloaddemo <- downloadHandler(
    filename = function() {
      "qPCR_demo_data.xlsx"
    },
    content = function(file) {
      file.copy("www/demo/qPCR demo.xlsx", file) 
    }
  )  
  
  # 3.4 动态UI 分析方法选择 ----
  output$dynamic_box <- renderUI({
    selected_method <- input$Id088
    if (selected_method == "delta delta Ct") {
      box(
        title = "Delta Delta Ct",
        color = "success",
        width = 12,
        height = 230,
        fluidRow(
          column(width = 6,
                 uiOutput('Control'),uiOutput('HKG')),
          column(width = 6, 
                 uiOutput('Treated'),uiOutput('GOI'))
        )
        
      )
    } else if (selected_method == "delta Ct") {
      box(
        title = "Delta Ct",
        color = "success",
        width = 12,
        height = 230,
        fluidRow(
          column(width = 6,
                 uiOutput('Control1'),uiOutput('Treated1')),
          column(width = 6, 
                 uiOutput('gene1'))
        )
      )
    } else if (selected_method == "standard curve") {
      box(
        title = "Standard curve",
        color = "success",
        width = 12,
        height = 230,
        fluidRow(
          column(width = 6,
                 uiOutput('Control2'),uiOutput('HKG2')),
          column(width = 6, 
                 uiOutput('Treated2'),uiOutput('GOI2'))
        )
        
      )
    }
  })
  
  ##获取动态UI的名称
  
  var_sample <- reactive({
    getName(df = inputdf(),object = "Sample")
  })
  
  var_gene <- reactive({
    getName(df = inputdf(),object = "Gene")
  })
  
  # 获取ddct的动态UI
  output$Control <- renderUI({
    selectInput('Control', 'Control', 
                choices = var_sample())
  })
  
  output$Treated <- renderUI({
    req(input$Control)
    selectInput('Treated', 'Treated',
                choices = setdiff(var_sample(), input$Control),
                multiple = T,
                selected = setdiff(var_sample(), input$Control)[1])
  })
  
  
  output$HKG <- renderUI({
    selectInput('HKG', 'Reference gene',
                choices = var_gene())
  })
  
  output$GOI <- renderUI({
    req(input$HKG)
    selectInput('GOI', 'Target gene',
                choices = setdiff(var_gene(), input$HKG))
  })
  
  # 获取dct的动态UI
  output$Control1 <- renderUI({
    selectInput('Control1', 'Control', 
                choices = var_sample())
  })
  
  output$Treated1 <- renderUI({
    req(input$Control1)
    selectInput('Treated1', 'Treated',
                choices = setdiff(var_sample(), input$Control1),
                multiple = T,
                selected = setdiff(var_sample(), input$Control1)[1])
  })
  
  output$gene1 <- renderUI({
    selectInput('gene1', 'Target gene',
                choices = var_gene())
  })
  
  # 获取标准曲线的动态UI
  output$Control2 <- renderUI({
    #if (rv$isCurve) {
      selectInput('Control2', 'Control', 
                choices = var_sample())
    #}
    
  })
  
  output$Treated2 <- renderUI({
    
    req(input$Control2)
    selectInput('Treated2', 'Treated',
                choices = setdiff(var_sample(), input$Control2),
                multiple = T,
                selected = setdiff(var_sample(), input$Control2)[1])
  })
  
  
  output$HKG2 <- renderUI({
    
    selectInput('HKG2', 'Reference gene',
                choices = var_gene())
  })
  
  output$GOI2 <- renderUI({
   
    req(input$HKG2)
    selectInput('GOI2', 'Target gene',
                choices = setdiff(var_gene(), input$HKG2))
  })
  
  # 3.5 动态UI 熔解曲线sample和gene选择 ----
  output$dynamic_plot_melt <- renderUI({
    #if (rv$isMelt) {
      box(
        title = "Plot parameters",
        color = "success",
        collapsible = F,
        width = 12,
        height = 300,
        fluidRow(column(12, uiOutput('melt_sample')),
                 column(12, uiOutput('melt_gene')))
      )
    #}
    
  })
  
  var_melt_sample <- reactive({
    unique(inputdf_sample()[["Sample"]])
  })
  
  var_melt_gene <- reactive({
    unique(inputdf_sample()[["Gene"]])
  })
  
  output$melt_sample <- renderUI({
    selectInput('melt_sample', 'Sample',
                choices = var_melt_sample())
  })
  
  output$melt_gene <- renderUI({
    selectInput('melt_gene', 'Gene',
                choices = var_melt_gene())
  })
  
  
  # 3.6 计算相对表达 ----
  relExpLevel <- reactive({
    selected_method <- input$Id088
    if (selected_method == "delta delta Ct") {
      req(input$GOI)
      req(input$Treated)
      pcr_analyze(df = inputdf(), method = 'delta_delta_ct',
                  Treated_samples = input$Treated,
                  Control = input$Control,
                  HKG = input$HKG, GOI = input$GOI)
      
    } else if (selected_method == "standard curve") {
      req(input$GOI2)
      req(input$Treated2)
      req(inputdf_standard())
      pcr_analyze(df = inputdf(), method = 'relative_curve', 
                  df_standard = inputdf_standard(),
                  Control = input$Control2, Treated_samples = input$Treated2,
                  HKG = input$HKG2, GOI = input$GOI2)
      
    } else if (selected_method == "delta Ct") {
      req(input$gene1)
      req(input$Treated1)
      pcr_analyze(df = inputdf(), method = 'delta_ct', 
                  Treated_samples = input$Treated1,
                  Control = input$Control1, gene = input$gene1)
      
    }
  })
  
  # 3.7 计算两两之间的p值 ----
  pVal_pos <- reactive({
    req(relExpLevel())
    get_P_pos(relExpLevel())
  })
  
  
  # 3.8 绘图展示相对表达 ----
  output$plot_relExp <- renderPlotly({
    req(relExpLevel())
    selected_method <- input$Id088
    if (selected_method == "delta delta Ct") {
      req(input$Treated)
      req(input$GOI)
      req(input$HKG)
      tryCatch({
        if(rv$isRel){
           plot <- pcr_plot_analyze(
            df_plot = relExpLevel(),
            label_pVal_y = pVal_pos(),
            ref_group = input$Control,
            p_label = input$plot_p_label,
            test_method = input$plot_test_method ) + 
          ggtitle("Reletive expression level")
        }else {
          plot <- pcr_plot_analyze(
            df_plot = relExpLevel(),
            label_pVal_y = pVal_pos(),
            ref_group = input$Control,
            p_label = input$plot_p_label,
            test_method = input$plot_test_method ) + 
            ggtitle("This is demo data") +
            theme(plot.title = element_text(colour = "red"))
        }
        ggplotly(plot)
      }, error = function(e) {
        message("An error occurred: ", e$message)
        return(NULL)
      })
    } else if (selected_method == "standard curve") {
      req(input$Control2)
      req(input$Treated2)
      req(input$HKG2)
      req(input$GOI2)
      tryCatch({
        
        if(rv$isRel){
          plot <- pcr_plot_analyze(
            df_plot = relExpLevel(),
            ref_group = input$Control2,
            label_pVal_y = pVal_pos(),
            p_label = input$plot_p_label,
            test_method = input$plot_test_method ) +
          ggtitle("Reletive expression level")
        }else{
          plot <- pcr_plot_analyze(
            df_plot = relExpLevel(),
            ref_group = input$Control2,
            label_pVal_y = pVal_pos(),
            p_label = input$plot_p_label,
            test_method = input$plot_test_method ) +
            ggtitle("This is demo data") +
            theme(plot.title = element_text(colour = "red"))
        }
        ggplotly(plot)
      }, error = function(e) {
        message("An error occurred: ", e$message)
        return(NULL)
      })
      
    } else if (selected_method == "delta Ct") {
      
      tryCatch({
        if(rv$isRel){
          plot <- pcr_plot_analyze(
            df_plot = relExpLevel(),
            ref_group = input$Control1,
            label_pVal_y = pVal_pos(),
            p_label = input$plot_p_label,
            test_method = input$plot_test_method
          ) +
            ggtitle("Reletive expression level")
        }else{
          plot <- pcr_plot_analyze(
            df_plot = relExpLevel(),
            ref_group = input$Control1,
            label_pVal_y = pVal_pos(),
            p_label = input$plot_p_label,
            test_method = input$plot_test_method
          ) +
            ggtitle("This is demo data") +
            theme(plot.title = element_text(colour = "red"))
        }
        ggplotly(plot)
      }, error = function(e) {
        message("An error occurred: ", e$message)
        return(NULL)
      })
    }
    
    
    
  })
  
  # 3.9 表格展示相对表达 ----
  output$table_rel <- DT::renderDataTable({
    req(relExpLevel())
    #if(rv$isRel){
      df <- relExpLevel()
      df$`Relative expression`  <- format(df$`Relative expression`, scientific =F, digit = 4)
      colnames(df)[3] <- "Rel. Exp."
      df
    #}
    
  },
  extensions = c('Buttons', 'Responsive','Scroller'),
  options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'excel'),
    deferRender = TRUE,
    scrollY = 300,
    scroller = TRUE,
    searching = FALSE
  ))
  
  # 3.10 表格展示两两之间p值 ----
  output$table_pVal <- DT::renderDataTable({
    req(relExpLevel())
    #if(rv$isRel){
      df <-  pcr_test(df_plot = relExpLevel(),method = input$plot_test_method)
      df$p_value <- format(df$p_value, scientific =T,digit = 4)
      df
    #}
    
  }, 
  extensions = c('Buttons', 'Responsive','Scroller'),
  options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'excel'),
    deferRender = TRUE,
    scrollY = 300,
    scroller = TRUE,
    searching = FALSE
  ))
  
  # 3.11 绘图展示标准曲线 ----
  output$plot_standard <- renderPlotly({
    if(rv$isCurve){
      tryCatch({
        plot <- pcr_standard(df = inputdf_standard(),plot = T) +
          ggtitle("Standard curve")
        ggplotly(plot)
      }, error = function(e) {
        message("An error occurred: ", e$message)
        return(NULL)
      })
    } else {
      tryCatch({
        plot <- pcr_standard(df = inputdf_standard(),plot = T) +
          ggtitle("This is demo data") + 
          theme(plot.title = element_text(color = "red"))
        ggplotly(plot)
      }, error = function(e) {
        message("An error occurred: ", e$message)
        return(NULL)
      })
    }
  })
  
  # 3.12 表格展示标准曲线 ----
  output$table_standard <- DT::renderDataTable({
    #if(rv$isCurve){
      df <- pcr_standard(df = inputdf_standard(),plot = F)
      df_rowNames <- row.names(df)
      df <- as.data.frame(lapply(df, format_scientific))
      colnames(df)[5] <- "% Efficiency"
      row.names(df) <- df_rowNames
      df
    #}
    
  }, 
  extensions = c('Buttons', 'Responsive','Scroller'),
  options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'excel'),
    deferRender = TRUE,
    scrollY = 300,
    scroller = TRUE,
    searching = FALSE
  ))
  
  # 3.13 表格展示熔解曲线样本 ----
  output$table_sample <- DT::renderDataTable({
    #if(rv$isMelt){
      inputdf_sample()
    #}
    
  }, 
  extensions = c('Responsive','Scroller'),
  options = list(
    deferRender = TRUE,
    scrollY = 300,
    scroller = TRUE,
    searching = F
  ))
  
  # 3.14 绘图展示熔解曲线 ----
  output$plot_melt <- renderPlotly({
    if (rv$isMelt) {
      tryCatch({
        plot <- pcr_melt(
            df_melt = inputdf_melt(),
            df_sample = inputdf_sample(),
            sample = input$melt_sample,
            gene = input$melt_gene ) +
          ggtitle("Melting Curve")
        ggplotly(plot)
      }, error = function(e) {
        message("An error occurred: ", e$message)
        return(NULL)
      })
    } else {
      tryCatch({
        plot <- pcr_melt(
          df_melt = inputdf_melt(),
          df_sample = inputdf_sample(),
          sample = input$melt_sample,
          gene = input$melt_gene ) +
          ggtitle("This is Demo data") +
          theme(plot.title = element_text(colour = "red"))
        ggplotly(plot)
      }, error = function(e) {
        message("An error occurred: ", e$message)
        return(NULL)
      })
    }
    
  })
  
  #3.15 图片 ----
  output$img1 <- renderImage({
    list(src = "www/picture/p1.png",
         alt = "Description of the image",
         width = "300px",  
         height = "200px"  
    )
  })
  
  output$img_abs <- renderImage({
    list(src = "www/picture/p2.png",
         alt = "Description of the image",
         width = "300px",  
         height = "200px"  
    )
  })
  
  output$img_dCt <- renderImage({
    list(src = "www/picture/p3.png",
         alt = "Description of the image",
         width = "600px",  
         height = "320px"  
    )
  })
  
  output$img_ddCt <- renderImage({
    list(src = "www/picture/p4.png",
         alt = "Description of the image",
         width = "630px",  
         height = "280px"  
    )
  })
  
  output$img_sc <- renderImage({
    list(src = "www/picture/p5.png",
         alt = "Description of the image",
         width = "600px",  
         height = "300px"  
    )
  })
  
  output$img_home <- renderImage({
    list(src = "www/picture/p6.png",
         alt = "Description of the image",
         width = "300px",  
         height = "270px"  
    )
  })
  
 
  
  
  
  
  
  
  
}

#### 4. 启动shiny App ####
shinyApp(ui = ui, server = server)