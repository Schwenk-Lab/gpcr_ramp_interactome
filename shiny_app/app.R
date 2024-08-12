### --- GPCR-RAMP complex detection app --- ###
## Leo Dahl
# Some modules are reused from the GPCR Ab validation app due to their related nature


### Packages ###
library(shiny)
library(tidyverse)
library(DT)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(plotly)
library(bslib)
library(sass)

### Modules ###
source("modules.R")

### Load data ###
dat_list <- readRDS("gpcr_cxdet_data.rds")

ui <- fluidPage(
  title = "GPCR-RAMP interactome app",
  
  navbarPage(
    theme = bs_theme(version = 3, bootswatch = "flatly") %>%
      # Extra rules for tables and such
      bs_add_rules(sass_file("www/shiny_app.scss")),
    collapsible = T,
    title = NULL,
    
    # Landing page
    tabPanel(
      strong("GPCR-RAMP interactome Study"),
      HTML("Companion app to<h2><i>Multiplexed Mapping of the Interactome of GPCRs with Receptor Activity-Modifying Proteins</i></h2>"),
      HTML("Ilana B. Kotliar, Annika Bendes, Leo Dahl, Yuanhuang Chen, Marcus Saarinen, Emilie Ceraudo, Tea Dodig-Crnković, Mathias Uhlén, Per Svenningsson, Jochen M. Schwenk, Thomas P. Sakmar"), br(),
      HTML("<a href='https://doi.org/10.1126/sciadv.ado9959' target='_blank' rel='noopener noreferrer'><strong>https://doi.org/10.1126/sciadv.ado9959</strong></a>"),
      h3("About"),
      HTML("<p><strong>Receptor activity–modifying proteins (RAMPs)</strong> form complexes with G protein–coupled receptors (GPCRs) and may regulate their cellular trafficking and pharmacology. RAMP interactions have been identified for about 50 GPCRs, but only a few GPCR-RAMP complexes have been studied in detail. To elucidate a comprehensive GPCR-RAMP interactome, we created a library of 215 dual epitope-tagged (DuET) GPCRs representing all GPCR subfamilies and coexpressed each GPCR with each of the three RAMPs. Screening the GPCR-RAMP pairs with customized multiplexed suspension bead array (SBA) immunoassays, we identified 122 GPCRs that showed strong evidence for interaction with at least one RAMP. We screened for interactions in three cell lines and found 23 endogenously expressed GPCRs that formed complexes with RAMPs. Mapping the GPCR-RAMP interactome expands the current system-wide functional characterization of RAMP-interacting GPCRs to inform the design of selective therapeutics targeting GPCR-RAMP complexes.</p>"),
      HTML("<p><strong>This app</strong> provides complementary visualizations for the interactome analysis, such as per-GPCR summaries
           of the detected interactions and customisable overview heatmaps.</p>"),
      HTML("<p><strong>Data</strong> can be found at the
           <a href='https://doi.org/10.17044/scilifelab.24999260.v1' target='_blank' rel='noopener noreferrer'><strong>SciLifeLab Data Repository on FigShare</strong></a>
           and <strong> code </strong> is found at the 
           <a href='https://github.com/Schwenk-Lab/gpcr_ramp_interactome' target='_blank' rel='noopener noreferrer'><strong> Schwenk Lab Gitub account</strong></a> and the associated
           <a href='https://zenodo.org/doi/10.5281/zenodo.11371645' target='_blank' rel='noopener noreferrer'> <strong> Zenodo repository </strong></a>."),
      hr(),
      HTML("App developed by Leo Dahl"), br(),
      HTML("App version 1.0.2"), br(),
      HTML("<details><summary>Click for version history</summary>
           <strong>v1.0.2</strong>: <br>
           - Fix incorrect labels in cell outline legend for heatmap in GPCR-RAMP interactome tab. <br>
           <strong>v1.0.1</strong>: <br>
           - Update with DOI links after publication. <br>
           - Set GLP1R to be shown in the beginning. <br>
           - Make cell outlines thicker and update legend for heatmap in GPCR-RAMP interactome tab. <br>
           - Add version history.
           </details>"), br(),
      HTML("Using the bslib Flatly theme")
    ),
    
    # Complex detection, per GPCR
    tabPanel(
      "GPCR-RAMP interaction",
      
      sidebarLayout(
        sidebarPanel(
          gpcr_select_ui("per_gpcr_select"),
          HTML(
            "<hr class='sidebarhr'>
            <details><summary><strong>Click for description</strong></summary>
            <p>This page shows GPCR-RAMP interactions per GPCR. The GPCR to display can be selected from the drop-down
            menus to the left, either by GPCR name or by UniProt ID. The table to the right shows information about
            the selected GPCR together with links to external resources (Human Protein Atlas (HPA), Ensembl, UniProt).
            <p><strong>Detected GPCR-RAMP interactions</strong> are shown as a heatmap and a circle plot.
            <strong>The interaction heatmap</strong> shows interactions as pass/fail for each capture scheme. The cell
            outlines (if there are any) show if there is previous literature showing interaction (solid outline),
            no interaction (dotted outline), or if there are multiple sources with conflicting reports (dashed outline).
            <strong>The circle plot</strong> shows interactions as links between circle sectors. Links are colour coded
            per RAMP, with RAMP1 being light blue, RAMP2 being light green, and RAMP3 being orange. Each link represents
            one capture scheme and the direction of the link (shown as an arrow) indicates the capture, at the base of the
            arrow, and the detection, at the tip of the arrow.</p>
            <p><strong>The density plots</strong> show the density of the data points for each Ab targeting the GPCR directly
            (i.e. not the RAMP or an epitope tag, says 'No data' if no data is available). Data points representing the
            selected GPCR together with a RAMP are shown as lines, colour coded in the same way as the circle plot links. The
            Ab-specific threshold for interaction is shown as a dotted vertical line.</p>
            <p><strong>The endogenous interactions heatmap</strong> displays endogenous GPCR-RAMP interactions that were
            detected in three different cell lines, given that the GPCR was tested for (otherwise says 'No data').</p>
            </details>
            "
          )
        ),
        
        mainPanel(
          gpcr_info_ui("gpcr_info"),
          splitLayout(
            single_gpcr_hm_ui("cxdet_hm"),
            circle_plot_ui("cxdet_circle"),
            cellWidths = c("40%", "60%")
          ),
          dens_plot_ui("cxdet_dens"),
          endog_hm_ui("endog_hm")
        )
      )
    ),
    
    # Complex detection, summary
    tabPanel(
      "Overview heatmap",
      
      sidebarLayout(
        sidebarPanel(
          all_hm_custom_ui("all_hm_custom"),
          HTML("<hr class='sidebarhr'>
               <details><summary><strong>Click for description</strong></summary>
               <p><strong>The heatmap</strong> shows an overview of detected GPCR-RAMP interactions. By default,
               the heatmap shows a summary as the fraction of passing capture schemes for each GPCR from epitope
               capture and protein capture. For protein capture only anti-GPCR Abs previously validated (Dahl, 2023)
               against their target GPCR are shown, resulting in gaps for some GPCRs. The GPCRs without protein
               capture can be hidden if desired. 
               <p>It is possible to show the colour scale as a discrete scale based on the amount
               of evidence for interaction, where 'Strong' means that the fraction of passing capture schemes is
               higher than 2/3, 'Medium' means that the fraction is higher than 1/3 but lower than or equal to 2/3,
               and 'Weak' means the the fraction is lower than or equal to 1/3.</p>
               <p>Results only from epitope capture or protein capture can also be hidden entirely. Another way to
               display the heatmap is through raw interactions, which shows 'Pass' or 'Fail' for a given capture scheme.</p>
               <p>The heatmap can be customised further by selecting a subfamily/ligand family or changing the size
               of the heatmap, which is done through the controls that appear when clicking above where it says to 'Click here'.</p>
               </details>")
        ),
        
        mainPanel(
          all_hm_ui("all_hm")
        )
      )
    ),
    
    tabPanel(
      "Session information",
      session_info_ui("session_info")
    )
  )
)

server <- function(input, output, session) {
  
  ## Per GPCR interactions
  # Selected GPCR
  sel_gpcr <- gpcr_select_server("per_gpcr_select", dat_list)
  
  # Information about selected GPCR
  gpcr_info_server("gpcr_info", sel_gpcr, dat_list)
  
  # Circle plot of interactions
  circle_plot_server("cxdet_circle", sel_gpcr, dat_list$cxdet_epi, dat_list$cxdet_hpa)
  
  # Heatmap of interactions
  single_gpcr_hm_server("cxdet_hm", sel_gpcr, dat_list$cxdet_epi, dat_list$cxdet_hpa)
  
  # Density plots per anti-GPCR Ab
  dens_plot_server("cxdet_dens", sel_gpcr, dat_list)
  
  endog_hm_server("endog_hm", sel_gpcr, dat_list$endog_hits)
  
  ## Overview
  # Heatmap customisation
  hm_param <- all_hm_custom_server("all_hm_custom", dat_list)
  
  # Making the heatmap
  all_hm_server("all_hm", dat_list, hm_param)
  
  session_info_server("session_info")
}

shinyApp(ui, server)
