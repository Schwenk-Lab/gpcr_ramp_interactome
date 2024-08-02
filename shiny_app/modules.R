### --- GPCR-RAMP shiny app modules --- ###


### Per-GPCR complex detection information ###
## List of GPCRs for selection in analysis
gpcr_select_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("gpcr_select")),
    uiOutput(ns("uniprot_select"))
  )
}

gpcr_select_server <- function(id, data_in) {
  # data_in, dataset containing the long data structure
  moduleServer(id, function(input, output, session) {
    
    # GPCR names to choose from
    gpcrs <- data_in$long_dat_all %>%
      distinct(gpcr) %>%
      filter(!gpcr %in% c("empty", "mock")) %>%
      pull(gpcr) %>%
      sort() %>%
      setNames(NULL)
    # UniProt IDs in same order as GPCRs
    uniprots <- data_in$long_dat_all[match(gpcrs, data_in$long_dat_all$gpcr), "uniprot", drop = T]
    # Matching GPCR and uniprot for syncing
    gpcr_uniprot <- gpcrs %>% `names<-`(., uniprots)
    
    # Make GPCR input UI with GPCRs from data
    output$gpcr_select <- renderUI({
      selectInput(NS(id, "choose_gpcr"), "Select GPCR", gpcrs, selected = "GLP1R")
    })
    
    output$uniprot_select <- renderUI({
      selectInput(NS(id, "choose_uniprot"), "Select UniProt ID", uniprots, selected = "P43220")
    })
    
    # Sync the inputs
    observeEvent(input$choose_uniprot, {
      updateSelectInput(session, "choose_gpcr", "Select GPCR", gpcrs, gpcr_uniprot[input$choose_uniprot])
    })
    observeEvent(input$choose_gpcr, {
      updateSelectInput(session, "choose_uniprot", "Select UniProt ID", uniprots, names(gpcr_uniprot)[gpcr_uniprot == input$choose_gpcr])
    })
    
    return(reactive(input$choose_gpcr))
  })
}

## Table with GPCR information
gpcr_info_ui <- function(id) {
  ns <- NS(id)
  uiOutput(ns("gpcr_tbl"))
}

gpcr_info_server <- function(id, gpcr, data_in) {
  # gpcr, reactive string of the selected GPCR
  # data_in, data set containing the long data
  moduleServer(id, function(input, output, session) {
    
    # Get one row of info for current GPCR
    current_gpcr <- reactive({
      data_in$long_dat_all %>%
        filter(gpcr == req(gpcr())) %>%
        select(gpcr, ensg, uniprot, class, subfamily, family, contains("int_ref")) %>%
        distinct() %>%
        # Sort across everything so that NA ends up at the bottom, in case some row happens to have a missing value where it should not
        arrange(across(everything())) %>%
        slice(1)
    })
    
    # Get expected/not expected interactions from literature
    int_yes <- reactive({
      ramps <- current_gpcr() %>% select(gpcr, starts_with("int_ref")) %>%
        pivot_longer(cols = -gpcr) %>% filter(!is.na(value)) %>%
        pull(name) %>% str_extract_all(., "ramp\\d", simplify = F) %>% toupper() %>% paste0(., collapse = ", ")
      
      if (ramps == "") {ramps <- "None"}
      return(ramps)
    })
    int_no <- reactive({
      ramps <- current_gpcr() %>% select(gpcr, starts_with("no_int")) %>%
        pivot_longer(cols = -gpcr) %>% filter(!is.na(value)) %>%
        pull(name) %>% str_extract_all(., "ramp\\d", simplify = F) %>% toupper() %>% paste0(., collapse = ", ")
      
      if (ramps == "") {ramps <- "None"}
      return(ramps)
    })
    
    # GPCR table containing links to HPA, Ensembl, UniProt, different identifiers, info about family/subfamily
    output$gpcr_tbl <- renderUI({
      HTML(paste0(
        "<h1>", req(gpcr()), "</h1>",
        "<p>
        <a href='", paste0("https://www.proteinatlas.org/",
                           current_gpcr()$ensg, "-", req(gpcr())),
        "' class='gpcrlink' target='_blank' rel='noopener noreferrer'>HPA portal</a>
        <a href='", paste0("https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",
                           current_gpcr()$ensg),
        "' class='gpcrlink' target='_blank' rel='noopener noreferrer'>Ensembl</a>
        <a href='", paste0("https://www.uniprot.org/uniprot/",
                           str_remove(current_gpcr()$uniprot, "[[:space:]]")), # Remove space at the end of uniprot IDs (if present), breaks link
        "' class='gpcrlink' target='_blank' rel='noopener noreferrer'>UniProt</a>
        </p>",
        "<table style='width:100%' class='gpcrinfo'>
          <tr>
            <td>Ensembl gene ID</td>
            <td>", current_gpcr()$ensg,"</td>
          </tr>
          <tr>
            <td>UniProt ID</td>
            <td>", current_gpcr()$uniprot, "</td>
          </tr>
          <tr>
            <td>Class</td>
            <td>", current_gpcr()$class, "</td>
          </tr>
          <tr>
            <td>Subfamily</td>
            <td>", current_gpcr()$subfamily, "</td>
          </tr>
          <tr>
            <td>Ligand Family</td>
            <td>", current_gpcr()$family, "</td>
          </tr>
          <tr>
            <td>Interactors from literature</td>
            <td>", int_yes(), "</td>
          </tr>
          <tr>
            <td>Non-interactors from literature</td>
            <td>", int_no(), "</td>
          </tr>
        </table>"
      ))
    })
  })
}

## Circle plot of detected GPCR-RAMP interactions
circle_plot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Interaction circle plot"),
    plotOutput(ns("gpcr_ramp_circle"), height = "500px", width = "500px")
  )
}

circle_plot_server <- function(id, gpcr_in, cxdet_epi, cxdet_hpa) {
  # gpcr_in, reactive input of selected GPCR
  # cxdet_epi, complex detection results for epitope capture
  # cxdet_hpa, complex detection results for HPA (protein) capture
  
  moduleServer(id, function(input, output, session) {
    
    output$gpcr_ramp_circle <- renderPlot({
      # Check that a GPCR has been selected
      g <- req(gpcr_in())
      
      # Find Abs binding the selected GPCR
      gab <- cxdet_hpa$interactions %>%
        # Add a dot to not get other GPCRs whose names contain the selected GPCR as a substring
        filter(str_detect(gene_name_id, paste0(g, "\\."))) %>%
        # Get antibody name
        separate(gene_name_ab, into = c("gpcr", "ab"), sep = "_") %>%
        pull(ab) %>%
        unique()
      
      # Make data frame to use for plotting
      a <- rbind(
        cxdet_epi$interactions %>% filter(which_gpcr == g) %>% select(-c(signal, cutoff, gpcr_expr, int_exp)),
        cxdet_hpa$interactions %>% separate(gene_name_ab, into = c("gpcr", "capt"), sep = "_") %>%
          filter(capt %in% gab) %>% select(-c(gene_name_id, gpcr, cutoff)) %>%
          pivot_longer(cols = -capt, names_to = "which_ramp", values_to = "pass") %>%
          mutate(which_gpcr = g, det = "OLLAS", pass = as.numeric(pass)) %>%
          select(which_gpcr, which_ramp, capt, det, pass)
      ) %>%
        # Add unique IDs for each capture-detection scheme
        mutate(capt_name = case_when(capt %in% c("1D4", "FLAG", gab) ~ paste0(which_gpcr, "_", capt), T ~ paste0(which_ramp, "_", capt)),
               det_name = case_when(det == "OLLAS" ~ paste0(which_ramp, "_", det), det == "1D4" ~ paste0(which_gpcr, "_", det))) %>%
        # rearrange order of links a bit to reduce link crossings
        arrange(desc(which_ramp), desc(capt), desc(det))
      
      # Sectors to use in plot
      s <- c(paste0(g, "_", c("1D4", "FLAG", gab)),
             paste0(rep("RAMP1", 3), c("_HA", "_OLLAS", "_RAMP1")),
             paste0(rep("RAMP2", 3), c("_HA", "_OLLAS", "_RAMP2")),
             paste0(rep("RAMP3", 3), c("_HA", "_OLLAS", "_RAMP3")))
      s <- factor(s, levels = s)
      
      # Colours depending on associated molecule
      circ_col <- strsplit(as.character(s), "_") %>%
        sapply(., function(x) {x[1]})
      circ_col <- case_when(circ_col == g ~ "#F2F2F2",
                            circ_col == "RAMP1" ~ "#BFEFFF",
                            circ_col == "RAMP2" ~ "#B1E063",
                            circ_col == "RAMP3" ~ "#EEAD0E")
      
      # Links in plot, should be directional
      link_df <- a %>%
        filter(pass == 1) %>%
        select(which_ramp, capt_name, det_name) %>%
        mutate(colour = case_when(which_ramp == "RAMP1" ~ "#BFEFFF",
                                  which_ramp == "RAMP2" ~ "#B1E063",
                                  which_ramp == "RAMP3" ~ "#EEAD0E"))
      
      if (nrow(link_df) > 0) {
        link_counts <- table(c(link_df$capt_name, link_df$det_name))
        link_mids <- sapply(names(link_counts), function(x) {
          # Even number: start on integers like 1, 2, 3
          # Odd number: start on half numbers like 1.5, 2.5, 3.5
          # Get the middle points of each link
          n <- link_counts[x]
          link_mid <- seq(from = 6 - n/2 + 0.5, by = 1, length.out = n)
          
          data.frame(link = x, link_mid = link_mid)
        }, simplify = F) %>% rbindlist() %>% setDF() %>% mutate(link = make.names(link, unique = T))
        # Reorder to match links in link_df
        link_mids <- link_mids[match(make.names(c(link_df$capt_name, link_df$det_name), unique = T), link_mids$link), ]
        link_df$link_mid_capt <- link_mids[1:nrow(link_df), "link_mid"]
        link_df$link_mid_det <- link_mids[(nrow(link_df)+1):nrow(link_mids), "link_mid"]
      }
      
      # Get the middle of the GPCR section to be at 12 o'clock
      # 0 degrees is 3 o'clock -> 90 for 12 (start of first GPCR sector will be at 12, not the middle)
      # Then add the length of half of the number of sectors that are from the GPCR
      stdeg <- 90 + (360 / length(s)) * (sum(grepl(g, s)) / 2)
      
      circos.clear()
      
      # Max number of possible links is 12 for 1D4, set xlim between 0 and 12
      circos.par("track.height" = 0.1, "start.degree" = stdeg, "gap.degree" = 1, "track.margin" = c(0,0))
      par(mar = c(1, 1, 2, 1))
      circos.initialize(sectors = s, xlim = c(0, 12))
      
      title(main = g, cex.main = 2.5, line = 0)
      
      circos.track(sectors = s, ylim = c(0, 1), bg.col = circ_col, bg.border = "grey40", panel.fun = function(x, y) {
        # Add text for molecule associated with sector
        txt <- (CELL_META$sector.index %>%
                  strsplit(., "_") %>%
                  unlist())[2] %>%
          str_remove(paste0(g, "\\."))
        
        circos.text(CELL_META$xcenter, CELL_META$ycenter, txt, facing = "bending.inside", cex = 0.9, font = 2)
      })
      
      if (nrow(link_df) > 0) {
        # Add links
        for (i in 1:nrow(link_df)) {
          circos.link(link_df$capt_name[i],
                      c(link_df$link_mid_capt[i] - 0.5, link_df$link_mid_capt[i] + 0.5),
                      link_df$det_name[i],
                      c(link_df$link_mid_det[i] - 0.5, link_df$link_mid_det[i] + 0.5),
                      col = link_df$colour[i], arr.type = "big.arrow", arr.length = 0.02,
                      border = "grey50", lwd = 0.5, directional = 1)
        }
      }
      
      # Highlight sectors with text showing if it is a GPCR or which RAMP
      for (i in c(paste0("RAMP", 1:3))) {
        cells <- grep(paste0("^", i), s)
        # Sectors to highlight
        # If even number, put between middle two
        # If odd, put on middle one
        if (length(cells) %% 2 == 0) {
          h_sectors <- s[cells[c(length(cells)/2, length(cells)/2 + 1)]] %>% as.character()
        } else {
          h_sectors <- s[cells[ceiling(length(cells)/2)]] %>% as.character()
        }
        
        highlight.sector(h_sectors, 1, text = i, col = NULL, text.vjust = "7mm", facing = "bending.inside", cex = 1.5, font = 2)
      }
    })
    
  })
}

## Plot showing data used for protein capture complex detection
dens_plot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Ab density plots"),
    plotOutput(ns("dens_plot"), height = "600px")
  )
}

dens_plot_server <- function(id, gpcr_in, data_in) {
  # gpcr_in, reactive input of selected GPCR
  # data_in, list of data frames with data
  
  moduleServer(id, function(input, output, session) {
    # Density plot
    output$dens_plot <- renderPlot({
      g <- req(gpcr_in())
      
      # Check for validated Abs
      v_abs <- data_in$prot %>% filter(str_detect(gene_name_id, paste0(g, "\\."))) %>%
        filter(pass) %>% pull(gene_name_id)
      
      if (length(v_abs) > 0) {
        plt_dat <- data_in$long_dat_all %>% filter(gene_name_id %in% v_abs) %>%
          # Look at detection by binding OLLAS
          filter(is.na(exclude), detection_ab == "OLLAS") %>%
          unite("gene_name_ab", gene_name, antibody, sep = "_") %>%
          select(gpcr, gene_name_ab, interactor, value_rz, target)
        
        # Get data for points to mark (target GPCR + RAMP1/2/3)
        lines_df <- filter(plt_dat, target == "Target" & interactor != "mock")
        # Get Ab-specific threshold
        thr_df <- lines_df %>% group_by(gene_name_ab) %>% summarise(gene_name_ab = unique(gene_name_ab)) %>%
          mutate(threshold = data_in$cxdet_hpa$interactions[match(gene_name_ab, data_in$cxdet_hpa$interactions$gene_name_ab), "cutoff", drop = T],
                 colr = "Threshold")
        
        # Same text size everywhere
        txt_sz <- 16
        
        plt <- ggplot(plt_dat) +
          geom_density(aes(x = value_rz, y = after_stat(scaled)), size = 1) +
          geom_vline(data = lines_df, aes(xintercept = value_rz, colour = interactor, linetype = interactor), size = 1) +
          geom_vline(data = thr_df, aes(xintercept = threshold, colour = colr, linetype = colr), size = 1) +
          facet_wrap(~ gene_name_ab, scales = "free_x") +
          scale_colour_manual(values = c("RAMP1" = "#BFEFFF", "RAMP2" = "#B1E063", "RAMP3" = "#EEAD0E", "Threshold" = "black")) +
          scale_linetype_manual(values = c("RAMP1" = "solid", "RAMP2" = "solid", "RAMP3" = "solid", "Threshold" = "dotted")) +
          labs(x = "Signal", y = "Density", colour = NULL, linetype = NULL) +
          theme_classic() +
          theme(strip.text = element_text(size = txt_sz),
                legend.text = element_text(size = txt_sz),
                axis.text.x = element_text(size = txt_sz),
                axis.text.y = element_text(size = txt_sz),
                axis.title.x = element_text(size = txt_sz),
                axis.title.y = element_text(size = txt_sz))
      } else {
        plt <- ggplot() + labs(title = "No data")
      }
      
      return(plt)
    })
  })
}

## Heatmap of interactions for the selected GPCR
single_gpcr_hm_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Interaction heatmap"),
    plotOutput(ns("hm"), height = "500px", width = "300px")
  )
}

single_gpcr_hm_server <- function(id, gpcr_in, cxdet_epi, cxdet_hpa) {
  # gpcr_in, reactive input of selected GPCR
  # cxdet_epi, complex detection results for epitope capture
  # cxdet_hpa, complex detection results for HPA (protein) capture
  
  moduleServer(id, function(input, output, session) {
    output$hm <- renderPlot({
      g <- req(gpcr_in())
      
      # Get interactions from both capture schemes, reshape to fit together
      plt_dat <- rbind(
        cxdet_epi$interactions %>% filter(which_gpcr == g) %>% select(which_gpcr, which_ramp, capt, pass) %>%
          mutate(capt = case_when(str_detect(capt, "RAMP") ~ "RAMP", T ~ capt)) %>% pivot_wider(id_cols = capt, names_from = which_ramp, values_from = pass),
        cxdet_hpa$interactions %>% filter(str_detect(gene_name_id, paste0(g, "\\."))) %>% select(gene_name_ab, contains("RAMP")) %>%
          mutate(across(where(is.logical), as.numeric)) %>% rename("capt" = "gene_name_ab")
      ) %>% column_to_rownames("capt") %>% mutate(across(everything(), function(x) {ifelse(x, "Pass", "Fail")}))
      
      # Get interaction expectations from literature
      int_exp <- cxdet_epi$interactions %>% filter(which_gpcr == g) %>%
        select(which_gpcr, which_ramp, int_exp) %>% distinct() %>% arrange(which_ramp) %>%
        # Convert literature interactions to line types
        mutate(int_exp = case_when(int_exp == "Expected" ~ "solid",
                                   int_exp == "None" ~ "dotted",
                                   int_exp == "Uncertain" ~ "dashed"))
      
      # Colours for heatmap
      col_vec <- RColorBrewer::brewer.pal(3, "Set2")[-3]
      names(col_vec) <- c("Pass", "Fail")
      
      hm <- Heatmap(as.matrix(plt_dat), col = col_vec, column_gap = unit(3, "points"), row_gap = unit(1, "points"),
                    cluster_rows = F, cluster_columns = F, row_split = 1:nrow(plt_dat), column_split = colnames(plt_dat), # Row splits seems to order the splits in alphabetical order, supply a numerical vector instead of names
                    show_row_names = T, row_title = NULL, show_column_names = F, row_names_side = "left", show_heatmap_legend = F,
                    row_names_gp = gpar(fontsize = 16), column_names_gp = gpar(fontsize = 14), column_title_gp = gpar(fontsize = 20), column_title_rot = 90,
                    
                    # Draw cell border based on literature interactions
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      if (!is.na(int_exp[j, "int_exp", drop = T])) {
                        border_lty <- int_exp[j, "int_exp", drop = T]
                        grid.rect(x = x, y = y, width = width, height = height,
                                  gp = gpar(lty = border_lty, lwd = 2, fill = NA, col = "black"))
                      }
                    })
      # Legends for colours and cell outlines
      hm_legend <- list(
        Legend(labels = c("Pass", "Fail"), title = "Interaction", type = "grid", legend_gp = gpar(fill = brewer.pal(3, "Set2")[-3]),
               labels_gp = gpar(fontsize = 13), title_gp = gpar(fontsize = 16),
               direction = "horizontal"),
        Legend(labels = c("Expected", "No interaction", "Uncertain", "No literature"),
               title = "Interaction in\nliterature", type = "lines",
               legend_gp = gpar(col = "black", lty = c(1, 2, 3, 0)),
               labels_gp = gpar(fontsize = 13), title_gp = gpar(fontsize = 16),
               direction = "horizontal")
      )
      
      return(draw(hm, annotation_legend_list = hm_legend, annotation_legend_side = "bottom"))
    })
    
  })
}

## Small heatmap of hits from endogenous analysis for selected GPCR
endog_hm_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Endogenous interactions"),
    plotOutput(ns("hm"), height = "300px", width = "400")
  )
}

endog_hm_server <- function(id, gpcr_in, endog_hits) {
  # gpcr_in, reactive input of selected GPCR
  # endog_hits, hits from endogenous analysis
  
  moduleServer(id, function(input, output, session) {
    
    output$hm <- renderPlot({
      g <- req(gpcr_in())
      
      hm_dat <- endog_hits %>%
        filter(gene_name == g)
      
      # Nothing to show if the GPCR is not present in the data (no validated Abs or no endogenous hits)
      if (nrow(hm_dat) == 0) {return(ggplot() + labs(title = "No data"))}
      
      # Colours for heatmap
      col_vec <- RColorBrewer::brewer.pal(3, "Set2")[-3]
      names(col_vec) <- c("Yes", "No")
      
      hm_dat <- hm_dat %>%
        pivot_wider(id_cols = cell_type, names_from = detection_ab, values_from = hit) %>%
        column_to_rownames("cell_type") %>%
        as.matrix()
      
      hm_out <- Heatmap(hm_dat,
                        row_names_side = "left", cluster_rows = F,
                        name = "Interaction\ndetected",
                        column_names_side = "top", cluster_columns = F, show_heatmap_legend = F,
                        col = col_vec, column_gap = unit(3, "points"), row_gap = unit(1, "points"),
                        row_split = 1:nrow(hm_dat), column_split = colnames(hm_dat),
                        show_row_names = T, show_column_names = T, row_title = NULL, column_title = NULL,
                        row_names_gp = gpar(fontsize = 16), column_names_gp = gpar(fontsize = 18))
      lgd <- list(
        Legend(labels = c("Yes", "No"), title = "Interaction\ndetected", type = "grid", legend_gp = gpar(fill = brewer.pal(3, "Set2")[-3]),
               labels_gp = gpar(fontsize = 16), title_gp = gpar(fontsize = 16))
      )
      
      return(draw(hm_out, annotation_legend_list = lgd))
    })
  })
}


### Overview of interactions ###
## Heatmap of interactions, customisation from sidebar
all_hm_custom_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # Customising content
    # Show raw interactions or fraction of passing capture schemes (summary)
    radioButtons(ns("int_or_sum"), "Display", c("Summary", "All interactions")),
    # Cannot combine heatmaps if showing interactions per capture scheme (as the number of Abs differs between GPCRs). Ok for showing fractions of passing schemes
    conditionalPanel(paste0("input['", ns("int_or_sum"), "'] == 'Summary'"),
                     radioButtons(ns("capt_sum"), "Capture", c("Epitope", "Protein", "Both"), "Both"),
                     # Choice between fraction pass and evidence strength for summary heatmaps
                     radioButtons(ns("frac_or_str"), "Summarisation", c("Fraction pass", "Interaction evidence"), "Fraction pass")),
    conditionalPanel(paste0("input['", ns("int_or_sum"), "'] == 'All interactions'"),
                     radioButtons(ns("capt_int"), "Capture", c("Epitope", "Protein"))),
    
    # Option to hide GPCRs without any Abs (if showing both capture methods)
    conditionalPanel(paste0("input['", ns("int_or_sum"), "'] == 'Summary' &
                             input['", ns("capt_sum"), "'] == 'Both'"),
                     checkboxInput(ns("hide_no_ab"), "Hide GPCRs without Abs", F)),
    # More content customisation
    HTML("<details><summary><strong>Click for more content customisation</strong></summary>"),
    # Subfamily, ligand family
    radioButtons(ns("filter_by"), "Filter by", c("Subfamily", "Ligand family"), "Subfamily"),
    uiOutput(ns("fam_sel")),
    HTML("</details>"),
    
    # Customising the heatmap layout
    HTML("<details><summary><strong>Click for layout customisation</strong></summary>"),
    h5("Heatmap size (px)"),
    splitLayout(
      numericInput(ns("hm_height"), "Height", value = 800, step = 10),
      numericInput(ns("hm_width"), "Width", value = 600, step = 10)
    ),
    HTML("</details>")
  )
}

all_hm_custom_server <- function(id, data_in) {
  # data_in, data set containing the long data
  moduleServer(id, function(input, output, session) {
    # Selection of subfamily/ligand family. In a renderui call as it uses values from the data
    output$fam_sel <- renderUI({
      fam_choices <- switch(input$filter_by,
                            "Subfamily" = c("All", "Rhodopsin (\u03B1)", "Rhodopsin (\u03B2)", "Rhodopsin (\u03B3)", "Rhodopsin (\u03B4)", "GSAF", "Other"),
                            "Ligand family" = c("All", data_in$long_dat_all %>% pull(family) %>% str_to_title() %>% unique() %>% sort() %>% na.omit()))
      
      selectInput(NS(id, "fam_select"), ifelse(input$filter_by == "Subfamily", "Select subfamily", "Select ligand family"), fam_choices)
    })
    
    # Different customisation options for the heatmap
    return(reactive(list(
      "int_or_sum" = input$int_or_sum,
      "capture_int" = input$capt_int,
      "capture_sum" = input$capt_sum,
      "frac_or_str" = input$frac_or_str,
      "hide_no_ab" = input$hide_no_ab,
      "hm_height" = input$hm_height,
      "hm_width" = input$hm_width,
      "filter_by" = input$filter_by,
      "fam_select" = input$fam_select
    )))
    
  })
}

## Heatmap of interactions, plotting
all_hm_ui <- function(id) {
  ns <- NS(id)
  tagList(
    plotlyOutput(ns("hm"))
  )
}

all_hm_server <- function(id, data_in, hm_param) {
  # data_in, data set containing the long data and interactions
  # hm_param, reactive list of parameters for heatmap
  moduleServer(id, function(input, output, session) {
    plt_dat <- reactive({
      if (hm_param()$int_or_sum == "Summary") {
        # Option to show epitope capture, protein capture, or both
        dat <- switch(
          hm_param()$capture_sum,
          "Epitope" = data_in$cxdet_epi$interaction_evidence %>%
            # Convert categorical axis to numeric for it to work
            mutate(int_evid = case_when(int_evid == "Weak" ~ 0, int_evid == "Medium" ~ 1, int_evid == "Strong" ~ 2)) %>%
            # Show fraction or discrete intervals as summary
            rename(zcol = ifelse(hm_param()$frac_or_str == "Fraction pass", "fraction_pass", "int_evid")),
          
          "Protein" = data_in$cxdet_hpa$interaction_evidence %>%
            mutate(int_evid = case_when(int_evid == "Weak" ~ 0, int_evid == "Medium" ~ 1, int_evid == "Strong" ~ 2)) %>%
            rename(zcol = ifelse(hm_param()$frac_or_str == "Fraction pass", "fraction_pass", "int_evid")),
          
          "Both" = rbind(
            data_in$cxdet_epi$interaction_evidence %>% mutate(capt = "Epitope"),
            data_in$cxdet_hpa$interaction_evidence %>% mutate(capt = "Protein")
          ) %>% mutate(int_evid = case_when(int_evid == "Weak" ~ 0, int_evid == "Medium" ~ 1, int_evid == "Strong" ~ 2)) %>%
            rename(zcol = ifelse(hm_param()$frac_or_str == "Fraction pass", "fraction_pass", "int_evid"))
        )
        
        # If both, allow to hide GPCRs without Abs, to not have lots of NA rows in the protein capture part
        if (hm_param()$capture_sum == "Both" & hm_param()$hide_no_ab) {
          dat <- dat %>% group_by(which_gpcr, which_ramp) %>% filter(if(!any("Protein" %in% capt)) {F} else {T}) %>% ungroup()
        }
        
      } else {
        # Option to show epitope capture or protein capture
        dat <- switch(
          hm_param()$capture_int,
          "Epitope" = data_in$cxdet_epi$interactions %>% select(which_gpcr, which_ramp, capt, pass) %>%
            # Remove number from RAMP for RAMP1/2/3 capture to not get empty columns
            mutate(capt = case_when(str_detect(capt, "RAMP") ~ "RAMP", T ~ capt)) %>%
            rename(zcol = pass),
          
          "Protein" = data_in$cxdet_hpa$interactions %>% select(gene_name_ab, contains("RAMP")) %>%
            pivot_longer(cols = -gene_name_ab, names_to = "which_ramp", values_to = "pass") %>%
            rename(which_gpcr = gene_name_ab, zcol = pass) %>% mutate(zcol = as.numeric(zcol))
        )
      }
      
      # Filter by subfamily or ligand family if desired
      # Error if anything but 'All' is selected and a switch between subfamily and ligand family occurs. This is due to the fam_select
      # value being used for both selections, but the content for fam_select changing with the filter_by selection. As the switch occurs,
      # the old fam_select value is incorrectly used briefly, resulting in an error as are no overlapping entries apart from 'All'. Use workaround
      # from the Ab validation app
      if ((
        req(hm_param()$filter_by) == "Subfamily" &
        req(hm_param()$fam_select) %in% c("Rhodopsin (\u03B1)", "Rhodopsin (\u03B2)", "Rhodopsin (\u03B3)", "Rhodopsin (\u03B4)", "GSAF", "Other")
      ) | (
        req(hm_param()$filter_by) == "Ligand family" &
        !req(hm_param()$fam_select) %in% c("All", "Rhodopsin (\u03B1)", "Rhodopsin (\u03B2)", "Rhodopsin (\u03B3)", "Rhodopsin (\u03B4)", "GSAF", "Other")
      )) {
        gpcrnames <- data_in$long_dat_all %>% group_by(gpcr) %>%
          summarise(f = ifelse(hm_param()$filter_by == "Subfamily", subfamily, family) %>% unique()) %>%
          filter(f == hm_param()$fam_select) %>% pull(gpcr)
        
        if (hm_param()$int_or_sum == "All interactions" & hm_param()$capture_in == "Protein") {
          # Deal with case when y-axis contains antibody names
          # Names contain "_" as separators between target and antibody names, and may contain ";" as separators between multiple targets
          # Make regex looking for GPCRname_ or GPCRname; to deal with multi-target Abs
          dat <- filter(dat, str_detect(which_gpcr, paste0(paste0(gpcrnames, "_", collapse = "|"), "|", paste0(gpcrnames, ";", collapse = "|"))))
          
        } else {
          # If y-axis contains GPCR names
          dat <- filter(dat, which_gpcr %in% gpcrnames)
        }
      }
      
      # Set y-axis to be alphabetical, ascending from the top (opposite of default)
      dat$which_gpcr <- factor(dat$which_gpcr, levels = rev(sort(unique(dat$which_gpcr))))
      
      return(dat)
    })
    
    output$hm <- renderPlotly({
      
      # Colour schemes for plots, depending on two parameters
      if (hm_param()$int_or_sum == "Summary" & hm_param()$frac_or_str == "Fraction pass") {
        hm_col <- viridis(101)
        cscale <- data.frame(z = seq(0, 1, length.out = 101), col = hm_col)
        cbar <- list(tickvals = seq(0, 1, by = 0.25), ticktext = seq(0, 1, by = 0.25))
      } else if (hm_param()$int_or_sum == "Summary" & hm_param()$frac_or_str == "Interaction evidence") {
        hm_col <- rep(brewer.pal(3, "Set2")[3:1], each = 2)
        cscale <- data.frame(z = c(0, 1/3, 1/3, 2/3, 2/3, 1), col = hm_col)
        cbar <- list(tickvals = c(1/3, 1, 5/3), ticktext = c("Weak", "Medium", "Strong"))
      } else {
        hm_col <- rep(brewer.pal(3, "Set2")[2:1], each = 2)
        cscale <- data.frame(z = c(0, 0.5, 0.5, 1), col = hm_col)
        cbar <- list(tickvals = c(0.25, 0.75), ticktext = c("Fail", "Pass"))
      }
      
      if ((hm_param()$int_or_sum == "Summary" & hm_param()$capture_sum == "Both") |
          (hm_param()$int_or_sum == "All interactions" & hm_param()$capture_int == "Epitope")) {
        # Need to plot subplot
        plt <- plt_dat() %>%
          nest_by(which_ramp) %>% mutate(p = list(
            # Lambda function to make multi-line code work in the list() call
            (function() {
              plt <- plot_ly(data = data, x = ~ capt, y = ~ which_gpcr, z = ~ zcol, name = which_ramp, height = hm_param()$hm_height, width = hm_param()$hm_width, legendgroup = "leg") %>%
                # Manually add titles for the subplots
                layout(
                  annotations = list(
                    text = which_ramp, xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom", align = "center",
                    x = 0.5, y = 1, showarrow = F, font = list(size = 18, color = "black")),
                  yaxis = list(title = list(text = "GPCR"))
                  )
              
              # Use different colours and colour bars depending on selections
              if (hm_param()$int_or_sum == "Summary" & hm_param()$frac_or_str == "Fraction pass") {
                plt <- plt %>% plotly::add_heatmap(colorscale = cscale, colorbar = cbar, xgap = 1, ygap = 0.3, zmin = 0, zmax = 1)
              } else {
                plt <- plt %>% plotly::add_heatmap(colorscale = cscale, colorbar = cbar, xgap = 1, ygap = 0.3)
              }
              
              if (which_ramp != "RAMP1") {plt <- plt %>% hide_colorbar()}
              
              return(plt)
            })()
          )) %>% pull(p) %>%
          subplot(nrows = 1, shareY = T) %>%
          # Add x axis title to middle plot, seems to not work if putting it in the layout call in the loop and having titleX = F in subplot
          layout(xaxis2 = list(title = list(text = "Capture")))
        
      } else {
        # If no need for facets
        plt <- plot_ly(data = plt_dat(), x = ~ which_ramp, y = ~ which_gpcr, z = ~ zcol, height = hm_param()$hm_height, width = hm_param()$hm_width, legendgroup = "leg") %>%
          layout(xaxis = list(title = list(text = ""), side = "top", tickfont = list(size = 18, colour = "black"), ticks = ""),
                 yaxis = list(title = list(text = ifelse(hm_param()$capture_int == "Protein", "GPCR_Antibody", "GPCR"))))
        
        if ((hm_param()$int_or_sum == "All interactions" & hm_param()$capture_int == "Protein") |
            (hm_param()$int_or_sum == "Summary" & hm_param()$frac_or_str == "Interaction evidence")) {
          # Need discretised colour scale for categorical heatmaps
          plt <- plt %>% plotly::add_heatmap(colorscale = cscale, colorbar = cbar, xgap = 20, ygap = 0.3)
        } else {
          plt <- plt %>% plotly::add_heatmap(colorscale = cscale, colorbar = cbar, xgap = 20, ygap = 0.3, zmin = 0, zmax = 1)
        }
      }
      
      # Add a title to the colour bar, name depending on selections
      plt <- plt %>% colorbar(title = case_when(hm_param()$int_or_sum == "All interactions" ~ "Interaction",
                                                hm_param()$frac_or_str == "Fraction pass" ~ "Fraction",
                                                hm_param()$frac_or_str == "Interaction evidence" ~ "Evidence of\ninteraction"))
      
      return(plt)
    })
  })
}

### Session information ###
session_info_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("System information"),
    htmlOutput(ns("system_info")),
    h4("Added packages"),
    tableOutput(ns("other_pckg_versions"))
  )
}

session_info_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    s <- sessionInfo()
    
    output$system_info <- renderUI({
      paste0(
        s$R.version$version.string, "<br>",
        "Platform: ", s$platform, "<br>",
        "Running under: ", s$running
      ) %>% HTML()
    })
    
    # Function for making a table of packages, their versions, and the citations from the citation() function
    pckg_version_tbl <- function(pckg_names) {
      pckg_versions <- sapply(pckg_names, getNamespaceVersion)
      
      # Also get citations for each package
      pckg_citations <- sapply(pckg_names, function(x) {
        # Get citation in HTML format, edit links to open in new tabs
        citation(x) %>% format(style = "html") %>% str_replace_all("<a href", "<a target='_blank' rel='noopener noreferrer' href") %>%
          unlist() %>% paste(collapse = "")
      })
      
      tbl_out <- data.frame("Package" = pckg_names,
                            "Version" = pckg_versions,
                            "Citation" = pckg_citations)
      return(tbl_out)
    }
    
    # Other attached packages (not base)
    output$other_pckg_versions <- renderTable({
      pckg_version_tbl(sort(names(s$otherPkgs)))
    }, striped = T, sanitize.text.function = function(x) x)
    
  })
}


