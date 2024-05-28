### Miscellaneous functions for GPCR RAMP complex detection analysis

# Wrapper of the symnum function to make all tests use the same pvalue-to-star conversion rate
p_to_star <- function(p) {
  star_out <- symnum(p, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                     symbols = c("****", "***", "**", "*", ""), corr = F, legend = F) %>%
    as.character()
  
  return(star_out)
}

# Quartile plots for looking at complex detection vs RAMP expression
quartile_plots <- function(long_dat_all, cxdet_epitope, cxdet_hpa, which_expr = c("RAMP", "GPCR")) {
  # long_dat_all, long format data with anti-HA antibodies and OLLAS detection
  # cxdet_epitope, data frame with interaction evidence strength from epitope capture complex detection
  # cxdet_hpa, data frame with interaction evidence strength from HPA capture complex detection
  # which_expr, string indicating what to show expression of, RAMP or GPCR (takes first element)
  
  plts <- long_dat_all %>%
    filter(if (which_expr[1] == "RAMP") {gene_name == "M anti-HA"} else if (which_expr[1] == "GPCR") {gene_name == "M anti-FLAG"},
           if (which_expr[1] == "RAMP") {detection_ab == "OLLAS"} else if (which_expr[1] == "GPCR") {detection_ab == "1D4"},
           sample_type != "gpcr.p11", !interactor %in% c("mock", "empty")) %>%
    rename("which_gpcr" = "gpcr", "which_ramp" = "interactor") %>%
    select(unique_sample_name, which_gpcr, which_ramp, value) %>%
    # Get mean signal for each GPCR-RAMP combination, in case there are multiple samples
    group_by(which_gpcr, which_ramp) %>%
    summarise(mean_mfi = mean(value, na.rm = T), .groups = "keep") %>%
    ungroup() %>%
    # Add interaction evidence data
    left_join(cxdet_epitope %>% select(contains("which"), int_evid) %>% rename("epitope" = "int_evid"),
              by = c("which_gpcr", "which_ramp")) %>%
    left_join(cxdet_hpa %>% select(contains("which"), int_evid) %>% rename("hpa" = "int_evid"),
              by = c("which_gpcr", "which_ramp")) %>%
    # Make plots per capture method
    pivot_longer(cols = -c(which_gpcr, which_ramp, mean_mfi), names_to = "capture", values_to = "int_evid") %>%
    group_by(capture) %>%
    nest() %>%
    mutate(plt = map(data, ~ {
      # Colours slightly different between epitope and HPA capture
      bar_col <- switch(capture,
                        "epitope" = brewer.pal(3, "Set2") %>% setNames(c("Strong", "Medium", "Weak")),
                        "hpa" = brewer.pal(3, "Dark2") %>% setNames(c("Strong", "Medium", "Weak")))
      
      plt_out <- .x %>%
        filter(!is.na(int_evid)) %>%
        # Reorder evidence levels to make more sense
        mutate(int_evid = factor(int_evid, levels = c("Strong", "Medium", "Weak"))) %>%
        # Divide samples into quartiles by signal strength
        group_by(which_ramp) %>%
        mutate(qrt = as.character(ntile(mean_mfi, 4))) %>%
        group_by(which_ramp, qrt) %>%
        count(int_evid) %>%
        ggplot(aes(x = qrt, y = n, fill = int_evid)) +
        geom_bar(stat = "identity", show.legend = ifelse(capture == "epitope", T, F)) +
        scale_fill_manual(values = bar_col) +
        facet_wrap(~ which_ramp) +
        labs(x = paste0(which_expr[1], " expression quartile"),, y = "GPCR count", fill = "Evidence of\ninteraction",
             title = paste0(ifelse(capture == "hpa", "HPA", "Epitope"), " capture")) +
        coord_flip() +
        theme_minimal() +
        theme(axis.text = element_text(size = 14),
              axis.title = element_text(size = 14),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 14),
              strip.text = element_text(size = 14),
              plot.title = element_text(size = 16))
      
      return(plt_out)
    }))
  
  return(plts %>% pull(plt) %>% setNames(plts$capture))
}

quartile_plt_per_capt <- function(long_dat_all, cxdet_epitope, cxdet_hpa, which_expr = c("RAMP", "GPCR")) {
  # long_dat_all, long format data with anti-HA antibodies and OLLAS detection
  # cxdet_epitope, data frame with **raw interaction results** from epitope capture complex detection
  # cxdet_hpa, data frame with **interaction evidence strength** from HPA capture complex detection
  # which_expr, string indicating what to show expression of, RAMP or GPCR (takes first element)
  
  int_df <- rbind(
    cxdet_epitope %>% select(which_gpcr, which_ramp, capt, pass),
    cxdet_hpa %>% rename("pass" = "int_evid") %>% mutate(capt = "GPCR") %>% select(which_gpcr, which_ramp, capt, pass)
  )
  
  plts <- long_dat_all %>%
    filter(if (which_expr[1] == "RAMP") {gene_name == "M anti-HA"} else if (which_expr[1] == "GPCR") {gene_name == "M anti-FLAG"},
           if (which_expr[1] == "RAMP") {detection_ab == "OLLAS"} else if (which_expr[1] == "GPCR") {detection_ab == "1D4"},
           sample_type != "gpcr.p11", !interactor %in% c("mock", "empty")) %>%
    rename("which_gpcr" = "gpcr", "which_ramp" = "interactor") %>%
    select(unique_sample_name, which_gpcr, which_ramp, value) %>%
    # Get mean signal for each GPCR-RAMP combination, in case there are multiple samples
    group_by(which_gpcr, which_ramp) %>%
    summarise(mean_mfi = mean(value, na.rm = T), .groups = "keep") %>%
    ungroup() %>%
    right_join(int_df, by = c("which_gpcr", "which_ramp")) %>%
    # Bunch together RAMP capture as they will be separated by facets anyway
    mutate(capt = str_replace(capt, "RAMP\\d", "RAMP")) %>%
    group_by(capt) %>%
    nest() %>%
    mutate(plt = map(data, ~ {
      if (capt == "GPCR") {pass_lvls = c("Strong", "Medium", "Weak")} else {pass_lvls = c("Pass", "Fail")}
      
      .x %>%
        filter(!is.na(pass), !is.na(mean_mfi)) %>%
        mutate(pass = factor(pass, levels = pass_lvls)) %>%
        group_by(which_ramp) %>%
        # Add quartiles
        mutate(qrt = as.character(ntile(mean_mfi, 4)),
               qrt = factor(qrt, levels = as.character(1:4))) %>%
        group_by(which_ramp, qrt) %>%
        count(pass) %>%
        ggplot(aes(x = qrt, y = n, fill = pass)) +
        geom_bar(stat = "identity") +
        facet_wrap(~ which_ramp) +
        scale_fill_brewer(palette = ifelse(unique(capt) == "GPCR", "Dark2", "Set2"), direction = 1) +
        labs(x = NULL, y = "GPCR count",
             fill = ifelse(capt == "GPCR", "Interaction\nHPA capture", "Interaction\nepitope capture"), title = paste0(capt, " capture")) +
        coord_flip() +
        theme_minimal() +
        theme(axis.text = element_text(size = 10),
              axis.title = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              strip.text = element_text(size = 10),
              plot.title = element_text(size = 12, hjust = 0.5))
    }))
  
  plts %<>% pull(plt) %>% setNames(plts$capt)
  # plts <- wrap_plots(plts[c("OLLAS", "1D4", "HA", "FLAG", "RAMP", "GPCR")], ncol = 2, guides = "collect")
  
  plts <- 
    (ggplot(data.frame(a = paste0(which_expr[1], " expression quartile"), x = 1, y = 1)) +
       geom_text(aes(x = x, y = y, label = a), angle = 90, size = 5) +
       theme_void() +
       coord_cartesian(clip = "off")) +
    (plts$OLLAS / plts$HA / plts$RAMP) +
    (plts$`1D4` / plts$FLAG / plts$GPCR) +
    plot_layout(widths = c(1, 20, 20), guides = "collect")
  
  return(plts)
}

# Classify observations as true positive, true negative, false positive, false negative
classify <- function(obs, ref) {
  pn <- case_when(obs & ref ~ "TP",
                  !obs & !ref ~ "TN",
                  obs & !ref ~ "FP",
                  !obs & ref ~ "FN")
  return(list("sens" = sum(pn == "TP") / sum(pn %in% c("TP", "FN")),
              "spec" = sum(pn == "TN") / sum(pn %in% c("TN", "FP"))))
}

# Get expression and expected interactions for GPCRs, per GPCR (not per Ab) assayed, in the format used for epitope capture complex detection
cxdet_prep <- function(gpcrs, dat_list) {
  # GPCR expression
  expr_list <- dat_list$long_dat_all %>% group_by(gpcr) %>%
    summarise(int_exp = case_when(any(gpcr_expr_ns_flag_capture_1d4_detect == "ns") ~ "uncertain", T ~ "supported")) %>%
    filter(gpcr %in% gpcrs)
  
  int_exp <- dat_list$long_dat_all %>%
    filter(sample_type != "gpcr.p11") %>%
    filter(!gpcr %in% c("empty", "mock")) %>%
    select(gpcr, contains("int_ref")) %>%
    distinct() %>%
    # If multiple rows, keep the one with the most entries
    mutate(nna = rowSums(is.na(select(., -1)))) %>%
    group_by(gpcr) %>%
    filter(nna == min(nna)) %>%
    select(-nna) %>%
    pivot_longer(cols = -gpcr) %>%
    mutate(ramp = toupper(str_extract(name, "ramp\\d"))) %>%
    # Summarise literature interactions per GPCR and RAMP
    group_by(gpcr, ramp) %>%
    summarise(yes_int = !is.na(value[str_detect(name, "^int")]),
              no_int = !is.na(value[str_detect(name, "^no")]),
              .groups = "keep") %>%
    # Some are mixed references of yes and no for interaction, classify as uncertain
    mutate(int_exp = case_when(yes_int & !no_int ~ "Expected",
                               !yes_int & no_int ~ "None",
                               yes_int & no_int ~ "Uncertain")) %>%
    ungroup() %>%
    select(gpcr, ramp, int_exp) %>%
    rename(which_gpcr = gpcr, which_ramp = ramp)
  
  return(list("expr_list" = expr_list$int_exp %>% setNames(expr_list$gpcr),
              "int_exp" = int_exp))
}

gpcrramp_sensspec <- function(dat_in, gpcrs_in, expr_in, int_exp_in, mk_plts=F) {
  # Compute cutoffs based on approximate intersection of sensitivity and specificity for
  # known GPCR-RAMP interactions (and non-interactions) from literature
  #
  # dat_in, long format data
  # gpcrs_in, names of GPCRs to look at
  # mk_plts, T to make a plot for each RAMP+capture combination showing the sensitivity, specificity and
  # cutoff. F (default) to not make plots which saves time as the cutoff is decided as soon as the
  # sensitivity and specificity intersect
  #
  # Returns a list containing the elements "cutoffs" and "interactions",
  # containing the cutoffs per capture and RAMP interaction, and the
  # classifications of each GPCR-RAMP interaction for the supplied GPCRs
  
  gpcr_df <- data.frame("which_gpcr" = rep(gpcrs_in, each = 3),
                        "which_ramp" = rep(paste0("RAMP", 1:3), length(gpcrs_in)),
                        "gpcr_expr" = rep(expr_in[gpcrs_in], each = 3)) %>%
    left_join(., int_exp_in, by = c("which_gpcr", "which_ramp"))
  
  # Get median signal of RAMP-containing sample for each GPCR, RAMP, and capture scheme combination
  med_sig <- dat_in %>%
    filter(gpcr %in% gpcrs_in) %>%
    filter((gene_name %in% c("1D4", "M anti-FLAG") & detection_ab == "OLLAS") |
             (gene_name %in% c("Ollas", "M anti-HA") & detection_ab == "1D4") |
             (gene_name == "RAMP1" & interactor == "RAMP1" & detection_ab == "1D4") |
             (gene_name == "RAMP2" & interactor == "RAMP2" & detection_ab == "1D4") |
             (gene_name == "RAMP3" & interactor == "RAMP3" & detection_ab == "1D4"),
           # gpcr %in% (gpcr_df %>% filter(gpcr_expr == "expressed") %>% pull(which_gpcr)),
           interactor %in% paste0("RAMP", 1:3)) %>%
    # If sample is to be excluded, turn its values into NA
    mutate(value = case_when(is.na(exclude) ~ value)) %>%
    group_by(gpcr, interactor, gene_name, detection_ab) %>%
    summarise(signal = median(value), .groups = "keep") %>%
    ungroup() %>%
    rename("which_gpcr" = "gpcr", "which_ramp" = "interactor",
           "capt" = "gene_name", "det" = "detection_ab") %>%
    mutate(capt = case_when(grepl("M anti-", capt) ~ gsub("M anti-", "", capt),
                            capt == "Ollas" ~ "OLLAS",
                            T ~ capt))
  
  # Per RAMP and capture scheme, try out cutoffs (step = 0.01) between slightly
  # below min value and slightly above max value
  # Classify TP/TN/FP/FN and calculate sensitivity and specificity for multiple cutoffs
  # Approximate intersection of sens and spec by taking mid-point between
  # cutoff just before intersection and cutoff just after intersection
  
  # Or if there is no interest in making plots of sens and spec, instead use a while loop,
  # start at the min and go up until specificity is higher than sensitivity, then take
  # mid-point of last and previous cutoffs (faster)
  cutoffs_percapt_perramp <- sapply(c("1D4", "FLAG", "OLLAS", "HA", "RAMP"), function(captab) {
    sapply(paste0("RAMP", 1:3), function(ramp) {
      
      # If RAMP capture, only look at interaction with the captured RAMP
      if (captab == "RAMP") {captab <- ramp}
      detab <- ifelse(captab %in% c("1D4", "FLAG"), "OLLAS", "1D4")
      
      temp_df <- gpcr_df %>%
        filter(int_exp %in% c("Expected", "None"),
               which_ramp == ramp) %>%
        # Recode to logical for convenience in computing TP/TN/FP/FN
        mutate(int_exp = int_exp == "Expected") %>%
        mutate(capt = captab, det = detab) %>%
        left_join(., med_sig, by = c("which_gpcr", "which_ramp", "capt", "det")) %>%
        # Some samples may have failed, remove those (will have NA signal)
        filter(!is.na(signal))
      
      # Step to iterate by
      stp <- 0.01
      
      # If make plot, try cutoffs from slightly below the minimum value to slightly above the minimum
      if (mk_plts) {
        # Df to keep results
        cutoff_df <- data.frame("which_ramp" = ramp,
                                "capt" = captab,
                                "cutoff" = seq(from = min(temp_df$signal) - stp,
                                               to = max(temp_df$signal) + stp,
                                               by = stp),
                                "sens" = NA, "spec" = NA)
        
        for (i in 1:nrow(cutoff_df)) {
          sens_spec <- classify(obs = temp_df$signal > cutoff_df$cutoff[i],
                                ref = temp_df$int_exp)
          
          cutoff_df[i, c("sens", "spec")] <- c(sens_spec$sens, sens_spec$spec)
        }
        
        # Get ideal cutoff from between cutoffs before and after the sens and spec intersect
        cutoff <- median(c(
          cutoff_df %>% filter(sens > spec) %>% slice(nrow(.)) %>% pull(cutoff),
          cutoff_df %>% filter(sens < spec) %>% slice(1) %>% pull(cutoff)
        ))
        
        # Output df for cutoff table
        sens_spec <- classify(obs = temp_df$signal > cutoff,
                              ref = temp_df$int_exp)
        out_df <- data.frame("which_ramp" = as.character(ramp),
                             "capt" = as.character(captab),
                             "det" = detab,
                             "cutoff" = cutoff,
                             "sens" = sens_spec$sens,
                             "spec" = sens_spec$spec)
        
        return(list("plt_dat" = cutoff_df,
                    "cutoff" = out_df))
        
        # Otherwise, stop when sensitivity < specificity
      } else {
        
        cutoff <- min(temp_df$signal) - stp
        sens_spec <- classify(obs = temp_df$signal > cutoff,
                              ref = temp_df$int_exp)
        
        while (sens_spec$sens > sens_spec$spec) {
          cutoff <- cutoff + stp
          
          sens_spec <- classify(obs = temp_df$signal > cutoff,
                                ref = temp_df$int_exp)
        }
        
        # Get mid-point between previous and current cutoffs as the cutoff to use
        cutoff <- median(c(cutoff - stp, cutoff))
        sens_spec <- classify(obs = temp_df$signal > cutoff,
                              ref = temp_df$int_exp)
        
        return(list("cutoff" = data.frame("which_ramp" = as.character(ramp),
                                          "capt" = as.character(captab),
                                          "det" = detab,
                                          "cutoff" = cutoff,
                                          "sens" = sens_spec$sens,
                                          "spec" = sens_spec$spec)))
      }
    }, simplify = F)
  }, simplify = F)
  
  # Classify all GPCRs using the cutoffs
  cutoffs <- lapply(cutoffs_percapt_perramp, function(per_capt) {
    lapply(per_capt, function(per_ramp) {
      per_ramp$cutoff
    }) %>% rbindlist()
  }) %>% rbindlist()
  
  med_sig %<>%
    left_join(., cutoffs %>% select(-sens, -spec), by = c("which_ramp", "capt", "det")) %>%
    mutate(pass = as.numeric(signal > cutoff)) %>%
    left_join(., gpcr_df, by = c("which_gpcr", "which_ramp"))
  
  # If there is FZD in there, the FLAG and HA capture should be NA as the tags are flipped
  med_sig[str_detect(med_sig$which_gpcr, "FZD") & med_sig$capt %in% c("FLAG", "HA"), c("signal", "pass")] <- NA
  
  # Number of passing schemes
  int_n <- med_sig %>%
    group_by(which_gpcr, which_ramp) %>%
    summarise(n_pass = ifelse(all(is.na(pass)), NA, sum(pass, na.rm = T)))
  
  # Interaction evidence (strong, medium, weak)
  # Account for FZD and excluded samples that have either two or all values as NA
  int_evid <- med_sig %>%
    group_by(which_gpcr, which_ramp) %>%
    summarise(fraction_pass = ifelse(all(is.na(pass)), NA,
                                 sum(pass, na.rm = T) / sum(!is.na(pass))),
              int_evid = ifelse(all(is.na(pass)), NA,
                                case_when(fraction_pass > 2/3 ~ "Strong",
                                          fraction_pass > 1/3 ~ "Medium",
                                          fraction_pass <= 1/3 ~ "Weak")),
              .groups = "keep") %>%
    ungroup()
  
  # Output depending on plot making
  if (mk_plts) {
    
    plt_dat <- lapply(cutoffs_percapt_perramp, function(per_capt) {
      lapply(per_capt, function(per_ramp) {
        per_ramp$plt_dat
      }) %>% rbindlist()
    }) %>% rbindlist()
    
    # Plot sens & spec for each combination of RAMP interaction and capture.
    # For RAMPs, bunch together capture as RAMP instead of having one row each for
    # RAMP1, RAMP2, RAMP3 (as that makes a 3x3 grid with only the diagnoal filled in, rest are empty plots)
    plt_out <- ggplot(plt_dat %>%
                        pivot_longer(cols = c(sens, spec)) %>%
                        mutate(name = case_when(name == "sens" ~ "Sensitivity",
                                                name == "spec" ~ "Specificity"),
                               capt = case_when(grepl("RAMP", capt) ~ "RAMP",
                                                T ~ capt))) +
      geom_line(aes(x = cutoff, y = value, group = name, colour = name), size = 1) +
      geom_vline(data = cutoffs %>%
                   mutate(capt = case_when(grepl("RAMP", capt) ~ "RAMP",
                                           T ~ capt)),
                 aes(xintercept = cutoff, linetype = "Threshold"), size = 0.7) +
      ggh4x::facet_grid2(rows = vars(capt), cols = vars(which_ramp), scales = "free_x", drop = T, independent = "x") +
      scale_colour_brewer(palette = "Paired") +
      scale_linetype_manual(values = c("Threshold" = "dotted"), name = NULL) +
      scale_x_continuous(n.breaks = 3) +
      labs(x = "Robust Z-score", y = NULL, colour = NULL) +
      theme_classic() +
      theme(axis.line = element_line(size = 1),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 16),
            legend.text = element_text(size = 14),
            legend.title = NULL,
            strip.text = element_text(size = 14))
    
    return(list("plot_data" = plt_dat,
                "plots" = plt_out,
                "interactions" = med_sig,
                "interaction_n" = int_n,
                "interaction_evidence" = int_evid,
                "cutoffs" = cutoffs))
  } else {
    
    return(list("interactions" = med_sig,
                "interaction_n" = int_n,
                "interaction_evidence" = int_evid,
                "cutoffs" = cutoffs)) 
  }
}

# Function for making a heatmap of pass status for each GPCR + capture combination, per RAMP
pass_hm <- function(dat_in,
                    row_label_size=3,
                    draw_hm=T,
                    plt_title=NULL,
                    pd=unit(c(2, 2, 2, 2), "mm")) {
  # dat_in, cxdet_cd_separate data frame with the desired GPCRs to display
  # row_label_size, number for the row label size
  # draw_hm, logical indicating whether the output heatmap should be drawn (T) or
  # given in a list (F) together with its legend
  
  # Colours for pass/not pass
  col_vec <- RColorBrewer::brewer.pal(3, "Set2")[-3]
  names(col_vec) <- c("Pass", "Fail")
  
  hm <- sapply(paste0("RAMP", 1:3), function(ramp) {
    plt_dat <- dat_in %>%
      filter(which_ramp == ramp) %>%
      select(which_gpcr, capt, pass, int_exp) %>%
      # Shorter column names and categorical fill
      mutate(capt = gsub("[[:space:]].*", "", capt),
             pass = case_when(pass == 1 ~ "Pass", pass == 0 ~ "Fail")) %>%
      pivot_wider(id_cols = c("which_gpcr", "int_exp"),
                  names_from = "capt", values_from = "pass") %>%
      as.data.frame() %>%
      column_to_rownames("which_gpcr")
    
    # Expected interactions to use for annotating the heatmap
    ramp_int_exp <- plt_dat$int_exp
    
    Heatmap(
      plt_dat %>% select(`1D4`, FLAG, HA, OLLAS, !!ramp) %>% as.matrix(),
      cluster_rows = F, cluster_columns = F, col = col_vec, na_col = "white",
      show_heatmap_legend = F, layer_fun = function(j, i, x, y, width, height, fill) {
        if (any(!is.na(ramp_int_exp))) {
          border_lty <- case_when(ramp_int_exp[i] == "Expected" ~ "solid",
                                  ramp_int_exp[i] == "Uncertain" ~ "dashed",
                                  ramp_int_exp[i] == "None" ~ "dotted")
          draw_border <- !is.na(ramp_int_exp[i])
          grid.rect(x = x[draw_border], y = y[draw_border], width = width[draw_border], height = height[draw_border],
                    gp = gpar(lty = border_lty[draw_border], fill = NA, col = "black"))
        }
      }, column_title = ramp, column_title_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 4),
      show_row_names = ifelse(ramp == "RAMP1", T, F), row_names_side = "left", row_names_gp = gpar(fontsize = row_label_size)
    )
    
  }, simplify = F)
  
  # Legends for colours and cell outlines
  hm_legend <- list(
    Legend(labels = names(col_vec), title = "Pass", type = "grid", legend_gp = gpar(fill = col_vec)),
    Legend(labels = c("Expected", "None", "Uncertain"), title = "Interaction in\nliterature", type = "lines", legend_gp = gpar(col = "black", lty = c("solid", "dotted", "dashed")))
  )
  
  # Combine the RAMP heatmaps into one
  hm_list <- NULL
  for (i in 1:length(hm)) {
    hm_list <- hm_list + hm[[i]]
  }
  
  # Output
  if (draw_hm) {
    # Draw the heatmap with the legend
    draw(hm_list, annotation_legend_list = hm_legend, padding = pd,
         column_title = plt_title, column_title_gp = gpar(fontface = "bold", fontsize = 10))
  } else {
    return(list("heatmap" = hm_list,
                "legend" = hm_legend))
  }
}


cxdet_dens <- function(long_dat_all, det_res, var_measure = mad, n_var = 6, make_plots = F) {
  # long_dat_all, contains interaction data
  # det_res, Ab validation results
  
  # Pick out Abs to look at
  abs_to_check <- det_res %>% filter(pass) %>% pull(gene_name_id) %>% unique()
  
  cxdet_res <- long_dat_all %>%
    filter(detection_ab == "OLLAS" & sample_type != "gpcr.p11") %>%
    filter(gene_name_id %in% abs_to_check) %>%
    select(unique_sample_name, gpcr_interactor, gpcr, interactor, gene_name_id, gene_name, antibody, value, target, exclude) %>%
    # Don't use pos ctrl CALCRL as on-target
    mutate(target = case_when(gpcr_interactor == "CALCRL.RAMP1.posctrl" ~ "CALCRL.RAMP1.posctrl", T ~ target)) %>%
    # NA value for samples to exclude
    mutate(value = case_when(is.na(exclude) ~ value)) %>%
    # Names containing the Ab names
    mutate(gene_name_ab = paste0(gene_name, "_", antibody)) %>%
    group_by(gene_name_ab, gene_name_id) %>%
    nest() %>%
    # Cutoff per Ab
    mutate(cutoff = map(data, ~ {
      # Get density, peak and distance from peak position for computing variance measure and cutoff
      val <- .x %>% filter(!is.na(value)) %>% pull(value)
      trgt <- .x %>% filter(!is.na(value)) %>% pull(target)
      dens <- density(val)
      peak <- dens$x[which.max(dens$y)]
      
      pdist <- abs(val - peak)
      neg_prop <- 1 - sum(trgt == "Target") / sum(trgt != "Target")
      peak_var <- var_measure(val[pdist <= quantile(pdist, probs = neg_prop)])
      
      cutoff <- peak + n_var * peak_var
      
      return(cutoff)
    }) %>% unlist()) %>%
    # Interaction classification per Ab
    mutate(cls = map(data, ~ {
      .x %>%
        filter(target == "Target") %>%
        group_by(interactor) %>%
        summarise(mean_signal = mean(value),
                  int_pass = mean_signal > cutoff) %>%
        arrange(interactor) %>%
        pull(int_pass) %>%
        paste0(collapse = "_")
    })) %>%
    # Split interaction classifications into their own columns
    separate(cls, paste0("RAMP", 1:3), sep = "_", convert = T) %>%
    # Make plots if desired
    mutate(plt = map(data, ~ {
      if (!make_plots) {return(NULL)}
      # target_gpcr <- .x %>% filter(target == "Target") %>% pull(gpcr) %>% unique() %>% paste0(collapse = "/")
      
      .x %>%
        mutate(target = case_when(target == "Target" & interactor == "RAMP1" ~ "RAMP1",
                                  target == "Target" & interactor == "RAMP2" ~ "RAMP2",
                                  target == "Target" & interactor == "RAMP3" ~ "RAMP3",
                                  target == "Other" ~ "Other GPCR",
                                  T ~ target)) %>%
        arrange(target) %>%
        ggplot(aes(x = 0, y = value, colour = target)) +
        geom_beeswarm(priority = "none") +
        scale_colour_manual(values = c("#BFEFFF", "#B1E063", "#EEAD0E", "grey80", "red") %>%
                              setNames(c(paste0("RAMP", 1:3), "Other GPCR", "CALCRL.RAMP1.posctrl"))) +
        geom_hline(yintercept = cutoff, linetype = "dotted") +
        labs(y = "Robust Z-score", x = NULL, title = gene_name_ab, colour = "Sample") +
        theme_classic(16) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    })) %>%
    ungroup() %>%
    arrange(gene_name_ab)
  
  # Evidence strength for interaction per GPCR
  # For the multi-target Ab targeting HCAR2 + HCAR3, count it for both
  # For the multi-target Abs targeting P2RY11 + PPAN-P2RY11 and NPY4R + NPY4R2, only one of each target is present in the assay
  int_evid <- cxdet_res %>%
    select(gene_name_id, contains("RAMP")) %>%
    pivot_longer(cols = -gene_name_id, names_to = "which_ramp", values_to = "pass")
  
  hcar <- int_evid %>% filter(str_detect(gene_name_id, "HCAR2\\.HCAR3")) %>%
    mutate(gene_name_id = str_remove(gene_name_id, "HCAR3\\."))
  
  int_evid %<>%
    mutate(gene_name_id = case_when(str_detect(gene_name_id, "HCAR2\\.HCAR3\\.") ~ str_remove(gene_name_id, "HCAR2\\."),
                                    str_detect(gene_name_id, "NPY4R\\.NPY4R2\\.") ~ str_remove(gene_name_id, "NPY4R2\\."),
                                    str_detect(gene_name_id, "PPAN\\.P2RY11\\.P2RY11\\.") ~ str_remove(gene_name_id, "PPAN\\.P2RY11\\."),
                                    T ~ gene_name_id)) %>%
    rbind(., hcar) %>%
    arrange(gene_name_id, which_ramp) %>%
    mutate(which_gpcr = str_remove(gene_name_id, "\\.\\d+")) %>%
    group_by(which_gpcr, which_ramp)
  
  # Keep a table with n passing Abs
  int_n <- int_evid %>%
    summarise(n_pass = sum(pass))
  
  int_evid %<>%
    summarise(fraction_pass = sum(pass) / n(),
              int_evid = case_when(fraction_pass >= 2/3 ~ "Strong",
                                   fraction_pass < 2/3 & fraction_pass >= 1/3 ~ "Medium",
                                   fraction_pass < 1/3 ~ "Weak"),
              .groups = "keep") %>%
    ungroup()
  
  return(list("interactions" = cxdet_res,
              "interaction_n" = int_n,
              "interaction_evidence" = int_evid))
}


dens_pass_hm <- function(x, int_exp, col_name_bot = T) {
  # x, df containing interaction info in separate columns (named RAMP1, RAMP2, RAMP3) for each Ab (gene_name_ab)
  # int_exp, whether corresponding interactions are expected according to literature
  # col_name_bot, boolean for whether the column names should be at the bottom or not (if not they will be at the top)
  
  plt_dat <- x %>%
    # Format for plotting
    select(gene_name_ab, contains("RAMP")) %>%
    arrange(gene_name_ab) %>%
    column_to_rownames("gene_name_ab") %>%
    mutate(across(everything(), as.numeric)) %>%
    as.matrix()
  
  # Colours for cells
  col_vec <- RColorBrewer::brewer.pal(3, "Set2")[-3]
  names(col_vec) <- c("1", "0")
  
  hm <- Heatmap(plt_dat, cluster_rows = F, cluster_columns = F, col = col_vec, row_names_side = "left", show_heatmap_legend = F,
                name = "Interaction", row_names_gp = gpar(fontsize = 4), cell_fun = function(j, i, x, y, width, height, fill) {
                  if (!is.na(int_exp[match(rownames(plt_dat)[i], int_exp$gene_name_ab),
                                     match(colnames(plt_dat)[j], colnames(int_exp))])) {
                    cell_outline <- switch(int_exp[match(rownames(plt_dat)[i], int_exp$gene_name_ab),
                                                   match(colnames(plt_dat)[j], colnames(int_exp)), drop = T],
                                           "Expected" = "solid", "None" = "dotted", "Uncertain" = "dashed")
                    grid.rect(x = x, y = y, width = width, height = height,
                              gp = gpar(lty = cell_outline, fill = NA, col = "black"))
                  }
                }, row_split = int_exp[match(rownames(plt_dat), int_exp$gene_name_ab), "gene_name", drop = T],
                row_title = NULL, row_gap = unit(1, "points"), na_col = "white",
                column_names_gp = gpar(fontsize = 5), column_names_side = ifelse(col_name_bot, "bottom", "top"))
  
  hm_legend <- list(
    Legend(labels = names(col_vec), title = "Pass", type = "grid", legend_gp = gpar(fill = col_vec)),
    Legend(labels = c("Expected", "None", "Uncertain"), title = "Interaction in\nliterature", type = "lines", legend_gp = gpar(col = "black", lty = c("solid", "dotted", "dashed")))
  )
  
  return(list("heatmap" = hm, "legend" = hm_legend))
}


# UpSet plot to compare interaction evidence across RAMPs (not very good)
cxdet_ramp_upset <- function(x, class_order = c("strength", "ramp"), ...) {
  # x, data_frame containing the interaction evidence strengths
  # class_order, character saying which order the classes/sets should be in, ordered by evidence strength ("strength") or ordered by RAMP ("ramp)
  # Other arguments to the ComplexUpset upset function
  
  plt_dat <- x %>%
    select(which_gpcr, which_ramp, int_evid) %>%
    group_by(which_gpcr, which_ramp) %>%
    summarise(strong = int_evid == "Strong",
              medium = int_evid == "Medium",
              weak = int_evid == "Weak",
              .groups = "keep") %>%
    pivot_wider(id_cols = which_gpcr, names_from = which_ramp, values_from = c(strong, medium, weak)) %>%
    rename_with(.fn = function(x) {
      # Skip first element as it is to be left untouched
      strength <- str_extract(x[-1], "strong|medium|weak")
      ramp <- str_extract(x[-1], "RAMP\\d")
      return(c(x[1], paste0(ramp, "_", strength)))
    })
  
  # Make order of classes
  sets <- switch(class_order[1],
                 "strength" = paste0(paste0("RAMP", 3:1), "_", rep(c("weak", "medium", "strong"), each = 3)),
                 "ramp" = paste0(rep(paste0("RAMP", 3:1), each = 3), "_", c("weak", "medium", "strong")))
  
  plt <- upset(plt_dat, intersect = sets, sort_sets = F, wrap = T, ...)
  
  return(plt)
}

# Check hits from validation experiments
valid_hits <- function(valid_data, use_mod_z = T) {
  # valid_data, list containing the data from the validation experiments (as loaded by the data loading function)
  # use_mod_z, logistic, choice of ordinary Z scores or modified Z scores
  
  prot <- valid_data$data_list$prot
  well_info <- valid_data$data_list$well_info
  mfi <- valid_data$data_list$mfi %>%
    left_join(well_info %>% select(samp_name, detection_ab, cell_type, exclude, gpcr, ramp), by = "samp_name") %>%
    # Skip empty samples and samples to exclude
    filter(cell_type != "buffer", is.na(exclude)) %>%
    unite("gpcr_ramp", gpcr, ramp) %>%
    select(-exclude)
  
  well_info <- filter(well_info, samp_name %in% mfi$samp_name)
  
  # Divide data into data without CALCRL (pos ctrl) and data with only CALCRL
  mfi_clr <- mfi %>%
    filter(str_detect(gpcr_ramp, "CALCRL") |
             str_detect(samp_name, "buffer"))
  
  mfi_all <- mfi
  
  mfi <- mfi %>%
    filter(str_detect(gpcr_ramp, "CALCRL", negate = T)) %>%
    # Add interaction type tag showing if the data is for endogenous (nothing transfected) or semi-endogenous interactions (RAMP transfected)
    mutate(int_type = case_when(str_detect(gpcr_ramp, "RAMP") ~ "Semi-endogenous", T ~ "Endogenous")) %>%
    # For SK-N-MC and SHSY-5Y cells, look only at endogenous interactions
    filter(!(cell_type %in% c("SK-N-MC", "SH-SY5Y") & gpcr_ramp != "untransfected_untransfected")) %>%
    # Skip OLLAS detection for untransfected/mock-mock cells as they have no OLLAS to bind
    filter(!(gpcr_ramp %in% c("untransfected_untransfected", "mock_mock") & detection_ab == "OLLAS"))
  
  # Scaling function, depending on selection of z or mod z
  scale2 <- function(x, mod_z) {
    if (mod_z) {
      # Scale long format data so just scaling one column at a time, no need for loops
      x_out <- (x - median(x, na.rm = T)) / mad(x, na.rm = T)
    } else {
      x_out <- scale(x)
    }
    
    return(x_out)
  }
  
  hit_check <- mfi %>%
    # If semi-endogenous and RAMP detection in Expi cells, keep just the samples with matching RAMP transfection
    filter(case_when(cell_type == "Expi" & int_type == "Semi-endogenous" & detection_ab != "OLLAS" ~ str_detect(gpcr_ramp, detection_ab), T ~ T)) %>%
    pivot_longer(cols = -c(samp_name, detection_ab, cell_type, gpcr_ramp, int_type), names_to = "binder", values_to = "mfi") %>%
    # Skip tags, RAMPs, and empty samples
    filter(!str_detect(binder, "OLLAS|RAMP|1D4|anti\\.HA|anti\\.FLAG|EMPTY")) %>%
    # Skip HUVEC samples, will be treated differently due to low number of samples
    filter(cell_type != "HUVEC") %>%
    arrange(cell_type, detection_ab, gpcr_ramp) %>%
    group_by(cell_type, detection_ab, gpcr_ramp) %>%
    nest() %>%
    mutate(norm_dat = map(data, ~ {
      # Perform quantile normalisation
      # Transpose as columns are assumed to be samples
      norm_out <- .x %>% select(-int_type) %>% mutate(mfi = log2(mfi)) %>%
        pivot_wider(id_cols = "samp_name", names_from = "binder", values_from = "mfi") %>%
        column_to_rownames("samp_name") %>% t()
      norm_out <- normalizeBetweenArrays(object = norm_out, method = "quantile")
      # Transform back into original format
      norm_out <- t(norm_out) %>% as.data.frame() %>% rownames_to_column("samp_name") %>%
        pivot_longer(cols = -samp_name, names_to = "binder", values_to = "mfi") %>%
        # Add binder target names
        mutate(gene_name = prot[match(binder, prot$gene_name_id), "Gene", drop = T]) %>%
        mutate(mfi = 2**mfi)
      
      return(norm_out)
    })) %>%
    # Beeswarm plots showing the cutoff
    mutate(plt = map(norm_dat, ~ {
      # Convert to Z-scores and mark Z-score of 3 (mean + 3sd) or robust Z-score of 3.5
      plt_dat <- .x %>% mutate(mfi = scale2(mfi, use_mod_z)) %>%
        # Points above threshold marked in data, but skip labelling due to cluttering. Keep marks as they are used by other code for picking out hits
        group_by(binder) %>%
        summarise(mean_z = mean(mfi)) %>%
        mutate(point_lab = case_when(mean_z > ifelse(use_mod_z, 3.5, 3) ~ binder, T ~ "")) %>%
        ungroup()
      
      ggplot(plt_dat, aes(x = 0, y = mean_z)) +
        geom_beeswarm(size = 0.5) +
        geom_hline(yintercept = ifelse(use_mod_z, 3.5, 3), linetype = "dotted", colour = "red") +
        labs(x = NULL, y = ifelse(use_mod_z, "Mean robust Z-score", "Mean Z-score"),
             title = paste0(cell_type, " cells")) +
        theme_classic() +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    })) %>%
    # Tables with the values
    mutate(tbl = map(norm_dat, ~ {
      .x %>% mutate(mfi = scale2(mfi, use_mod_z)) %>%
        # Label points above threshold by Ab
        group_by(binder) %>%
        summarise(gene_name = unique(gene_name),
                  mean_z = mean(mfi)) %>%
        mutate(point_lab = case_when(mean_z > ifelse(use_mod_z, 3.5, 3) ~ binder, T ~ "")) %>%
        ungroup()
    })) %>%
    select(-data) %>%
    mutate(hits = map(tbl, ~ {
      .x %>% filter(point_lab != "") %>% pull(gene_name) %>% unique()
    }))
  
  return(hit_check)
} 



