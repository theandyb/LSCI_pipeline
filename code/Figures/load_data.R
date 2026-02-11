# Paths
base_dir <- "output"
single_pos_dir <- paste0(base_dir, "/single_pos")
control_control_dir <- paste0(base_dir, "/single_pos_cc")

subtypes <- c("AT_CG", "AT_GC", "AT_TA",
    "GC_TA", "GC_CG", "GC_AT",
    "cpg_GC_TA", "cpg_GC_CG", "cpg_GC_AT")

# Load position-level results
get_pos_results <- function(res_dir, pop, subtype){
  f_name <- paste0(res_dir, "/", pop, "/", subtype, ".csv")
  df <- read_csv(f_name, show_col_types = FALSE) %>%
    rowwise() %>%
    mutate(re = dev / (2 * (singletons + controls)),
           is_sig = pchisq(dev, 3, lower.tail = FALSE) < 0.05,
           pop = pop)
  df
}

# Load position-level results for all populations
get_pos_results_all_pops <- function(res_dir, subtype){
    pops <- c("AFR", "AMR", "EAS", "EUR", "SAS", "ALL")
    results <- vector(mode = "list", length = length(pops))

    for(i in 1:length(pops)){
        results[[i]] <- get_pos_results(res_dir, pops[[i]], subtype)
    }

    bind_rows(results)
}

# Load position-level control-control results
get_pos_cc <- function(control_control, pop, subtype){
  f_name <- paste0(control_control, "/", pop, "/", subtype, "_v2.csv")
  df <- read_csv(f_name, show_col_types = FALSE) %>%
           mutate(re = dev / (2 * (singletons + controls)),
           is_sig = dev > qchisq(0.95, 3),
           pop = pop)
  df
}

# Load position-level results for multiple subtypes
get_pos_results_all <- function(res_dir, pop, subtypes){
  results <- vector(mode = "list", length = length(subtypes))

  for(i in 1:length(subtypes)){
    results[[i]] <- get_pos_results(res_dir, pop, subtypes[[i]])
    results[[i]]$subtype <- subtypes[[i]]
  }

  bind_rows(results)
}

# Get nice text for subtypes
subtype_print_names <- function(st){
  if(str_starts(st, "AT")){
    return(paste0("A → ", str_sub(st, 4, 4)))
  } else if(str_starts(st, "GC")){
    return(paste0("C → ", str_sub(st, 5, 5)))
  } else if(str_starts(st, "all")){
    return(paste0("C* → ", str_sub(st, 9, 9)))
  } else{
    return(paste0("CpG → ", str_sub(st, 9, 9), "pG"))
  }
}

# Load position-level residuals
get_pos_resid <- function(res_dir, pop, subtype, rp){
  f_name <- paste0(res_dir, "/resid/", pop, "/", subtype, "_rp_", rp, ".csv")
  df <- read_csv(f_name, show_col_types = FALSE)
  df$rp <- rp

  n_tot <- sum(df$n)
  df$re_res <- (df$res^2) / (2*n_tot)
  df
}

# Load position-level residuals for multiple positions
get_pos_resid_all <- function(res_dir, pop, subtype, rp = c(-500:-1, 1:500)){
  results <- vector(mode = "list", length = length(rp))
  for(i in 1:length(rp)){
    if(rp[i] == 1 & str_starts(subtype, "cpg")) next
    results[[i]] <- get_pos_resid(res_dir, pop, subtype, rp[i])
  }
  df <- bind_rows(results)
  df
}

# Summary across all subtypes
load_all_pos_results <- function(subtypes, pop = "ALL"){
  final <- get_pos_results(single_pos_dir, pop, subtypes[1]) %>%
    mutate(subtype = subtypes[1])
  for(i in 2:length(subtypes)){
    final <- final %>%
      bind_rows({
        get_pos_results(single_pos_dir, pop, subtypes[i]) %>%
          mutate(subtype = subtypes[i])
      })
  }
  final <- final %>%
    mutate(from_n = ifelse(str_starts(subtype, "A"), "A", 
                           ifelse(str_starts(subtype, "cpg"), "CpG", 
                                  ifelse(str_starts(subtype, "all"), "C*", "C"))),
         to_n = ifelse(str_starts(subtype, "A"), str_sub(subtype, 4,4),
                       ifelse(str_starts(subtype, "cpg"),
                              paste0(str_sub(subtype, -1), "pG"),
                              str_sub(subtype, -1))))
  final
}