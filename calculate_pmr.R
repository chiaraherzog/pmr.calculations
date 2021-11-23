# PMR calculation
# Author: Chiara Herzog
# contact: chiara.herzog@uibk.ac.at
# Date 17 Nov 2021
# Type: gblock or sss1


calculate_pmr <- function(folder, output, experimentname,
                          write.results = TRUE,
                          type = "gblock"){
  
  # libraries
  require(dplyr)
  require(tidyverse)
  require(xlsx)
  
  # find files
  files <- list.files(folder)
  curve <- grep("Standard Curve", files)
  results <- grep("Results", files)
  
  if(length(curve) == 0){
    stop("Error: file with standard curves not found")
  }
  
  if(length(results) == 0){
    stop("Error: results file not found in provided folder")
  }
  
  # check if type is available
  if(type == "gblock"){
    cat("Starting calculation using gblocks...\n")
  } else if(type == "sss1"){
    cat("Starting calculation using Sss1-treated DNA...\n")
  } else {
    stop("Please select either type 'gblock' or 'sss1'.")
  }
  
  output_folder <- list.files(output)
  
  if(paste(experimentname, ".csv", sep = "") %in% output_folder == TRUE && write.results == TRUE){
    stop("File with the same name already exists in the output folder. Please choose a unique name.")
  }
  
  # Generate standard curve from COL2A1 values
  curve <- read.table(file = paste(folder, files[curve], sep = ""),
                      sep = ",",
                      header = TRUE)
  
  curve <- suppressWarnings(curve %>%
                              mutate(ct = as.numeric(Cq.Mean)))
  
  if(type == "gblock"){
    plot <- curve %>%
      filter(Target == "COL2A1" & ((Sample %in% c("STD1", "STD2", "STD3","STD4") | (Sample %in% c("STD_1", "STD_2", "STD_3", "STD_4")) | (Sample %in% c("Std 1", "Std 2", "Std 3", "Std 4"))))) %>%
      mutate(x = case_when(Sample %in% c("STD1", "STD_1", "Std 1") ~ 6,
                           Sample %in% c("STD2", "STD_2", "Std 2") ~ 5,
                           Sample %in% c("STD3", "STD_3", "Std 3") ~ 4,
                           Sample %in% c("STD4", "STD_4", "Std 4") ~ 3)) %>%
      mutate(ct = ifelse(Cq == "Undetermined", NA, as.numeric(Cq)))
  } else {
    plot <- curve %>%
      filter(Target == "COL2A1" & ((Sample %in% c("STD1", "STD2", "STD3","STD4") | (Sample %in% c("STD_1", "STD_2", "STD_3", "STD_4")) | (Sample %in% c("Std 1", "Std 2", "Std 3", "Std 4"))))) %>%
      mutate(x = case_when(Sample %in% c("STD1", "STD_1", "Std 1") ~ 250,
                           Sample %in% c("STD2", "STD_2", "Std 2") ~ 62.5,
                           Sample %in% c("STD3", "STD_3", "Std 3") ~ 15.625,
                           Sample %in% c("STD4", "STD_4", "Std 4") ~ 3.9)) %>%
      mutate(ct = ifelse(Cq == "Undetermined", NA, as.numeric(Cq)))
  }
  
  fit <- lm(ct~x, data = plot)
  anno <- paste0("R2 = ", signif(summary(fit)$r.squared, 3),
                 "\ny = ", round(fit$coefficients[2],2) , "x + ", round(fit$coefficients[1],2))
  y = max(plot$ct) - 0.05*max(plot$ct)
  
  plot <- plot %>%
    ggplot(aes(x = x,
               y = ct)) +
    geom_smooth(method = "lm",
                size = 0.5,
                se = FALSE,
                alpha = 0.3,
                colour = "gray40") +
    geom_point(aes(colour = Sample),
               size = 0.75) +
    theme(legend.title = element_blank(),
          legend.position = "top",
          legend.key = element_blank()) +
    xlab("Concentration") +
    ylab("Ct value")  +
    annotate("text",
             x = 5,
             y = y,
             label = anno,
             hjust = 0,
             size = 3)
  
  ggsave(plot,
         file = paste0(output,experimentname, ".png"),
         width = 5,
         height = 4) 
  cat("Curve plot saved under", paste0(output,experimentname, ".png\n"))
  
  
  model <- curve %>%
    filter(Target == "COL2A1" & ((Sample %in% c("STD1", "STD2", "STD3","STD4") | (Sample %in% c("STD_1", "STD_2", "STD_3", "STD_4")) | (Sample %in% c("Std 1", "Std 2", "Std 3", "Std 4"))))) %>%
    group_by(Sample) %>%
    summarise(ct = mean(ct),
              quantity = log10(mean(Quantity))) %>%
    mutate(quantity = case_when(Sample %in% c("STD_1", "STD1", "Std 1") ~ 6,
                                Sample %in% c("STD_2", "STD2", "Std 2") ~ 5,
                                Sample %in% c("STD_3", "STD3", "Std 3") ~ 4,
                                Sample %in% c("STD_4", "STD4", "Std 4") ~ 3))
  
  std <- lm(ct ~ quantity, data = model)
  intercept <- unname(std$coefficients[1])
  slope <- unname(std$coefficients[2])
  
  # Get & format data
  data <- read.table(file = paste(folder, files[results], sep = ""),
                     sep = ",",
                     header = TRUE)
  
  data <- data %>%
    mutate(rep = case_when(grepl("A|C|E|G|I|K|M|O", Well.Position) == TRUE ~ 1,
                           grepl("B|D|F|H|J|L|N|P", Well.Position) == TRUE ~ 2)) %>%
    mutate(ct = as.numeric(ifelse(Cq == "Undetermined", NA, Cq.Mean)))
  
  # If one rep is undetermined, set both to NA (undetermined)
  samples <- unique(data$Sample)
  targets <- unique(data$Target)
  reprocess <- data.frame(sample = "",
                          target = "")
  
  for (i in 1:length(samples)){
    for(j in 1:length(targets)){
      tmp <- data %>%
        filter(Sample == samples[i] & Target == targets[j]) %>%
        pull(ct)
      
      und <- sum(is.na(tmp))
      ind1 <- na.omit(data$Sample==samples[i] & data$Target==targets[j])
      ct <- ifelse(und<2, na.omit(tmp), NA)
      
      if(und > 0 && und < 2){
        x <- data.frame(sample = unique(data[ind1,]$Sample),
                        target = unique(data[ind1,]$Target))
        reprocess <- rbind(reprocess,x)
        cat("Sample ", samples[i], " has one undetermined Ct value for target ", targets[j], ". One PMR was carried forward for analysis for this target but reprocessing is recommended.\n", sep = "")
        data[ind1,]$ct <- ct
      } else if(und == 2){
        data[ind1,]$ct <- ct
      }
    }
  }
  
  reprocess <- reprocess %>%
    slice(-1)  %>%
    pivot_wider(id_cols = sample,
                values_from = target,
                names_from = target) %>%
    unite("targets", RALYL:ImC,
          na.rm = TRUE,
          sep = ", ")
  
  data <- data %>%
    pivot_wider(id_cols = Sample,
                names_from = Target,
                values_from = ct,
                values_fn = mean)
  
  # Warning: not amp
  not_amp <- data %>%
    filter(is.na(COL2A1) & Sample != "H2O")
  
  if(nrow(not_amp) > 0){
    cat(paste("COL2A1 was not amplified in ", nrow(not_amp), " samples (", paste(not_amp$Sample, collapse = ", "), ").\n", sep = ""))
  }  
  
  
  # Warning: high Cq
  high_col2a1 <- data %>%
    filter(COL2A1 > 30 & ! Sample %in% c("STD4")) %>%
    arrange(desc(Sample)) %>%
    mutate(id = paste(Sample, sep = " rep ")) 
  
  if(nrow(high_col2a1) > 0){
    cat(paste("COL2A1 had a Cq of > 30 in ", nrow(high_col2a1), " wells (", length(unique(high_col2a1$Sample)), " samples):\n", paste(high_col2a1$id, collapse = ", "), "\n", sep = ""))
  }  
  
  
  # Calculate input amounts
  data <- data %>%
    mutate(input_ZSCAN12 = 10^((ZSCAN12-intercept)/slope),
           input_GYPC1 = 10^((GYPC1-intercept)/slope),
           input_GYPC2 = 10^((GYPC2-intercept)/slope),
           input_DPP6 = 10^((DPP6-intercept)/slope),
           input_RALYL = 10^((RALYL-intercept)/slope),
           input_GSX1 = 10^((GSX1-intercept)/slope),
           input_COL2A1 = 10^((COL2A1-intercept)/slope),
           input_EpC = 10^((EpC-intercept)/slope),
           input_ImC = 10^((ImC-intercept)/slope)) %>%
    mutate(ref_ZSCAN12 = input_ZSCAN12/input_COL2A1,
           ref_GYPC1 = input_GYPC1/input_COL2A1,
           ref_GYPC2 = input_GYPC2/input_COL2A1,
           ref_DPP6 = input_DPP6/input_COL2A1,
           ref_RALYL = input_RALYL/input_COL2A1,
           ref_GSX1 = input_GSX1/input_COL2A1,
           ref_COL2A1 = input_COL2A1/input_COL2A1,
           ref_EpC = input_EpC/input_COL2A1,
           ref_ImC = input_ImC/input_COL2A1)
  
  gblock <- data %>%
    filter(Sample %in% c("gBlock", "gBLOCK", "gBlock (+)")) 
  
  
  # Calculate PMRs
  results <- data %>%
    mutate(ZSCAN12 = (ref_ZSCAN12/gblock$ref_ZSCAN12)*100,
           GYPC1 = (ref_GYPC1/gblock$ref_GYPC1)*100,
           GYPC2 = (ref_GYPC2/gblock$ref_GYPC2)*100,
           DPP6 = (ref_DPP6/gblock$ref_DPP6)*100,
           RALYL = (ref_RALYL/gblock$ref_RALYL)*100,
           GSX1 = (ref_GSX1/gblock$ref_GSX1)*100,
           EpC = (ref_EpC/gblock$ref_EpC)*100,
           ImC = (ref_ImC/gblock$ref_ImC)*100) %>%
    select(Sample, ZSCAN12, GYPC1, GYPC2, DPP6, RALYL, GSX1, EpC, ImC)
  
  results[is.na(results)] <- 0 # set all NAs to 0
  
  if(write.results==TRUE){
    results <- as.data.frame(results)
    write.xlsx(results,
               file = paste(output, experimentname, ".xlsx", sep = ""),
               sheetName = "Results",
               row.names = FALSE)
    
    if(!is_empty(reprocess)){
      write.xlsx(as.data.frame(reprocess),
                 file = paste(output, experimentname, ".xlsx", sep = ""),
                 sheetName = "Reprocessing recommended",
                 row.names = FALSE,
                 col.names = FALSE,
                 append = TRUE)
    }
    cat("Output file created.")
  }
  
  return(results)
}
