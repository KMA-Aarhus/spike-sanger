####################################################################
# This script takes a 


library(tidyverse)


args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

# Read cli-args
#batch = args[1]
csc_results = args[1]
metadata_file = args[2]


devel = F
# rm(list = ls()); devel = T


if (devel) {
    
    batch = "210216"
    metadata_file = paste0("~/GenomeDK/clinmicrocore/sanger/input/", batch, ".xls")
    csc_results_file = paste0("~/GenomeDK/clinmicrocore/sanger/output/", batch, "/csc/", batch, "_results.csv")
    
}




#metadata = read_tsv(metadata_file)
metadata = readxl::read_excel(metadata_file, sheet = 1, na = c("", "?"))

metadata %>% glimpse
metadata = metadata %>% 
    select(kma_sample = `Originalnr.`, Primer, ef_sample = `Rør nr`) %>% 
    filter(!is.na(kma_sample))





results = read_csv(csc_results_file, col_types = cols(.default = "c")) 
results %>% glimpse
results = results %>% 
    mutate(ef_sample = str_remove(sample, "_.+$")) %>% 
    select(ef_sample, ef_raw_sample = sample, comment, everything()) #%>% 
    
    #pivot_longer(c(everything(), -sample, -comment))






joined = metadata %>% left_join(results)


# Inform about missing data
missing_data = anti_join(metadata, results)
write(paste0("Number of KMA samples missing from EF data: ", nrow(missing_data), "\n",
            paste(c(" ~KMA:", missing_data$kma_sample),
                  c("~EF:", missing_data$ef_sample),
                  collapse = "\n"), "\n"),
      stderr())




# Write the joined data out










# Format the data for mads.