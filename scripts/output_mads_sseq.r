###################################################################################
# This script takes an xls with the columns "Originalnr.", "Primer" and "Rør nr",
# and a table with the results of the csc-program.
# The tables are joined and ultimately output i a format compatible with mads.


library(tidyverse)


args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())

# Read cli-args
#batch = args[1]
metadata_file = args[1]
csc_results = args[2]



devel = F
# rm(list = ls()); devel = T


if (devel) {
    
    batch = "210216"
    
    metadata_file = paste0("~/GenomeDK/clinmicrocore/sanger/input/", batch, ".xls")
    csc_results_file = paste0("~/GenomeDK/clinmicrocore/sanger/output/", batch, "/csc/", batch, "_results.csv")
    
}




#metadata = read_tsv(metadata_file)
metadata_read = readxl::read_excel(metadata_file, sheet = 1, na = c("", "?"))

metadata_read %>% glimpse
metadata = metadata_read %>% 
    select(kma_raw_sample_name = `Originalnr.`, Primer, ef_sample_name = `Rør nr`) %>% 
    filter(!is.na(kma_raw_sample_name)) %>% 

    # Mark the sample type
    # Should be synced with pipe19/scripts/parse_path.r:
    mutate(type = case_when(str_detect(tolower(kma_raw_sample_name), "positiv|pos$|seqpos") ~ "positive_control",
                            str_detect(tolower(kma_raw_sample_name), "negativ|h2o|^empty|blank|tom|^neg|neg$") ~ "negative_control",
                            str_detect(tolower(kma_raw_sample_name), "^afd|^00") ~ "other",
                            TRUE ~ "sample")) %>% 
    
    # convert the numbers to sample types
    mutate(sample_name_prefix = str_sub(kma_raw_sample_name, 1, 2), # first two characters
           sample_name_suffix = str_sub(kma_raw_sample_name, -6),   # last six characters
           #raw_full_name = paste0(batch, ".", plate, ".", moma_serial, "_", kma_raw_sample_name),
           sample_name_prefix_converted = recode(sample_name_prefix,
                                                 "87" = paste0("L"),
                                                 "88" = paste0("R"),
                                                 "89" = paste0("I"),
                                                 "90" = paste0("V"),
                                                 "96" = paste0("P"))) %>% 
    
    rowwise() %>% 
    mutate(kma_ya_sample_name = if_else(type == "sample", # "year_agnostic_sample_name"
                                    paste0(sample_name_prefix_converted, sample_name_suffix),
                                    paste0(kma_raw_sample_name))) %>% 
    ungroup() %>% 
    
    select(kma_ya_sample_name, kma_raw_sample_name, ef_sample_name, primer = Primer, type)


    





results_read = read_csv(csc_results_file, col_types = cols(.default = "c")) 
results_read %>% glimpse
results = results_read %>% 
    mutate(ef_sample_name = str_remove(sample, "_.+$")) %>% 
    select(ef_sample_name, ef_raw_sample_name = sample, comment, everything()) %>% 
    


    # Prefix variant-position-columns with "pos_"
    # Be aware that this call is dependent on the number, and ordering of columns
    # The forth column should be the leftmost variant column.
    rename_at(.vars = vars(4:last_col()),
              .funs = ~ paste0("pos_", .x))


# Join the data together
joined = metadata %>% 
    left_join(results, by = "ef_sample_name") %>% 
    arrange(kma_ya_sample_name)
    
    # Mark the overlapping samples
    #group_by(kma_ya_sample_name, type, comment) %>% 
    #mutate(rank = row_number(primer)) 

    
    

# Check out discrepancies:
joined %>% 
    summarize_at(vars(starts_with("pos_"), ef_sample_name), ~ list(unique(.x))) %>% View



# If you are happy with the discrepancies, you may proceed:
joined %>% 
    mutate_at(vars(starts_with("pos_")), ~  na.omit(.x)) %>% View
    



# Inform about missing data
missing_data = anti_join(metadata, results)
write(paste0("Number of KMA samples missing from EF data: ", nrow(missing_data), "\n",
            paste(c(" ~KMA:", missing_data$kma_raw_sample_name),
                  c("~EF:", missing_data$ef_sample),
                  collapse = "\n"), "\n"),
      stderr())




# Write the joined data out
joined %>% write_tsv(paste0("~/GenomeDK/clinmicrocore/sanger/output/", batch, "_joined.tsv")) %>% 

    









# Format the data for mads.

# 1: open the example sent to lene/helge
# 2: open the varianttable from ssi
joined %>% 
    
    
    
    
    
    