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
# devel = T


if (devel) {
    rm(list = ls())
    devel = T
    
    arg_batch = "210214"
    
    metadata_file = paste0("~/GenomeDK/clinmicrocore/sanger/input/", arg_batch, ".xls")
    csc_results_file = paste0("~/GenomeDK/clinmicrocore/sanger/output/", arg_batch, "/csc/", arg_batch, "_results.csv")
    
}




#metadata = read_tsv(metadata_file)
metadata_read = readxl::read_excel(metadata_file, sheet = 1, na = c("", "?"))

metadata_read %>% glimpse
metadata = metadata_read %>% 
    select(kma_raw_sample_name = `Originalnr.`, primer_raw = Primer, ef_sample_name = `Rør nr`) %>% 
    filter(!is.na(kma_raw_sample_name)) %>% 
    mutate(primer = case_when(str_detect(primer_raw, "^L-") ~ "left_primer",
                              str_detect(primer_raw, "^R-") ~ "right_primer"),
           plate = str_extract(primer_raw, "\\d+$")) %>%

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
    
    select(kma_ya_sample_name, kma_raw_sample_name, ef_sample_name, plate, primer, type)



    



# Define position of interest
poi = c("Sseq_E484K", "Sseq_N501Y", "Sseq_Q677H", "Sseq_P681H")

# Read csc-results
results_read = read_csv(csc_results_file, col_types = cols(.default = "c")) 
results_read %>% glimpse
results = results_read %>% 
    mutate(ef_sample_name = str_remove(sample, "_.+$")) %>% 
    select(ef_sample_name, ef_raw_sample_name = sample, comment, everything()) %>% 
    


    # Prefix variant-position-columns with "Sseq_"
    # Be aware that this call is dependent on the number, and ordering of columns
    # The forth column should be the leftmost variant column.
    rename_at(.vars = vars(4:last_col()),
              .funs = ~ paste0("Sseq_", .x)) %>% 
  
    # To make matters simpler, we're only interested in keeping the poi columns.
    select(ef_sample_name, ef_raw_sample_name, comment, all_of(poi))



# Join the data together
joined = metadata %>% 
    left_join(results, by = "ef_sample_name") %>% 
    arrange(kma_ya_sample_name)
    
    # Mark the overlapping samples
    #group_by(kma_ya_sample_name, type, comment) %>% 
    #mutate(rank = row_number(primer)) 




# Make sure to add all columns



# This dummy will make sure that missing columns are added. 
# Should only be a problem when csc doesn't support the variants that SSI defines necessary.
dummy = tibble(nas = c(NA), cols_ = poi) %>% 
    pivot_wider(names_from = cols_, values_from = nas) %>% 
    mutate(kma_ya_sample_name = "dummy") %>% 
    select(kma_ya_sample_name, everything())
    
# add dummy records for keeping all relevant variants 
# dummy = tibble(kma_ya_sample_name = "dummy-sample",
#                type = "dummy",
#                position = poi)



# Trying a new approach. Long-pivoting the joined-table immediately:
twin_called = joined %>% 
    select(kma_ya_sample_name, plate, type, primer, starts_with("Sseq_")) %>%
  
    # Add dummy data to ensure that all poi columns will be present
    bind_rows(tibble(kma_ya_sample_name = "dummy", type = "dummy", primer = "left_primer")) %>% 
    left_join(dummy, by = "kma_ya_sample_name", suffix = c("", "_dummy")) %>% 
    select(-ends_with("_dummy")) %>% 

    pivot_longer(starts_with("Sseq_"), names_to = "position", values_to = "call") %>% 
  

    
  
    # Recode the calls 
    mutate(call = recode(call, "0" = "negativ", "1" = "positiv")) %>% # hold ikke-kaldte som NA, så man senere kan coalesce den første kaldte værdi.
    
    
    # Have a look
    arrange(kma_ya_sample_name, position) %>% 
  
    # pivot left and right out to two columns
    #mutate(both = paste(position, primer)) %>% 
    pivot_wider(id_cols = c(kma_ya_sample_name, plate, position, type), names_from = primer, values_from = call) %>% 
     
    # Have a look
    #View
  
    # Now call the variants using both directions
    mutate(twin_call = case_when(`left_primer` == `right_primer` ~ `left_primer`,
                                 is.na(`left_primer`) | is.na(`right_primer`) ~ coalesce(`left_primer`, `right_primer`, "inkonklusiv_begge"),
                                 TRUE ~ "inkonklusiv_forskellig"))
  
    # 
    # 
    # group_by(kma_ya_sample_name, name) %>% 
    # summarize(paste(value, collapse = ", ")) %>% 
    # 
    # View





# TODO: Warn about non-passing negative controls






# Warn about differing cals between the two directions
inc_differing = twin_called %>% filter(twin_call == "inkonklusiv_forskellig")
if (nrow(inc_differing) > 0) {
  write(paste("Warning: discrepancies beween the two directions in the following samples:\n",
              paste(inc_differing %>% pull(kma_ya_sample_name), collapse = "\t")), stderr())
}



# 
# definition = tribble(
#   ~Sseq_variantbeskrivelse, ~Sseq_smitsomhed, ~Sseq_E848K, ~



out = twin_called %>% 

  
    # Homogenize inconclusives
    mutate(twin_call = if_else(str_detect(twin_call, "inkonklusiv"), "inkonklusiv", twin_call)) %>%
  
    # Pivot back to wide to combine multiple positions per sample
    pivot_wider(id_cols = c(kma_ya_sample_name, plate, type), names_from = position, values_from = twin_call) %>% 
    filter(type != "dummy") %>% 
  
    # Pick out the positions of interest
    select(`sample-id` = kma_ya_sample_name, type, all_of(poi)) %>% 
  
    
    

    # Now add the overall calls
    mutate(Sseq = case_when(Sseq_E484K == "negativ" & Sseq_N501Y == "positiv" & Sseq_Q677H == "negativ" & Sseq_P681H == "positiv" ~ "Spike: Forenelig med B.1.1.7|Variant med øget smitsomhed",
                            Sseq_E484K == "positiv" & Sseq_N501Y == "positiv" & Sseq_Q677H == "negativ" & Sseq_P681H == "negativ" ~ "Spike: Forenelig med B.1.351 eller P.1|Variant med øget smitsomhed og nedsat følsomhed for antistoffer",
                            Sseq_E484K == "positiv" & Sseq_N501Y == "negativ" & Sseq_Q677H == "positiv" & Sseq_P681H == "negativ" ~ "Spike: Forenelig med B.1.525|Variant med nedsat følsomhed for antistoffer",
                            Sseq_E484K == "positiv" & Sseq_N501Y == "negativ" &                           Sseq_P681H == "negativ" ~ "Spike: E484K mutation|Variant med nedsat følsomhed for antistoffer",
                            Sseq_E484K == "negativ" & Sseq_N501Y == "positiv" &                           Sseq_P681H == "negativ" ~ "Spike: N501Y mutation|Variant med øget smitsomhed",
                            Sseq_E484K == "negativ" & Sseq_N501Y == "negativ" & Sseq_Q677H == "negativ" & Sseq_P681H == "negativ" ~ "Spike: Forenelig med oprindelig variant|Variant med formodet normal smitsomhed",
                            TRUE                                                                                                  ~ "Spike: inkonklusiv|Prøven er ikke sekventerbar")) %>% 
    separate(Sseq, c("Sseq_variantbeskrivelse", "Sseq_smitsomhed"), sep = "\\|") %>%
  
    
    


    mutate(MDSU = "32092")
  
    # Have a look before pivoting
    
      

#out %>% format_tsv


out %>% 
  
  pivot_longer(c(Sseq_variantbeskrivelse, Sseq_smitsomhed, starts_with("Sseq_"))) %>% 
  
  filter(type == "sample") %>% 
  select(-type) %>% 

  #arrange(`sample-id`, name) %>% # DO NOT SORT
  write.table(paste0("~/GenomeDK/clinmicrocore/pipe19/batch/mads/output/32092_Sseq_", arg_batch, ".csv"), quote = F, sep = ";", fileEncoding = "cp1252", row.names = F)   # TODO: set output path for args


    





identity









































    
    
    
    