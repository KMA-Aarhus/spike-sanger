###############################################################################
# This script rummages through all the .xls files, and join the content to a 
# eurofins sanger run.

library(tidyverse)


args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste(1:length(args), args), stderr())
write("", stderr())

# Read cli-args
arg_batch = args[1]
csc_results_file = args[2]
final_out_file = args[3]
xls_path_ = "input/total"


#metadata_file = args[1] # deprecated as it pools all of them.


devel = F
# devel = T


if (devel) {
    rm(list = ls())
    devel = T
    xls_path_ = "~/GenomeDK/clinmicrocore/spike-sanger/input/total"
    
    
    #arg_batch = "210214"
    #arg_batch = "210216"
    #arg_batch = "11107265503-1"
    
    #metadata_file = paste0("~/GenomeDK/clinmicrocore/spike-sanger/input/", arg_batch, ".xls") # not used anymore as the files are now pooled together.
    #csc_results_file = paste0("~/GenomeDK/clinmicrocore/spike-sanger/output/", arg_batch, "/csc/", arg_batch, "_results.csv")
    #final_out_file = paste0("~/GenomeDK/clinmicrocore/spike-sanger/mads_out/32092_Sseq_", arg_batch, ".csv")
    
    
    # taken directly from the snakefile:
    #arg_batch = "11107289500-1"
    #csc_results_file = "~/GenomeDK/clinmicrocore/spike-sanger/output/11107289500-1/csc/11107289500-1_results.csv"
    #final_out_file = "~/GenomeDK/clinmicrocore/spike-sanger/mads_out/32092_Sseq_11107289500-1.csv"
    
    
    # Another batch:
    # arg_batch = "11107267912-1"
    # csc_results_file = "~/GenomeDK/clinmicrocore/spike-sanger/output/11107267912-1/csc/11107267912-1_results.csv"
    # final_out_file = "~/GenomeDK/clinmicrocore/spike-sanger/mads_out/32092_Sseq_11107267912-1.csv"
    
    # Yet another
    arg_batch = "11107294571-1"
    csc_results_file = "~/GenomeDK/clinmicrocore/spike-sanger/output/11107294571-1/csc/11107294571-1_results.csv"
    final_out_file = "~/GenomeDK/clinmicrocore/spike-sanger/mads_out/32092_Sseq_11107294571-1.csv"
    
    
    
}









## Before we do anything, we want to concatenate all the metadata-files. Other

collect_all_metadata = function(xls_path) {
  write("reading dirs:", stdout())
  xls_in = dir(xls_path, full.names = T)
  write(paste("xls_in:", xls_in), stdout())
  
  df = tibble()
  for (i in xls_in) {
    
    write(paste0("Reading:\n  ", i), stderr())
    
    temp = readxl::read_excel(i, col_types = "text")
    
    write(paste0("Names:\n  ", paste(names(temp), collapse = ", "), stderr()))
    
    #write("name", stderr())
    #write(paste0(glimpse(df)), stderr())
    temp = temp %>% 
      select(ef_sample_name = `Stregkode nr`, kma_raw_sample_name = `Originalnr.`, primer_raw = `Primer-bakke`) %>% 
      mutate(table_file = basename(i))
      #drop_na()
    
    # Warn if any of the needed columns are missing.
    # And also if excessive NA values exists.
  
    
    write(paste0("Samples read:\n  ", dim(temp)[1]), stderr())
    
    
    df = temp %>% bind_rows(df)
    
    write("", stdout())
    
  }
  
  df
}

write("Collecting all metadata ...", stderr())
metadata_read = collect_all_metadata(xls_path_) 

#metadata %>% select(






# Below is the old metadata ingestion:
#metadata = read_tsv(metadata_file)
#metadata_read = readxl::read_excel(metadata_file, sheet = 1, na = c("", "?"))

metadata_read %>% glimpse
write("Transforming metadata ...", stderr())
metadata = metadata_read %>% 
    #select(kma_raw_sample_name = `Originalnr.`, primer_raw = Primer, ef_sample_name = `Rør nr`) %>% 
    #filter(!is.na(kma_raw_sample_name)) %>% 
    drop_na(ef_sample_name) %>% 
    mutate(primer = case_when(str_detect(primer_raw, "^L-") ~ "left_primer",
                              str_detect(primer_raw, "^R-") ~ "right_primer"),
           plate = str_extract(primer_raw, "\\d+$")) %>%

    # Mark the sample type
    # Should be synced with pipe19/scripts/parse_path.r:
    mutate(type = case_when(str_detect(tolower(kma_raw_sample_name), "positiv|pos$|seqpos") ~ "positive_control",
                            str_detect(tolower(kma_raw_sample_name), "negativ|h2o|^empty|blank|tom|^neg|neg$") ~ "negative_control",
                            str_detect(tolower(kma_raw_sample_name), "^afd|^00") ~ "other",
                            TRUE ~ "sample")) %>% 
    
    # Consider warning about the samples that are missed here:
    drop_na(kma_raw_sample_name) %>% 
    
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
    
    select(kma_ya_sample_name, kma_raw_sample_name, ef_sample_name, plate, primer, type, table_file)
# Old metadata ingestion done.


    



# Define Positions Of Interest.
# Please make sure all of these start with "Sseq_"
poi = c("Sseq_E484K", "Sseq_N501Y", "Sseq_Q677H", "Sseq_P681H")


# This dummy will make sure that missing columns are added. 
# Should only be a problem when csc doesn't support the variants that SSI defines as necessary.
dummy = tibble(nas = c(NA), cols_ = poi) %>% 
  pivot_wider(names_from = cols_, values_from = nas) %>% 
  mutate(ef_sample_name = "dummy") %>% 
  select(ef_sample_name, everything())



# Read csc-results
write("reading csc-results ...", stderr())
results_read = read_csv(csc_results_file, col_types = cols(.default = "c")) 
results_read %>% glimpse
results = results_read %>% 
    mutate(ef_sample_name = str_remove(sample, "_.+$")) %>% # Clean the sample name.
    select(ef_sample_name, ef_raw_sample_name = sample, comment, everything()) %>% 
  
  
    


    # Prefix variant-position-columns with "Sseq_"
    # Be aware that this call is dependent on the number, and ordering of columns
    # The forth column should be the leftmost variant column.
    rename_at(.vars = vars(4:last_col()),
              .funs = ~ paste0("Sseq_", .x)) %>% 
  
   
    # add the dummy, and delete it again. 
    # seems stupid, but remember that it adds NA-columns which will be very useful later!
    bind_rows(dummy) %>% 
    filter(ef_sample_name != "dummy") %>% 
  
    # To make matters simpler, we're only interested in keeping the poi columns.
    select(ef_sample_name, comment, all_of(poi))



# Join the data together
write("joining meta and results ...")
joined = results %>% 
    left_join(metadata, by = "ef_sample_name") %>% 
    arrange(kma_ya_sample_name) %>% 
  
  
  
  
    select(ef_sample_name, kma_ya_sample_name, comment, plate, primer, type, all_of(poi))
  
  
    
    # Quick and dirty way of removing duplicated names:
    #distinct(kma_ya_sample_name, plate, primer, .keep_all = T)
    # I have chosen to a use a different quick and dirty method further down: use the "first()" summary function for picking unique values in pivot wider for twin_calling

    # # Set a rank that shows duplicated names
    # group_by(kma_ya_sample_name) %>%
    # mutate(dup_rank = row_number())
    # 

  
  

# This block is deprecated since duplicates are removed catastrophically.
# Now where we have joined the the results and metadata, we can check that there is ONLY TWO samples per sample name.
# I have had some problems with duplicated names for negative controls
# dup_names = joined %>% 
#   count(kma_ya_sample_name) %>% 
#   filter(n > 2) %>% 
#   pull(kma_ya_sample_name)
# write(paste("Warning: duplicated names for:", dup_names), stderr())
#   
    



    
    # Mark the overlapping samples
    #group_by(kma_ya_sample_name, type, comment) %>% 
    #mutate(rank = row_number(primer)) 




# Make sure to add all columns




    
# add dummy records for keeping all relevant variants 
# dummy = tibble(kma_ya_sample_name = "dummy-sample",
#                type = "dummy",
#                position = poi)

#joined %>% write_tsv("joined_before_twincalling.tsv")


# Trying a new approach. Long-pivoting the joined-table immediately:
#Kunne man ikke på en eller anden måde, få twin_called til at indeholde de to ef_sample_names der indgår i mads-prøven?
write("twin calling ...", stderr())
twin_called = joined %>% 
    #select(kma_ya_sample_name, plate, type, primer, starts_with("Sseq_")) %>%
  
  
    # Because the dummy is added before joining the tables, we now have to come in and add some information.
    # We need to ensure the presence of both primers, otherwise the later pivot will fail.
    bind_rows(tibble(kma_ya_sample_name = "dummy2", primer = c("left_primer", "right_primer"))) %>% 
  
    # Alternatively to "starts_with()", one could use "all_of()"
    pivot_longer(starts_with("Sseq_"), names_to = "position", values_to = "call") %>% 
  
    

  

    
  
    # Recode the calls before splitting the directions out into distinct columns.
    mutate(call = recode(call, "0" = "negativ", "1" = "positiv")) %>% # hold ikke-kaldte som NA, så man senere kan coalesce den første kaldte værdi.

    #write_tsv("temp3_before_pivot_just_called.tsv")    
  
    
  
    # We do not need the dummy any more.
    #filter(ef_sample_name != "dummy") %>%   
    # I moved this up earlier. Hence the commenting.
    
  
    # Have a look
    #arrange(kma_ya_sample_name, position) %>% View

  
    # pivot left and right out to two columns
    #mutate(both = paste(position, primer)) %>% 
    arrange(is.na(call)) %>% 
    pivot_wider(id_cols = c(kma_ya_sample_name, plate, position, type), names_from = primer, values_from = call, values_fn = first) %>% 
    
    
    #write_tsv("temp4_just_pivoted.tsv")
      
    # Immediately after this pivot, we can remove the second dummy
    filter(kma_ya_sample_name != "dummy2") %>% 
    
    
  
    # Now call the variants using both directions
    mutate(twin_call = case_when(`left_primer` == `right_primer` ~ `left_primer`,     # Hvis de er enige, så skriv det de er enige om
                                 is.na(`left_primer`) | is.na(`right_primer`) ~ coalesce(`left_primer`, `right_primer`, "inkonklusiv_begge"), # Hvis en er NA, så skriv den anden. Når begge er NA skal der stå "inkonklusiv begge".
                                 TRUE ~ "inkonklusiv_forskellig")) %>%                # Ellers, ved uenighed, skal der stå inkonklusiv forskellig".
    
  
    #drop_na()  
    identity
  
    # 
    # 
    # group_by(kma_ya_sample_name, name) %>% 
    # summarize(paste(value, collapse = ", ")) %>% 
    # 
    # View





# I'm saving this table to disk as well, because it is nice to have something to debug with.
#twin_called %>% write_tsv(paste0(dirname(final_out_file), "/debug1_Sseq_", arg_batch, ".tsv"))




# TODO: Warn about non-passing negative controls


# Warn about differing cals between the two directions
inc_differing = twin_called %>% filter(twin_call == "inkonklusiv_forskellig")
if (nrow(inc_differing) > 0) {
  write(paste("Warning: discrepancies beween the two directions in the following samples:\n",
              paste(inc_differing %>% pull(kma_ya_sample_name), collapse = "\t")), stderr())
} else {
  write(paste("Info: No disrepancies when comparing the two directions between the samples."))
}



write("applying conditionals ...", stderr())
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
    
      

# I'm saving this table to disk as well, because it is nice to have something to debug with.
out %>% write_tsv(paste0(dirname(final_out_file), "/debug2_Sseq_", arg_batch, ".tsv"))



write("writing file out", stderr())
out %>% 
  
  pivot_longer(c(Sseq_variantbeskrivelse, Sseq_smitsomhed, starts_with("Sseq_"))) %>% 
  
  filter(type == "sample") %>% 
  select(-type) %>% 

  # DO NOT SORT
  write.table(final_out_file, quote = F, sep = ";", fileEncoding = "cp1252", row.names = F)   # TODO: set output path for args
  

    
# Overvej at der måske også skal printes en debug tabel, hvor man kan se hvad der går galt med hvilke prøver..
# Bare så det er lidt nemmere at tjekke, at alt er i orden.
# For kørslen af batch *267912-1 er der fx. ingen right_primer's til stede?!?





    
    