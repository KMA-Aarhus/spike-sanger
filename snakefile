
# snakemake --profile configs/slurm all

print("<!-- ")



import os
from datetime import datetime
import time
import re
from shutil import copyfile

import pandas as pd
import re

configfile: "config.json"

input_table = config["input_table"]
out_base = config["out_base"]


df = pd.read_table(input_table, sep = "\t", dtype = str, comment = "#")

print(df)
print("//")
print()

print(" -->")





onerror:
    print("An error occurred")
    shell("mail -s 'spike-sanger pipeline error' kobel@pm.me < {log}")




rule all:
    input:
        expand(["output/{batch}/csc/{batch}_results.csv",
                "output/{batch}/32092_Sseq_{batch}.csv",
                "output/{batch}/32092_Sseq_{batch}_mail_sent.flag"],
                batch = df["batch"])




rule start:
    input:
        input_table
    output:
        ab1_dir = directory("output/{batch}/ab1"),
        final = "output/{batch}/csc/{batch}_results.csv"
    params:
        reference = config["reference"],
        cp_in = lambda wildcards: df[df["batch"] == wildcards.batch]["raw_path"].values[0]
    container: "docker://kblin/covid-spike-classification"
    shell:
        """
        
        mkdir -p {output.ab1_dir}
        cp {params.cp_in} {output.ab1_dir}


        covid-spike-classification \
            --reference {params.reference} \
            -i ab1 \
            --outdir output/{wildcards.batch}/csc \
            {output.ab1_dir}

        mv output/{wildcards.batch}/csc/results.csv {output.final}


        """



# Now the data is generated, and we just need to generate some reports for mads.
rule mads_out:
    input:
        "output/{batch}/csc/{batch}_results.csv"
    output:
        "output/{batch}/32092_Sseq_{batch}.csv"
    container: "docker://rocker/tidyverse"
    shell:
        """

    
        Rscript scripts/output_mads_sseq.r {wildcards.batch} {input} {output}


        """

rule send_mail:
    input: "output/{batch}/32092_Sseq_{batch}.csv"
    output: "output/{batch}/32092_Sseq_{batch}_mail_sent.flag"
    shell:
        """
        
        
        mail -v -s "Automail: mads Sseq-svar" -a {input} carkob@rm.dk <<< "Autogenereret mads svar for Eurofins sanger-sekventeringsbatch: {wildcards.batch}\nMutationssignaturerne er senest opdateret d. 16. feb. 2021.\n\nGenereret med $(pwd)/snakefile."

        touch {output}

        """


