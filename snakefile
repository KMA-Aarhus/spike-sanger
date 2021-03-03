
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
                "mads_out/32092_Sseq_{batch}.csv"],
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
    shell:
        """
        
        mkdir -p {output.ab1_dir}
        cp {params.cp_in} {output.ab1_dir}


        singularity run docker://kblin/covid-spike-classification \
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
        "mads_out/32092_Sseq_{batch}.csv"
    shell:
        """

        singularity run docker://rocker/tidyverse \
            Rscript scripts/output_mads_sseq.r {wildcards.batch} {input} {output}


        """


