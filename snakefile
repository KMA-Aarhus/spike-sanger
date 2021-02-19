
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


out_base = config["out_base"]

df = pd.read_table("input_table.tsv", sep = "\t", dtype = str, comment = "#")

print(df)
print("//")
print()

print(" -->")





onerror:
    print("An error occurred")
    shell("mail -s 'sanger pipeline error' kobel@pm.me < {log}")




rule all:
    input:
        expand("output/{batch}/csc/{batch}_results.csv",
                batch = df["batch"])




rule start:
    input:
        "input/{batch}.xls"
    output:
        ab1_dir = directory("output/{batch}/ab1"),
        final = "output/{batch}/csc/{batch}_results.csv"
    params:
        reference = config["reference"],
        cp_in = lambda wildcards: df[df["batch"] == wildcards.batch]["raw_path"].values[0]

                #lambda wildcards: df[df["full_name"]==wildcards.sample][["forward_path", "reverse_path"]].values[0].tolist()
    shell:
        """
        
        mkdir -p {output.ab1_dir}
        cp {params.cp_in} {output.ab1_dir}


        singularity run docker://cmkobel/covid-spike-classification \
            covid-spike-classification \
                --reference {params.reference} \
                -i ab1 \
                --outdir output/{wildcards.batch}/csc \
                {output.ab1_dir}

        mv output/{wildcards.batch}/csc/results.csv {output.final}


        """

