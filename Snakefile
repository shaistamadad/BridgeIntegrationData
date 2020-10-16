datasets = ["pbmc_3k", "pbmc_10k", "pbmc_humanbrain_3k"]

rule download:
    input: "datasets/{dset}.txt"
    output: touch("data/{dset}/download.done")
    message: "Download datasets"
    shell:
        """
        wget -i {input} -P data/{wildcards.dset}
        """

rule process:
    input: "data/{dset}/download.done"
    output: "objects/{dset}.rds"
    message: "Process data"
    shell: "Rscript code/dataset_processing/process_{wildcards.dset}.R"

