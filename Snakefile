datasets = ["pbmc_atac", "pbmc_rna", "adult_brain_ss1",
            "adult_brain_ss2", "allen_rna", "adult_brain_sci",
            "human_brain_rna", "human_brain_atac", "pbmc_multiomic"]

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
