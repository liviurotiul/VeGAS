FROM mambaorg/micromamba:1.5.8

# Create env with both snakemake and snippy
RUN micromamba create -y -n snakemake-snippy -c conda-forge -c bioconda \
    snakemake snippy && \
    micromamba clean --all --yes

SHELL ["micromamba", "run", "-n", "snakemake-snippy", "/bin/bash", "-c"]

# Optional: show versions
RUN snakemake --version && snippy --version

ENTRYPOINT ["micromamba", "run", "-n", "snakemake-snippy", "snakemake"]
