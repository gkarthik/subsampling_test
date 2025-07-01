FROM condaforge/mambaforge:23.11.0-0

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    bash \
    openjdk-17-jre-headless \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow
RUN curl -fsSL https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/nextflow \
    && chmod +x /usr/local/bin/nextflow \
    && nextflow self-update 23.04.0 \
    && rm -rf /root/.nextflow

# Copy and create conda environment
COPY env.yml /tmp/env.yml
RUN mamba env create -f /tmp/env.yml && \
    mamba clean -a -y && \
    rm /tmp/env.yml

# Set environment path
ENV PATH="/opt/conda/envs/subsampling-test/bin:$PATH"

# Activate environment by default
SHELL ["conda", "run", "--no-capture-output", "-n", "subsampling-test", "/bin/bash", "-c"]

WORKDIR /workspace 