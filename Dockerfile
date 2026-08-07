FROM mambaorg/micromamba:latest

LABEL org.opencontainers.image.source="https://github.com/metagenlab/zAMPExplorer"
LABEL org.opencontainers.image.description="zAMPExplorer Shiny app for reproducible 16S microbiome downstream analysis"

COPY --chown=$MAMBA_USER:$MAMBA_USER . /pkg

RUN micromamba config set extract_threads 1 && \
    micromamba install -n base -y -f /pkg/env.yml && \
    micromamba clean -afy

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN Rscript /pkg/install_dependencies.R && R CMD INSTALL /pkg

EXPOSE 3838
WORKDIR /results

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "R", "-e", "options(shiny.host='0.0.0.0', shiny.port=3838, browser=FALSE); library(zAMPExplorer); zAMPExplorer::zAMPExplorer_app()"]
