# staphit2-python: pinned Python + libraries for staphit2's bin/ scripts.
# Built and pushed by .github/workflows/staphit2-python-image.yml; tag = VERSION.
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
USER root
LABEL org.opencontainers.image.source="https://github.com/alarawms/staphit2" \
      org.opencontainers.image.description="Python environment for staphit2 bin/ scripts"
