# Philomenos container: the `philomenos` command plus every external tool, each
# in its own micromamba env created from envs/*.yml (philomenos, snippy,
# gubbins, pairsnp, quast, qc, mask; envs/qc.yml holds skani and CheckM2).
#
#   docker build -t philomenos .
#   docker run --rm philomenos smoke --tools-from /opt/philomenos/tools.container.yaml --outdir /tmp/smoke
#
# Real runs: mount a working directory and run as your own user so outputs are
# yours. Do not use `run --detach` inside a container (the container stops when
# its main process exits); use `docker run -d` instead:
#
#   docker run -d --name philo --user "$(id -u):$(id -g)" -v "$PWD":/work philomenos \
#       run -c config.yaml
#   docker exec philo philomenos status runs/<name>
#
# Container config: tools.envs take env NAMES (not paths). The file
# /opt/philomenos/tools.container.yaml, written below, holds:
#
#   tools:
#     conda: conda          # a shim that forwards `conda run` to micromamba
#     envs:
#       snippy: snippy
#       gubbins: gubbins
#       iqtree: gubbins
#       pairsnp: pairsnp
#       quast: quast
#       skani: qc
#       checkm2: qc
#       mask: mask
#
# Copy that `tools:` block into your run config. CheckM2 (qc.checkm2: true)
# also needs its database, downloaded once into a mounted folder:
#   docker run --rm -v "$PWD":/work --entrypoint micromamba philomenos \
#       run -n qc checkm2 database --download --path /work/checkm2_db

FROM mambaorg/micromamba:1.5.10

ENV MAMBA_ROOT_PREFIX=/opt/conda \
    MPLCONFIGDIR=/tmp/matplotlib \
    XDG_CACHE_HOME=/tmp/.cache \
    PYTHONDONTWRITEBYTECODE=1

# One env per tool group, so their dependency pins never conflict.
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/ /tmp/envs/
RUN set -eux; \
    for env in philomenos snippy gubbins pairsnp quast qc mask; do \
        micromamba create -y -f "/tmp/envs/${env}.yml"; \
    done; \
    micromamba clean --all --yes; \
    rm -rf /tmp/envs

# The package itself, into the philomenos env.
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md /opt/philomenos/src/
COPY --chown=$MAMBA_USER:$MAMBA_USER src/ /opt/philomenos/src/src/
RUN micromamba run -n philomenos pip install --no-deps --no-cache-dir /opt/philomenos/src

USER root
# `conda` shim: Philomenos calls `conda run -n <env> --no-capture-output <cmd>`;
# micromamba streams output already and does not know that flag.
RUN set -eux; \
    printf '%s\n' \
        '#!/usr/bin/env bash' \
        '# Forward conda calls to micromamba, dropping --no-capture-output.' \
        'args=()' \
        'for a in "$@"; do [[ "$a" == "--no-capture-output" ]] || args+=("$a"); done' \
        'exec micromamba "${args[@]}"' \
        > /usr/local/bin/conda; \
    chmod 0755 /usr/local/bin/conda; \
    printf '%s\n' \
        '# Tool locations inside the philomenos container (env names, not paths).' \
        'tools:' \
        '  conda: conda' \
        '  envs:' \
        '    snippy: snippy' \
        '    gubbins: gubbins' \
        '    iqtree: gubbins' \
        '    pairsnp: pairsnp' \
        '    quast: quast' \
        '    skani: qc' \
        '    checkm2: qc' \
        '    mask: mask' \
        > /opt/philomenos/tools.container.yaml; \
    mkdir -p /work && chmod 0777 /work
USER $MAMBA_USER

WORKDIR /work
ENTRYPOINT ["/opt/conda/envs/philomenos/bin/philomenos"]
CMD ["--help"]
