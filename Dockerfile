FROM mambaorg/micromamba:1.5.10

WORKDIR /app
ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.txt /tmp/requirements.txt

RUN micromamba install -y -n base -c conda-forge python=3.11 rdkit scipy gemmi && \
    micromamba clean --all --yes

RUN python -m pip install --no-cache-dir -r /tmp/requirements.txt

COPY --chown=$MAMBA_USER:$MAMBA_USER src /app/src
COPY --chown=$MAMBA_USER:$MAMBA_USER test_input_positive_control.json /app/test_input_positive_control.json
COPY --chown=$MAMBA_USER:$MAMBA_USER test_input_pazopanib_kdr.json /app/test_input_pazopanib_kdr.json
COPY --chown=$MAMBA_USER:$MAMBA_USER README.md /app/README.md

ENV PYTHONPATH=/app/src
ENV PYTHONUNBUFFERED=1

CMD ["python", "-u", "/app/src/handler.py"]
