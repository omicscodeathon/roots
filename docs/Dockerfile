FROM python:3.13-slim

# build-essential covers any dependency without a prebuilt wheel for this
# platform/Python combo; harmless to keep even if everything resolves to
# wheels.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# ligand_prep.py's only third-party imports: pandas, rdkit, meeko     scipy     gemmi
# (subprocess/sys/os/csv/time/json/multiprocessing/pathlib/datetime/
# concurrent.futures are all stdlib, nothing to install for those).
RUN pip install --no-cache-dir \
    pandas \
    rdkit \
    meeko     scipy     gemmi

# step3/, step4/, and ligand_prep.py itself are bind-mounted at runtime
# (see the -v flags in run_with_notify.sh's docker run call), so nothing
# is COPY'd into the image -- this keeps rebuilds fast since editing the
# script never invalidates this layer.

CMD ["python3", "/workspace/ligand_prep.py"]
