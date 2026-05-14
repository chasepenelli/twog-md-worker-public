# TWOG MD Smoke RunPod Worker

This worker owns the first TWOG MD compute lane. It is not a therapeutic efficacy engine. It proves the worker contract: protein PDB input, ligand SMILES input, deterministic ligand 3D generation, ligand PDBQT preparation from SDF, optional smoke docking, structured stage diagnostics, and safe smoke-scale behavior.

## Contract

RunPod sends a job shaped as:

```json
{
  "input": {
    "protein_pdb": "ATOM ...",
    "compound_smiles": "CCO",
    "target_name": "target",
    "compound_name": "compound",
    "simulation_steps": 10,
    "temperature": 300.0,
    "enable_docking": false,
    "protein_source": "provenance",
    "ligand_source": "provenance",
    "preparation_method": "provenance"
  }
}
```

The response always includes:

- `status`
- `worker_version`
- `stages`
- `artifacts`
- `warnings`
- `errors`

Structured worker failures return `status="failed"` inside the output payload with the failing stage and diagnostic fields.

`enable_docking=true` is a separate tier. It requires a fresh TWOG expert-review packet and approval because the packet hash changes. When enabled, the worker prepares a receptor PDBQT, runs a small AutoDock Vina smoke dock, and returns the receptor/docked ligand artifacts plus the best parsed affinity when Vina reports one. If docking is enabled but the worker image lacks the docking executable, receptor preparation, or a non-empty docking artifact, the worker returns a structured `docking` failure rather than silently passing the job.

Optional docking controls:

- `docking_center`: object with `x`, `y`, `z`
- `docking_box_size`: number or object with `x`, `y`, `z`
- `docking_exhaustiveness`
- `docking_num_modes`
- `docking_cpu`

## Local Tests

From this directory:

```bash
python src/handler.py test_input_positive_control.json
python src/handler.py test_input_pazopanib_kdr.json
```

The container workflow runs both commands before publishing the image.

## Docker

```bash
docker build --platform linux/amd64 -t twog-md-worker:test .
docker run --rm twog-md-worker:test python /app/src/handler.py /app/test_input_positive_control.json
docker run --rm twog-md-worker:test python /app/src/handler.py /app/test_input_pazopanib_kdr.json
docker run --rm twog-md-worker:test python -c "import shutil; assert shutil.which('vina'); assert shutil.which('mk_prepare_receptor.py')"
```

## Hosted Endpoint

The private repo build publishes:

```text
ghcr.io/chasepenelli/twog-md-worker:smoke-v1
```

RunPod uses the credential-free public mirror image:

```text
ghcr.io/chasepenelli/twog-md-worker-public:smoke-v1
```

Create a new RunPod serverless endpoint named `twog-md-smoke-v1`, set `workersMin=0`, and set `workersMax=2` for smoke validation. Update `HSA_RUNPOD_ENDPOINT_ID` in GitHub Actions and Dagster+ only after the endpoint is created.
