# TWOG MD Smoke RunPod Worker

This worker owns the first TWOG MD compute lane. It is not a therapeutic efficacy engine. It proves the worker contract: protein PDB input, ligand SMILES input, deterministic ligand 3D generation, ligand PDBQT preparation from SDF, structured stage diagnostics, and safe smoke-scale behavior.

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
```

## Hosted Endpoint

The intended image is:

```text
ghcr.io/chasepenelli/twog-md-worker:smoke-v1
```

Create a new RunPod serverless endpoint named `twog-md-smoke-v1`, set `workersMin=0`, and set `workersMax=2` for smoke validation. Update `HSA_RUNPOD_ENDPOINT_ID` in GitHub Actions and Dagster+ only after the endpoint is created.

