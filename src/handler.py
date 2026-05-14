"""TWOG-owned RunPod worker for smoke-scale MD/ligand-prep validation.

This worker intentionally starts with a conservative scope: validate the input
contract, sanitize the protein PDB, generate a 3D ligand from SMILES, and prepare
ligand PDBQT from SDF/MOL data rather than from ligand PDB. Docking and true MD
are explicit later stages and are reported as skipped unless enabled by worker
configuration.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import tempfile
import traceback
from typing import Any


WORKER_VERSION = "twog-md-smoke-v1"
TAIL_LIMIT = 4000


@dataclass
class StageFailure(Exception):
    stage: str
    message: str
    details: dict[str, Any] = field(default_factory=dict)

    def __str__(self) -> str:
        return self.message


def handler(job: dict[str, Any]) -> dict[str, Any]:
    """RunPod handler entrypoint."""

    payload = job.get("input") if isinstance(job, dict) and isinstance(job.get("input"), dict) else job
    if not isinstance(payload, dict):
        payload = {}
    return run_worker(payload)


def run_worker(payload: dict[str, Any]) -> dict[str, Any]:
    result = _base_result(payload)
    try:
        input_payload = _validate_input(payload)
        result["stages"].append(_stage("input_validation", "completed", required_fields=sorted(_required_fields())))

        with tempfile.TemporaryDirectory(prefix="twog_md_") as tmp:
            workdir = Path(tmp)
            protein_path = _sanitize_protein_pdb(input_payload["protein_pdb"], workdir)
            result["stages"].append(
                _stage(
                    "protein_prep",
                    "completed",
                    artifact="protein.pdb",
                    line_count=_line_count(protein_path),
                    byte_count=protein_path.stat().st_size,
                )
            )
            result["artifacts"]["protein_pdb"] = _artifact_summary(protein_path)

            ligand_info = _prepare_ligand_3d(input_payload["compound_smiles"], input_payload["compound_name"], workdir)
            result["stages"].append(_stage("ligand_3d", "completed", **ligand_info["stage_details"]))
            result["artifacts"]["ligand_sdf"] = _artifact_summary(ligand_info["sdf_path"])
            if ligand_info.get("mol_path"):
                result["artifacts"]["ligand_mol"] = _artifact_summary(ligand_info["mol_path"])

            ligand_pdbqt = _prepare_ligand_pdbqt(ligand_info["sdf_path"], workdir)
            result["stages"].append(_stage("ligand_pdbqt", "completed", **ligand_pdbqt["stage_details"]))
            result["artifacts"]["ligand_pdbqt"] = _artifact_summary(ligand_pdbqt["pdbqt_path"])

            result["stages"].append(_docking_stage(input_payload))
            result["stages"].append(_md_smoke_stage(input_payload))
    except StageFailure as exc:
        result["status"] = "failed"
        result["errors"].append(
            {
                "stage": exc.stage,
                "message": exc.message,
                **exc.details,
            }
        )
        result["stages"].append(_stage(exc.stage, "failed", message=exc.message, **exc.details))
    except Exception as exc:  # pragma: no cover - defensive RunPod diagnostic path.
        result["status"] = "failed"
        result["errors"].append(
            {
                "stage": "unhandled_exception",
                "message": str(exc),
                "traceback": traceback.format_exc()[-TAIL_LIMIT:],
            }
        )
        result["stages"].append(_stage("unhandled_exception", "failed", message=str(exc)))
    return result


def _base_result(payload: dict[str, Any]) -> dict[str, Any]:
    return {
        "status": "completed",
        "worker_version": WORKER_VERSION,
        "input_summary": {
            "compute_job_id": payload.get("compute_job_id"),
            "queue_item_id": payload.get("queue_item_id"),
            "target_name": payload.get("target_name"),
            "compound_name": payload.get("compound_name"),
            "simulation_steps": payload.get("simulation_steps"),
            "temperature": payload.get("temperature"),
            "ph": payload.get("ph"),
            "force_field": payload.get("force_field"),
            "solvent_model": payload.get("solvent_model"),
        },
        "stages": [],
        "artifacts": {},
        "warnings": [],
        "errors": [],
    }


def _required_fields() -> set[str]:
    return {
        "protein_pdb",
        "compound_smiles",
        "target_name",
        "compound_name",
        "simulation_steps",
        "temperature",
        "protein_source",
        "ligand_source",
        "preparation_method",
    }


def _validate_input(payload: dict[str, Any]) -> dict[str, Any]:
    missing = [field for field in sorted(_required_fields()) if not str(payload.get(field) or "").strip()]
    if missing:
        raise StageFailure("input_validation", "Missing required input fields.", {"missing_fields": missing})
    protein_pdb = str(payload["protein_pdb"])
    if "ATOM" not in protein_pdb and "HETATM" not in protein_pdb:
        raise StageFailure("input_validation", "protein_pdb must contain at least one ATOM or HETATM record.")
    if "TER" not in protein_pdb and "END" not in protein_pdb:
        raise StageFailure("input_validation", "protein_pdb must contain TER or END.")
    try:
        simulation_steps = int(payload.get("simulation_steps") or 0)
    except (TypeError, ValueError) as exc:
        raise StageFailure("input_validation", "simulation_steps must be an integer.") from exc
    if simulation_steps < 1 or simulation_steps > 1000:
        raise StageFailure(
            "input_validation",
            "simulation_steps must be between 1 and 1000 for smoke runs.",
            {"simulation_steps": simulation_steps},
        )
    return dict(payload)


def _sanitize_protein_pdb(protein_pdb: str, workdir: Path) -> Path:
    retained: list[str] = []
    for raw_line in protein_pdb.splitlines():
        line = raw_line.rstrip()
        record = line[:6].strip().upper()
        if record in {"ATOM", "TER", "END"}:
            retained.append(line)
    if not any(line.startswith("ATOM") for line in retained):
        raise StageFailure("protein_prep", "No protein ATOM records were retained after sanitization.")
    if not retained[-1].startswith("END"):
        retained.append("END")
    path = workdir / "protein.pdb"
    path.write_text("\n".join(retained) + "\n", encoding="utf-8")
    return path


def _prepare_ligand_3d(smiles: str, compound_name: str, workdir: Path) -> dict[str, Any]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as exc:  # pragma: no cover - covered in container, not local repo env.
        raise StageFailure(
            "ligand_3d",
            "RDKit is not available in the worker environment.",
            {"exception": repr(exc)},
        ) from exc

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise StageFailure("ligand_3d", "compound_smiles could not be parsed by RDKit.", {"compound_smiles": smiles})
    mol = Chem.AddHs(mol)
    embed_code = AllChem.EmbedMolecule(mol, randomSeed=17)
    if embed_code != 0:
        raise StageFailure("ligand_3d", "RDKit failed to embed a 3D conformer.", {"embed_code": embed_code})
    method = "MMFF94"
    if AllChem.MMFFHasAllMoleculeParams(mol):
        optimize_code = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    else:
        method = "UFF"
        optimize_code = AllChem.UFFOptimizeMolecule(mol, maxIters=500)

    mol.SetProp("_Name", compound_name)
    sdf_path = workdir / "ligand.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()
    mol_path = workdir / "ligand.mol"
    mol_block = Chem.MolToMolBlock(mol)
    mol_path.write_text(mol_block, encoding="utf-8")
    return {
        "sdf_path": sdf_path,
        "mol_path": mol_path,
        "stage_details": {
            "compound_name": compound_name,
            "atom_count": int(mol.GetNumAtoms()),
            "conformer_count": int(mol.GetNumConformers()),
            "optimization_method": method,
            "optimization_code": int(optimize_code),
            "intermediate_format": "sdf",
        },
    }


def _prepare_ligand_pdbqt(sdf_path: Path, workdir: Path) -> dict[str, Any]:
    command = _find_command("mk_prepare_ligand.py")
    if command is None:
        raise StageFailure("ligand_pdbqt", "mk_prepare_ligand.py is not available in the worker image.")
    pdbqt_path = workdir / "ligand.pdbqt"
    completed = _run_subprocess([command, "-i", str(sdf_path), "-o", str(pdbqt_path)], stage="ligand_pdbqt")
    if not pdbqt_path.exists() or pdbqt_path.stat().st_size == 0:
        raise StageFailure(
            "ligand_pdbqt",
            "Ligand PDBQT preparation completed without producing a non-empty artifact.",
            {"command": completed["command"], "stdout_tail": completed["stdout_tail"], "stderr_tail": completed["stderr_tail"]},
        )
    return {
        "pdbqt_path": pdbqt_path,
        "stage_details": {
            "command": completed["command"],
            "return_code": completed["return_code"],
            "stdout_tail": completed["stdout_tail"],
            "stderr_tail": completed["stderr_tail"],
            "input_format": "sdf",
            "output_format": "pdbqt",
        },
    }


def _docking_stage(payload: dict[str, Any]) -> dict[str, Any]:
    if not _truthy(payload.get("enable_docking")):
        return _stage(
            "docking",
            "skipped",
            reason="Docking is disabled in smoke-v1 unless enable_docking=true is supplied.",
        )
    vina = _find_command("vina")
    if vina is None:
        return _stage("docking", "failed", message="vina executable is not available in the worker image.")
    return _stage("docking", "skipped", reason="Docking command wiring is intentionally deferred until prep smoke passes.")


def _md_smoke_stage(payload: dict[str, Any]) -> dict[str, Any]:
    return _stage(
        "md_smoke",
        "skipped",
        reason="OpenMM MD execution is intentionally deferred until ligand/receptor preparation passes.",
        simulation_steps=int(payload.get("simulation_steps") or 0),
    )


def _run_subprocess(command: list[str], *, stage: str) -> dict[str, Any]:
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    result = {
        "command": command,
        "return_code": completed.returncode,
        "stdout_tail": (completed.stdout or "")[-TAIL_LIMIT:],
        "stderr_tail": (completed.stderr or "")[-TAIL_LIMIT:],
    }
    if completed.returncode != 0:
        raise StageFailure(stage, "Subprocess failed.", result)
    return result


def _find_command(name: str) -> str | None:
    return shutil.which(name)


def _stage(stage: str, status: str, **details: Any) -> dict[str, Any]:
    return {"stage": stage, "status": status, **{key: value for key, value in details.items() if value is not None}}


def _artifact_summary(path: Path) -> dict[str, Any]:
    return {
        "path": str(path),
        "file_name": path.name,
        "byte_count": path.stat().st_size if path.exists() else 0,
        "line_count": _line_count(path) if path.exists() else 0,
    }


def _line_count(path: Path) -> int:
    return len(path.read_text(encoding="utf-8", errors="replace").splitlines())


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y", "on"}
    return bool(value)


def main() -> None:
    if len(sys.argv) > 1:
        job = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
        print(json.dumps(handler(job), indent=2, sort_keys=True))
        return
    try:
        import runpod
    except Exception as exc:
        raise SystemExit(f"runpod package is required when no local test input is supplied: {exc}") from exc
    runpod.serverless.start({"handler": handler})


if __name__ == "__main__":
    main()

