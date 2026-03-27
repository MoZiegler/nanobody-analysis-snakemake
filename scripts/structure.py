# from __future__ import annotations

import json
import logging
import random
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

import requests
from Bio import SeqIO
from Bio.PDB import PDBParser


class StructurePredictor(ABC):
    @abstractmethod
    def predict(self, sequence: str) -> dict:
        raise NotImplementedError


class MockPredictor(StructurePredictor):
    def predict(self, sequence: str) -> dict:
        return {
            "backend": "mock",
            "confidence": round(random.uniform(0.65, 0.95), 4),
            "stability": round(random.uniform(0.40, 0.95), 4),
            "notes": "Synthetic placeholder prediction.",
        }


class ESMFoldPredictor(StructurePredictor):
    def __init__(self, api_url: str, timeout_seconds: int) -> None:
        self.api_url = api_url
        self.timeout_seconds = timeout_seconds

    @staticmethod
    def _mean_plddt_from_pdb(pdb_text: str) -> float:
        parser = PDBParser(QUIET=True)
        from io import StringIO

        structure = parser.get_structure("esmfold", StringIO(pdb_text))
        b_factors: list[float] = []
        for atom in structure.get_atoms():
            b_factors.append(float(atom.get_bfactor()))

        if not b_factors:
            raise ValueError("No atoms found in ESMFold PDB response.")
        return sum(b_factors) / len(b_factors)

    def predict(self, sequence: str) -> dict:
        if not sequence:
            raise ValueError("Cannot predict structure for empty sequence.")

        response = requests.post(
            self.api_url,
            data=sequence,
            timeout=self.timeout_seconds,
        )
        response.raise_for_status()
        pdb_text = response.text.strip()
        if not pdb_text or "ATOM" not in pdb_text:
            raise ValueError("ESMFold API returned an invalid or empty PDB payload.")

        mean_plddt = self._mean_plddt_from_pdb(pdb_text)
        confidence = max(0.0, min(mean_plddt / 100.0, 1.0))
        # In this lightweight pipeline, stability is proxied from confidence.
        stability = confidence

        return {
            "backend": "esmfold",
            "confidence": round(confidence, 4),
            "stability": round(stability, 4),
            "plddt_mean": round(mean_plddt, 2),
            "pdb_text": pdb_text,
            "notes": "Real ESMFold inference via ESM Atlas API.",
        }


class RemotePdbPredictor(StructurePredictor):
    def __init__(
        self,
        backend_name: str,
        api_url: str,
        timeout_seconds: int,
        request_mode: str,
        json_sequence_key: str,
        response_pdb_key: str,
        response_confidence_key: str,
    ) -> None:
        if not api_url:
            raise ValueError(f"{backend_name} backend requires a non-empty API URL.")

        self.backend_name = backend_name
        self.api_url = api_url
        self.timeout_seconds = timeout_seconds
        self.request_mode = request_mode.strip().lower()
        self.json_sequence_key = json_sequence_key
        self.response_pdb_key = response_pdb_key
        self.response_confidence_key = response_confidence_key

    @staticmethod
    def _mean_plddt_from_pdb(pdb_text: str) -> float:
        parser = PDBParser(QUIET=True)
        from io import StringIO

        structure = parser.get_structure("remote_pdb", StringIO(pdb_text))
        b_factors: list[float] = [float(atom.get_bfactor()) for atom in structure.get_atoms()]
        if not b_factors:
            raise ValueError("No atoms found in PDB response.")
        return sum(b_factors) / len(b_factors)

    @staticmethod
    def _get_nested_value(payload: Any, dotted_key: str) -> Any:
        current = payload
        for part in dotted_key.split("."):
            if not isinstance(current, dict) or part not in current:
                return None
            current = current[part]
        return current

    def predict(self, sequence: str) -> dict:
        if not sequence:
            raise ValueError(f"Cannot predict structure for empty sequence ({self.backend_name}).")

        if self.request_mode == "json":
            response = requests.post(
                self.api_url,
                json={self.json_sequence_key: sequence},
                timeout=self.timeout_seconds,
            )
        else:
            response = requests.post(
                self.api_url,
                data=sequence,
                timeout=self.timeout_seconds,
            )

        response.raise_for_status()

        pdb_text: str | None = None
        confidence_raw: Any = None
        content_type = response.headers.get("content-type", "").lower()

        if "application/json" in content_type:
            payload = response.json()
            pdb_candidate = self._get_nested_value(payload, self.response_pdb_key)
            if isinstance(pdb_candidate, str):
                pdb_text = pdb_candidate.strip()
            confidence_raw = self._get_nested_value(payload, self.response_confidence_key)
        else:
            pdb_text = response.text.strip()

        if not pdb_text or "ATOM" not in pdb_text:
            raise ValueError(f"{self.backend_name} API returned an invalid or empty PDB payload.")

        mean_plddt = self._mean_plddt_from_pdb(pdb_text)
        confidence = max(0.0, min(mean_plddt / 100.0, 1.0))
        if confidence_raw is not None:
            try:
                confidence = max(0.0, min(float(confidence_raw), 1.0))
            except (TypeError, ValueError):
                pass

        return {
            "backend": self.backend_name,
            "confidence": round(confidence, 4),
            "stability": round(confidence, 4),
            "plddt_mean": round(mean_plddt, 2),
            "pdb_text": pdb_text,
            "notes": f"Real {self.backend_name} inference via configured API.",
        }


class AlphaFold3Predictor(RemotePdbPredictor):
    def __init__(self, workflow_config: dict[str, Any]) -> None:
        super().__init__(
            backend_name="alphafold3",
            api_url=str(workflow_config.get("alphafold3_api_url", "")).strip(),
            timeout_seconds=int(workflow_config.get("alphafold3_timeout_seconds", 300)),
            request_mode=str(workflow_config.get("alphafold3_request_mode", "json")),
            json_sequence_key=str(workflow_config.get("alphafold3_json_sequence_key", "sequence")),
            response_pdb_key=str(workflow_config.get("alphafold3_response_pdb_key", "pdb")),
            response_confidence_key=str(workflow_config.get("alphafold3_response_confidence_key", "confidence")),
        )


class RosettaPredictor(RemotePdbPredictor):
    def __init__(self, workflow_config: dict[str, Any]) -> None:
        super().__init__(
            backend_name="rosetta",
            api_url=str(workflow_config.get("rosetta_api_url", "")).strip(),
            timeout_seconds=int(workflow_config.get("rosetta_timeout_seconds", 300)),
            request_mode=str(workflow_config.get("rosetta_request_mode", "json")),
            json_sequence_key=str(workflow_config.get("rosetta_json_sequence_key", "sequence")),
            response_pdb_key=str(workflow_config.get("rosetta_response_pdb_key", "pdb")),
            response_confidence_key=str(workflow_config.get("rosetta_response_confidence_key", "confidence")),
        )


def configure_logger(log_path: Path, level_name: str) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=getattr(logging, level_name.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )
    return logging.getLogger(__name__)


def get_predictor(name: str, workflow_config: dict[str, Any]) -> StructurePredictor:
    normalized = name.strip().lower()
    if normalized == "mock":
        return MockPredictor()
    if normalized == "esmfold":
        api_url = str(workflow_config.get("esmfold_api_url", "https://api.esmatlas.com/foldSequence/v1/pdb/"))
        timeout_seconds = int(workflow_config.get("esmfold_timeout_seconds", 120))
        return ESMFoldPredictor(api_url=api_url, timeout_seconds=timeout_seconds)
    if normalized in {"alphafold", "alphafold3"}:
        return AlphaFold3Predictor(workflow_config)
    if normalized == "rosetta":
        return RosettaPredictor(workflow_config)
    raise ValueError(f"Unknown structure backend: {name}")


def predict_with_fallback(
    predictor: StructurePredictor,
    sequence: str,
    logger: logging.Logger,
    requested_backend: str,
) -> dict:
    try:
        result = predictor.predict(sequence)
        result["backend_requested"] = requested_backend
        result["backend_used"] = result.get("backend", requested_backend)
        return result
    except Exception as exc:
        logger.warning(
            "Backend '%s' failed (%s). Falling back to mock predictor.",
            requested_backend,
            exc,
        )
        fallback = MockPredictor().predict(sequence)
        fallback["backend_requested"] = requested_backend
        fallback["backend_used"] = "mock"
        fallback["fallback_reason"] = str(exc)
        return fallback


def main() -> None:
    input_path = Path(str(snakemake.input[0]))
    output_path = Path(str(getattr(snakemake.output, "predictions", snakemake.output[0])))
    pdb_dir = Path(str(getattr(snakemake.output, "pdb_dir", snakemake.config.get("structure_pdb_dir", "data/pdb"))))
    backend_name = str(getattr(snakemake.params, "backend", "mock"))
    log_path = Path(str(snakemake.log[0]))
    log_level = str(snakemake.config.get("log_level", "INFO"))

    logger = configure_logger(log_path, log_level)
    logger.info("Starting structure step with backend='%s'", backend_name)

    try:
        try:
            predictor = get_predictor(backend_name, dict(snakemake.config))
            logger.info("Using predictor backend: %s", backend_name)
        except Exception as exc:
            logger.warning("Requested backend '%s' is invalid (%s). Using mock.", backend_name, exc)
            predictor = MockPredictor()

        predictions: list[dict] = []
        pdb_dir.mkdir(parents=True, exist_ok=True)

        for record in SeqIO.parse(str(input_path), "fasta"):
            pred = predict_with_fallback(predictor, str(record.seq), logger, backend_name)
            pred["sequence_id"] = record.id

            pdb_text = pred.pop("pdb_text", None)
            if isinstance(pdb_text, str) and pdb_text.strip():
                pdb_path = pdb_dir / f"{record.id}.pdb"
                pdb_path.write_text(pdb_text, encoding="utf-8")
                pred["pdb_file"] = str(pdb_path.as_posix())

            predictions.append(pred)

        if not predictions:
            raise ValueError("No sequences found for structure prediction.")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(predictions, indent=2), encoding="utf-8")
        logger.info("Structure step completed for %d sequences.", len(predictions))
    except Exception as exc:
        logger.exception("Structure step failed: %s", exc)
        raise


if __name__ == "__main__":
    main()
