# from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def configure_logger(log_path: Path, level_name: str) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=getattr(logging, level_name.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )
    return logging.getLogger(__name__)


def normalize_hydrophobicity(value: float) -> float:
    normalized = (value + 4.5) / 9.0
    return float(np.clip(normalized, 0.0, 1.0))


def get_optional_md_path(md_input: Any) -> Path | None:
    if isinstance(md_input, str) and md_input:
        return Path(md_input)
    if isinstance(md_input, (list, tuple)) and md_input:
        return Path(str(md_input[0]))
    return None


def main() -> None:
    features_path = Path(str(snakemake.input.features))
    structures_path = Path(str(snakemake.input.structures))
    md_path = get_optional_md_path(snakemake.input.md)
    output_path = Path(str(snakemake.output[0]))
    log_path = Path(str(snakemake.log[0]))
    log_level = str(snakemake.config.get("log_level", "INFO"))

    logger = configure_logger(log_path, log_level)
    logger.info("Starting scoring step.")

    try:
        features_df = pd.read_csv(features_path)
        structures_data = json.loads(structures_path.read_text(encoding="utf-8"))
        structures_df = pd.DataFrame(structures_data)

        if features_df.empty or structures_df.empty:
            raise ValueError("Features or structures input is empty.")

        merged = features_df.merge(structures_df, on="sequence_id", how="inner", validate="one_to_one")
        if merged.empty:
            raise ValueError("No matching sequence IDs between features and structures.")

        merged["hydrophobicity_norm"] = merged["hydrophobicity"].astype(float).map(normalize_hydrophobicity)
        merged["charge_balance"] = merged["charge_balance"].astype(float)
        merged["stability"] = merged["stability"].astype(float)

        merged["developability_score"] = (
            0.4 * merged["stability"]
            + 0.3 * merged["hydrophobicity_norm"]
            + 0.3 * merged["charge_balance"]
        )

        if md_path is not None and md_path.exists():
            logger.info("Integrating optional MD data from %s", md_path)
            md_data = json.loads(md_path.read_text(encoding="utf-8"))
            md_df = pd.DataFrame(md_data)
            if not md_df.empty and "sequence_id" in md_df.columns:
                merged = merged.merge(md_df, on="sequence_id", how="left")
                merged["md_stability_adjustment"] = merged["md_stability_adjustment"].fillna(0.0).astype(float)
                merged["developability_score"] = np.clip(
                    merged["developability_score"] + 0.1 * merged["md_stability_adjustment"],
                    0.0,
                    1.0,
                )
        else:
            logger.info("Optional GROMACS data not provided; scoring without MD adjustments.")

        merged = merged.sort_values("developability_score", ascending=False).reset_index(drop=True)
        merged.insert(0, "rank", merged.index + 1)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(output_path, index=False)
        logger.info("Scoring completed for %d sequences.", len(merged))
    except Exception as exc:
        logger.exception("Scoring failed: %s", exc)
        raise


if __name__ == "__main__":
    main()
