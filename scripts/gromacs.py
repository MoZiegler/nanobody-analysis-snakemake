# from __future__ import annotations

import json
import logging
import random
from pathlib import Path


def configure_logger(log_path: Path, level_name: str) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=getattr(logging, level_name.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )
    return logging.getLogger(__name__)


def main() -> None:
    input_path = Path(str(snakemake.input[0]))
    output_path = Path(str(snakemake.output[0]))
    log_path = Path(str(snakemake.log[0]))
    log_level = str(snakemake.config.get("log_level", "INFO"))

    logger = configure_logger(log_path, log_level)
    logger.info("Starting optional GROMACS placeholder simulation.")

    try:
        structures = json.loads(input_path.read_text(encoding="utf-8"))
        if not isinstance(structures, list) or not structures:
            raise ValueError("Structure input is empty or invalid.")

        md_results: list[dict[str, float | str]] = []
        for item in structures:
            seq_id = str(item.get("sequence_id", "unknown"))
            md_results.append(
                {
                    "sequence_id": seq_id,
                    "md_rmsd": round(random.uniform(1.0, 3.5), 4),
                    "md_stability_adjustment": round(random.uniform(-0.08, 0.08), 4),
                    "md_backend": "gromacs-mock",
                }
            )

        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(md_results, indent=2), encoding="utf-8")
        logger.info("Optional GROMACS placeholder completed for %d sequences.", len(md_results))
    except Exception as exc:
        logger.exception("GROMACS placeholder step failed: %s", exc)
        raise


if __name__ == "__main__":
    main()
