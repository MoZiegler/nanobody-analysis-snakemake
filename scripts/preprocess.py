# from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def configure_logger(log_path: Path, level_name: str) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=getattr(logging, level_name.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )
    return logging.getLogger(__name__)


def clean_sequence(seq: str) -> str:
    return "".join(aa for aa in seq.upper() if aa in VALID_AA)


def clean_records(records: Iterable[SeqRecord], logger: logging.Logger) -> list[SeqRecord]:
    cleaned: list[SeqRecord] = []
    for record in records:
        seq_clean = clean_sequence(str(record.seq))
        if not seq_clean:
            logger.warning("Skipping %s because no valid amino acids remained after cleaning.", record.id)
            continue
        new_record = SeqRecord(Seq(seq_clean), id=record.id, description=record.description)
        cleaned.append(new_record)
    return cleaned


def main() -> None:
    input_path = Path(str(snakemake.input[0]))
    output_path = Path(str(snakemake.output[0]))
    log_path = Path(str(snakemake.log[0]))
    log_level = str(snakemake.config.get("log_level", "INFO"))

    logger = configure_logger(log_path, log_level)
    logger.info("Starting preprocessing: %s -> %s", input_path, output_path)

    try:
        records = list(SeqIO.parse(str(input_path), "fasta"))
        if not records:
            raise ValueError(f"Input FASTA is empty: {input_path}")

        cleaned = clean_records(records, logger)
        if not cleaned:
            raise ValueError("All sequences were removed during cleaning.")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(cleaned, str(output_path), "fasta")
        logger.info("Preprocessing completed. Retained %d/%d sequences.", len(cleaned), len(records))
    except Exception as exc:
        logger.exception("Preprocessing failed: %s", exc)
        raise


if __name__ == "__main__":
    main()
