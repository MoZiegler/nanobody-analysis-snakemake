# from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import molecular_weight


AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
HYDRO_SCALE = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}


def configure_logger(log_path: Path, level_name: str) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=getattr(logging, level_name.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )
    return logging.getLogger(__name__)


def amino_acid_frequencies(sequence: str) -> dict[str, float]:
    length = len(sequence)
    if length == 0:
        return {f"freq_{aa}": 0.0 for aa in AA_ORDER}
    return {f"freq_{aa}": sequence.count(aa) / length for aa in AA_ORDER}


def average_hydrophobicity(sequence: str) -> float:
    if not sequence:
        return 0.0
    return sum(HYDRO_SCALE.get(aa, 0.0) for aa in sequence) / len(sequence)


def net_charge(sequence: str) -> float:
    pos = sequence.count("K") + sequence.count("R") + 0.1 * sequence.count("H")
    neg = sequence.count("D") + sequence.count("E")
    return float(pos - neg)


def charge_balance_from_net(net: float) -> float:
    return float(1.0 / (1.0 + abs(net)))


def extract_features(sequence_id: str, sequence: str) -> dict[str, float | str]:
    seq_len = len(sequence)
    hydro = average_hydrophobicity(sequence)
    charge = net_charge(sequence)
    features: dict[str, float | str] = {
        "sequence_id": sequence_id,
        "length": seq_len,
        "hydrophobicity": hydro,
        "net_charge": charge,
        "charge_balance": charge_balance_from_net(charge),
        "molecular_weight": float(molecular_weight(sequence, seq_type="protein")),
    }
    features.update(amino_acid_frequencies(sequence))
    return features


def main() -> None:
    input_path = Path(str(snakemake.input[0]))
    output_path = Path(str(snakemake.output[0]))
    log_path = Path(str(snakemake.log[0]))
    log_level = str(snakemake.config.get("log_level", "INFO"))

    logger = configure_logger(log_path, log_level)
    logger.info("Starting feature extraction: %s -> %s", input_path, output_path)

    try:
        rows: list[dict[str, float | str]] = []
        for record in SeqIO.parse(str(input_path), "fasta"):
            seq = str(record.seq).upper()
            rows.append(extract_features(record.id, seq))

        if not rows:
            raise ValueError("No sequences found for feature extraction.")

        output_path.parent.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)
        logger.info("Feature extraction completed for %d sequences.", len(df))
    except Exception as exc:
        logger.exception("Feature extraction failed: %s", exc)
        raise


if __name__ == "__main__":
    main()
