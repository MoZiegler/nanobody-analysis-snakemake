import logging
from pathlib import Path

import pandas as pd
from Bio import AlignIO


def configure_logger(log_path: Path, level_name: str) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=getattr(logging, level_name.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )
    return logging.getLogger(__name__)


def most_common_residue(column: str) -> str:
    residues = [aa for aa in column if aa != "-"]
    if not residues:
        return "-"
    counts: dict[str, int] = {}
    for aa in residues:
        counts[aa] = counts.get(aa, 0) + 1
    return max(counts.items(), key=lambda item: item[1])[0]


def main() -> None:
    input_path = Path(str(snakemake.input[0]))
    output_path = Path(str(snakemake.output[0]))
    log_path = Path(str(snakemake.log[0]))
    log_level = str(snakemake.config.get("log_level", "INFO"))

    logger = configure_logger(log_path, log_level)
    logger.info("Starting MSA summary step: %s -> %s", input_path, output_path)

    try:
        alignment = AlignIO.read(str(input_path), "fasta")
        n_sequences = len(alignment)
        aln_len = alignment.get_alignment_length()

        if n_sequences == 0 or aln_len == 0:
            raise ValueError("Alignment is empty and cannot be summarized.")

        consensus = "".join(most_common_residue(alignment[:, i]) for i in range(aln_len))

        rows: list[dict[str, float | int | str]] = []
        for record in alignment:
            aligned_seq = str(record.seq)
            gap_count = aligned_seq.count("-")
            nongap_count = aln_len - gap_count

            matches = 0
            for aa, con in zip(aligned_seq, consensus):
                if aa != "-" and con != "-" and aa == con:
                    matches += 1
            consensus_identity = (matches / nongap_count) if nongap_count else 0.0

            rows.append(
                {
                    "sequence_id": record.id,
                    "alignment_length": aln_len,
                    "nongap_length": nongap_count,
                    "gap_count": gap_count,
                    "gap_fraction": gap_count / aln_len,
                    "consensus_identity": consensus_identity,
                }
            )

        summary_df = pd.DataFrame(rows).sort_values("sequence_id").reset_index(drop=True)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(output_path, index=False)
        logger.info("MSA summary completed for %d sequences.", len(summary_df))
    except Exception as exc:
        logger.exception("MSA summary failed: %s", exc)
        raise


if __name__ == "__main__":
    main()
