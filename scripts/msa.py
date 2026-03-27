import logging
from pathlib import Path

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def configure_logger(log_path: Path, level_name: str) -> logging.Logger:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=getattr(logging, level_name.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )
    return logging.getLogger(__name__)


def add_sequence_to_alignment(
    aligned_sequences: dict[str, str],
    reference_id: str,
    new_id: str,
    new_sequence: str,
) -> dict[str, str]:
    """Progressively merge one sequence into an existing alignment using a reference-guided pairwise alignment."""
    old_reference_alignment = aligned_sequences[reference_id]
    ungapped_reference = old_reference_alignment.replace("-", "")

    pairwise = pairwise2.align.globalms(
        ungapped_reference,
        new_sequence,
        2.0,
        -1.0,
        -5.0,
        -0.5,
        one_alignment_only=True,
    )
    if not pairwise:
        raise ValueError(f"Pairwise alignment failed for sequence '{new_id}'.")

    new_reference_alignment = pairwise[0].seqA
    new_sequence_alignment = pairwise[0].seqB

    merged: dict[str, list[str]] = {seq_id: [] for seq_id in aligned_sequences}
    merged[new_id] = []

    old_idx = 0
    new_idx = 0

    while old_idx < len(old_reference_alignment) or new_idx < len(new_reference_alignment):
        old_char = old_reference_alignment[old_idx] if old_idx < len(old_reference_alignment) else None
        new_char = new_reference_alignment[new_idx] if new_idx < len(new_reference_alignment) else None

        if old_char == "-" and new_char == "-":
            for seq_id, seq in aligned_sequences.items():
                merged[seq_id].append(seq[old_idx])
            merged[new_id].append("-")
            old_idx += 1
            new_idx += 1
            continue

        if old_char == "-":
            for seq_id, seq in aligned_sequences.items():
                merged[seq_id].append(seq[old_idx])
            merged[new_id].append("-")
            old_idx += 1
            continue

        if new_char == "-":
            for seq_id in aligned_sequences:
                merged[seq_id].append("-")
            merged[new_id].append(new_sequence_alignment[new_idx])
            new_idx += 1
            continue

        if old_char is None and new_char is not None:
            if new_char != "-":
                raise ValueError("Inconsistent alignment state while appending new sequence.")
            for seq_id in aligned_sequences:
                merged[seq_id].append("-")
            merged[new_id].append(new_sequence_alignment[new_idx])
            new_idx += 1
            continue

        if new_char is None and old_char is not None:
            if old_char != "-":
                raise ValueError("Inconsistent alignment state while extending existing alignment.")
            for seq_id, seq in aligned_sequences.items():
                merged[seq_id].append(seq[old_idx])
            merged[new_id].append("-")
            old_idx += 1
            continue

        if old_char is None and new_char is None:
            break

        for seq_id, seq in aligned_sequences.items():
            merged[seq_id].append(seq[old_idx])
        merged[new_id].append(new_sequence_alignment[new_idx])
        old_idx += 1
        new_idx += 1

    merged_strings = {seq_id: "".join(chars) for seq_id, chars in merged.items()}
    lengths = {len(seq) for seq in merged_strings.values()}
    if len(lengths) != 1:
        raise ValueError("Merged alignment has inconsistent sequence lengths.")

    return merged_strings


def progressive_msa(records: list[SeqRecord], logger: logging.Logger) -> dict[str, str]:
    if not records:
        raise ValueError("No sequences were provided for MSA.")

    reference_id = records[0].id
    alignment: dict[str, str] = {reference_id: str(records[0].seq).upper()}

    for record in records[1:]:
        sequence = str(record.seq).upper()
        logger.info("Aligning sequence '%s' into progressive MSA.", record.id)
        alignment = add_sequence_to_alignment(alignment, reference_id, record.id, sequence)

    return alignment


def main() -> None:
    input_path = Path(str(snakemake.input[0]))
    output_path = Path(str(snakemake.output[0]))
    log_path = Path(str(snakemake.log[0]))
    log_level = str(snakemake.config.get("log_level", "INFO"))

    logger = configure_logger(log_path, log_level)
    logger.info("Starting MSA step: %s -> %s", input_path, output_path)

    try:
        records = list(SeqIO.parse(str(input_path), "fasta"))
        if not records:
            raise ValueError("Input FASTA is empty after preprocessing.")

        alignment = progressive_msa(records, logger)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        aligned_records = [
            SeqRecord(Seq(aligned_seq), id=seq_id, description="aligned")
            for seq_id, aligned_seq in alignment.items()
        ]
        SeqIO.write(aligned_records, str(output_path), "fasta")

        alignment_length = len(next(iter(alignment.values())))
        logger.info(
            "MSA completed for %d sequences (alignment length: %d).",
            len(aligned_records),
            alignment_length,
        )
    except Exception as exc:
        logger.exception("MSA step failed: %s", exc)
        raise


if __name__ == "__main__":
    main()
