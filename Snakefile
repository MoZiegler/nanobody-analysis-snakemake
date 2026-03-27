from pathlib import Path

configfile: "config.yaml"


def as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


USE_GROMACS = as_bool(config.get("use_gromacs", False))
RUN_MSA = as_bool(config.get("run_msa", True))
SCORING_MD_INPUT = "data/md_results.json" if USE_GROMACS else []
ALL_TARGETS = ["results/report.csv"]

if RUN_MSA:
    ALL_TARGETS.append("data/msa.fasta")
    ALL_TARGETS.append("data/msa_summary.csv")


rule all:
    input:
        ALL_TARGETS


rule preprocess:
    input:
        "data/input.fasta"
    output:
        "data/cleaned.fasta"
    log:
        "logs/preprocess.log"
    script:
        "scripts/preprocess.py"


rule features:
    input:
        "data/cleaned.fasta"
    output:
        "data/features.csv"
    log:
        "logs/features.log"
    script:
        "scripts/features.py"


rule structure:
    input:
        "data/cleaned.fasta"
    output:
        predictions="data/structures.json",
        pdb_dir=directory("data/pdb")
    params:
        backend=config.get("structure_backend", "mock")
    log:
        "logs/structure.log"
    script:
        "scripts/structure.py"


if RUN_MSA:
    rule msa:
        input:
            "data/cleaned.fasta"
        output:
            "data/msa.fasta"
        log:
            "logs/msa.log"
        script:
            "scripts/msa.py"


if RUN_MSA:
    rule msa_summary:
        input:
            "data/msa.fasta"
        output:
            "data/msa_summary.csv"
        log:
            "logs/msa_summary.log"
        script:
            "scripts/msa_summary.py"


if USE_GROMACS:
    rule gromacs_simulation:
        input:
            "data/structures.json"
        output:
            "data/md_results.json"
        log:
            "logs/gromacs.log"
        script:
            "scripts/gromacs.py"


rule scoring:
    input:
        features="data/features.csv",
        structures="data/structures.json",
        md=SCORING_MD_INPUT
    output:
        "results/report.csv"
    log:
        "logs/scoring.log"
    script:
        "scripts/scoring.py"
