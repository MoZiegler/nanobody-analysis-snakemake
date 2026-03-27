# Nanobody In-Silico Pipeline (Snakemake, HPC-ready)

A reproducible, modular, and extensible Snakemake workflow for nanobody (VHH) sequence analysis, structure prediction, and developability scoring.

---

## Tech Stack

Python · Snakemake · Docker · HPC (SLURM) · ESMFold · BioPython

---

## Quickstart

```bash
git clone <repo>
cd nanobody-pipeline
python3.11 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
snakemake --cores 4 --config structure_backend=esmfold
```

---

## Example Output

Top-ranked sequences (example):

| rank | sequence_id | developability_score |
|------|------------|----------------------|
| 1    | seq_3      | 0.82                 |
| 2    | seq_1      | 0.78                 |
| 3    | seq_5      | 0.74                 |

---

## Highlights

* Modular protein structure prediction pipeline (ESMFold [1] / AlphaFold [2] / Rosetta [3])
* Reproducible Snakemake [4] workflow
* HPC-ready execution (SLURM / cluster support)
* Docker + devcontainer support
* Fault-tolerant backend system with fallback mechanisms
* Optional molecular dynamics integration (GROMACS [5] placeholder)

---

## Context

Nanobodies (VHH domains) are widely used in therapeutic and diagnostic applications [6].
Efficient computational pipelines are required to process large candidate sets, predict structure, and prioritize sequences for experimental validation.

This project demonstrates a modular, reproducible **in-silico pipeline for early-stage nanobody screening and prioritization**, reflecting real-world computational protein engineering workflows.

---

## Design Goals

* Modular backend abstraction for structure prediction (ESMFold, AlphaFold, Rosetta)
* Robustness via fallback mechanisms
* Reproducibility across local, containerized, and HPC environments
* Extensibility towards real-world protein design pipelines
* Clear separation of orchestration (Snakemake) and computation (Python modules)

---

## Why This Matters

This pipeline reflects core challenges in computational protein engineering:

- orchestrating heterogeneous tools (AI models, MD simulations)
- ensuring reproducibility across environments
- scaling to large candidate sets
- enabling rapid iteration in early-stage design

The architecture is designed to mirror real-world in-silico platforms.

---

## Pipeline Overview

```
FASTA
  ↓
Preprocess
  ↓
MSA (optional)
  ↓
Feature Extraction
  ↓
Structure Prediction (pluggable backend)
  ↓
Optional MD (GROMACS)
  ↓
Scoring
  ↓
Ranked Report
```

---

## Features

* Snakemake pipeline orchestration from raw FASTA to ranked report.
* Optional multiple sequence alignment (MSA) output for downstream inspection.
* Modular structure prediction backend system (Strategy + Factory).
* Optional HPC-friendly execution (`--cluster ...`) and CLI config overrides.
* Reproducible local lightweight setup with `venv`, plus optional Docker and devcontainer.
* Optional GROMACS placeholder simulation (Docker-only mode) that integrates into scoring.

---

## AI Integration

The pipeline supports AI-based structure prediction via ESMFold,
a transformer-based protein language model.

This enables integration of modern deep learning approaches into classical bioinformatics workflows.

---

## Structure Prediction Abstraction

Structure prediction is implemented via a pluggable backend system:

* Unified interface (`StructurePredictor`)
* Runtime backend selection via config/CLI
* Graceful degradation via fallback to mock backend
* API-based integration for external tools

This design mirrors real-world computational biology platforms where multiple tools are orchestrated and swapped depending on availability, cost, or compute constraints.

---

## HPC Compatibility

The pipeline is designed to run on HPC clusters using Snakemake's native cluster support:

```bash
snakemake --cluster "sbatch -A project -t 01:00:00"
```

This enables:

* distributed execution
* scalability to large sequence sets
* integration with SLURM-based environments

---

## Molecular Dynamics Integration (GROMACS)

An optional MD step is included as a placeholder for structure refinement workflows.

While implemented as a lightweight mock, the interface is designed to integrate real tools such as GROMACS for:

* energy minimization
* stability estimation
* trajectory-based scoring

This reflects typical post-structure prediction workflows in computational protein engineering.

---

## Design Tradeoffs

Heavy computational tools (AlphaFold, Rosetta, GROMACS) are not bundled directly.

Instead, this pipeline focuses on:

* clean integration interfaces
* reproducibility
* portability

This allows seamless replacement of placeholder logic with production-grade tools in real environments.

---

## Setup

### Local (venv + pip)

```bash
python3.11 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

---

### Optional Docker

Use Docker when you want the optional GROMACS step (`use_gromacs=true`).

```bash
docker build -f docker/Dockerfile -t nanobody-pipeline .
docker run --rm nanobody-pipeline
```

---

### Optional VS Code Devcontainer

Open the folder in VS Code and choose **Reopen in Container**.
Also supports the optional GROMACS step (`use_gromacs=true`).

---

## Usage

### Local mock run

```bash
snakemake --cores 4
```

### Use structure prediction backend

```bash
snakemake --cores 4 --config structure_backend=esmfold
```

### HPC execution

```bash
snakemake --cluster "sbatch -A project -t 01:00:00"
```

---

## Architecture

Pipeline steps:

1. preprocess
2. msa (optional)
3. msa_summary (optional)
4. features
5. structure
6. gromacs_simulation (optional)
7. scoring

---

## Configuration

`config.yaml` defaults:

```yaml
structure_backend: "mock"
use_gromacs: false
run_msa: true
esmfold_api_url: "https://api.esmatlas.com/foldSequence/v1/pdb/"
esmfold_timeout_seconds: 120
alphafold3_api_url: "" # provide URL if you have a local AlphaFold3 instance
alphafold3_timeout_seconds: 300
alphafold3_request_mode: "json"
alphafold3_json_sequence_key: "sequence"
alphafold3_response_pdb_key: "pdb"
alphafold3_response_confidence_key: "confidence"
rosetta_api_url: "" # provide URL if you have a local Rosetta instance
rosetta_timeout_seconds: 300
rosetta_request_mode: "json"
rosetta_json_sequence_key: "sequence"
rosetta_response_pdb_key: "pdb"
rosetta_response_confidence_key: "confidence"
structure_pdb_dir: "data/pdb"
log_level: "INFO"
```

You can override values from CLI using `--config`.

---

## Outputs

* `results/report.csv`: ranked developability report
* `data/msa.fasta`: optional alignment
* `data/msa_summary.csv`: alignment statistics
* `data/pdb/*.pdb`: predicted structures
* `logs/*.log`: execution logs

---

## Developability Score

The pipeline computes a lightweight composite developability score:

* stability (0.4)
* hydrophobicity (0.3)
* charge balance (0.3)

This approximates early-stage developability screening used in antibody engineering pipelines [7], [8], [9], [10].

---

## Notebook

`notebooks/analysis.ipynb` provides:

* score visualization
* ranking inspection
* MSA analysis
* structure visualization

---

## Scalability

The pipeline is designed to scale to large sequence sets:

- parallel execution via Snakemake
- backend abstraction allows distribution across compute nodes
- compatible with HPC schedulers (e.g. SLURM)

This enables processing of large candidate libraries in realistic settings.

---

## Notes

* Structure prediction and MD simulation are intentionally lightweight placeholders to keep the pipeline runnable everywhere.
* Local `venv` mode is a lightweight Snakemake workflow; optional GROMACS execution is Docker-only.
* Heavy scientific tools (ESMFold, AlphaFold, Rosetta, GROMACS) can be integrated by replacing stub logic while preserving interfaces.

---

## Limitations

- Developability score is a simplified heuristic
- Structure-based stability is approximated
- MSA implementation is lightweight (not MAFFT-level)
- MD step is a placeholder

The focus of this project is pipeline architecture, not biological accuracy.

---

## Future Work

- Integration of real MD simulations (GROMACS / OpenMM)
- Improved scoring using structural features (SASA, RMSD)
- Integration with protein design models
- Distributed execution across HPC clusters

---

## Summary

This project demonstrates how modern protein engineering workflows can be:

* modular
* reproducible
* scalable
* backend-agnostic

and integrated into a clean, extensible pipeline architecture.

---

## References

[1]	Z. Lin et al., “Evolutionary-scale prediction of atomic-level protein structure with a language model,” Science, vol. 379, no. 6637, pp. 1123–1130, Mar. 2023, doi: 10.1126/science.ade2574.

[2]	J. Abramson et al., “Accurate structure prediction of biomolecular interactions with AlphaFold 3,” Nature, vol. 630, no. 8016, pp. 493–500, Jun. 2024, doi: 10.1038/s41586-024-07487-w

[3]	J. K. Leman et al., “Macromolecular modeling and design in Rosetta: recent methods and frameworks,” Nat Methods, vol. 17, no. 7, pp. 665–680, Jul. 2020, doi: 10.1038/s41592-020-0848-2.

[4]	F. Mölder et al., “Sustainable data analysis with Snakemake,” Sep. 23, 2025, F1000Research. doi: 10.12688/f1000research.29032.3.

[5]	M. Abraham et al., “GROMACS 2026.1 Manual,” Mar. 2026, doi: 10.5281/zenodo.18886967.

[6]	S. Muyldermans, “Applications of Nanobodies,” Annual Review of Animal Biosciences, vol. 9, no. Volume 9, 2021, pp. 401–421, Feb. 2021, doi: 10.1146/annurev-animal-021419-083831.

[7]	J. Kyte and R. F. Doolittle, “A simple method for displaying the hydropathic character of a protein,” J Mol Biol, vol. 157, no. 1, pp. 105–132, May 1982, doi: 10.1016/0022-2836(82)90515-0

[8]	T. M. Lauer, N. J. Agrawal, N. Chennamsetty, K. Egodage, B. Helk, and B. L. Trout, “Developability index: a rapid in silico tool for the screening of antibody aggregation propensity,” J Pharm Sci, vol. 101, no. 1, pp. 102–115, Jan. 2012, doi: 10.1002/jps.22758.

[9]	M. I. J. Raybould et al., “Five computational developability guidelines for therapeutic antibody profiling,” Proc Natl Acad Sci U S A, vol. 116, no. 10, pp. 4025–4030, Mar. 2019, doi: 10.1073/pnas.1810576116.

[10] J. Jumper et al., “Highly accurate protein structure prediction with AlphaFold,” Nature, vol. 596, no. 7873, pp. 583–589, Aug. 2021, doi: 10.1038/s41586-021-03819-2.