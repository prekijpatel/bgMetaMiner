# bgMetaMiner

**Retrieve, normalize, explore, and visualize prokaryotic genome metadata from NCBI.**

bgMetaMiner is a metadata-mining tool designed to simplify the retrieval, normalization, and exploration of metadata associated with prokaryotic genome assemblies hosted by NCBI.

It provides automated metadata extraction and normalization, including **NLP-based normalization of isolation-source metadata**, and an interactive dashboard for exploring large genome collections.

![Graphical abstract](https://github.com/prekijpatel/bgMetaMiner/blob/main/img/graphical%20abstract.png)

> **Release status**
>
> bgMetaMiner is under active development. The current version introduces NLP-based normalization and several changes to the metadata-processing workflow.
>
> Windows users can use the prebuilt executable, while Linux/Ubuntu users can run bgMetaMiner directly from source using the supplied Conda environment.

---

## Features

bgMetaMiner can be used to:

* retrieve prokaryotic genome metadata directly from NCBI;
* process locally downloaded NCBI genome metadata files;
* extract and organize assembly, sequencing, geographic, isolation-source, and related metadata;
* normalize geographic information and other heterogeneous metadata fields;
* normalize **host, source category, source, and sample information using an NLP-based model**;
* standardize sequencing technologies and assembly-related attributes;
* explore metadata interactively using an integrated dashboard;
* filter and visualize large genome collections using multiple metadata attributes;
* export normalized metadata for downstream genomic and epidemiological analyses.

For detailed explanations of the workflow, output fields, dashboard, and normalization procedures, please see the [bgMetaMiner Wiki](https://github.com/prekijpatel/bgMetaMiner/wiki).

---

# Installation

bgMetaMiner can be used in two ways:

| Platform           | Recommended installation                |
| ------------------ | --------------------------------------- |
| **Windows**        | Prebuilt executable                     |
| **Linux / Ubuntu** | Source distribution + Conda environment |

> A GPU is **not required**. The distributed environment uses CPU-based PyTorch for NLP inference.

---

## Windows

For most Windows users, the prebuilt executable is the simplest way to use bgMetaMiner.

### 1. Download bgMetaMiner

Download the latest Windows release from:

[**bgMetaMiner Releases**](https://github.com/prekijpatel/bgMetaMiner/releases)

### 2. Run bgMetaMiner

Launch the downloaded `bgMetaMiner.exe`.

No separate Python installation is required when using the prebuilt executable.

---

## Linux / Ubuntu

Linux users can run bgMetaMiner directly from the source distribution using the supplied `environment.yml`.

### Prerequisites

You will need:

* a 64-bit Linux system;
* [Conda](https://docs.conda.io/) or Miniconda installed;
* an internet connection during environment creation and for downloading metadata from NCBI.

### 1. Download and extract bgMetaMiner

Download the Linux/source `.tar.gz` archive from the [bgMetaMiner Releases](https://github.com/prekijpatel/bgMetaMiner/releases) page.

Then extract it:

```bash
tar -xzf bgMetaMiner-<version>.tar.gz
cd bgMetaMiner-<version>
```

### 2. Create the Conda environment

The archive contains an `environment.yml` with the tested versions of the required dependencies.

Create the environment using:

```bash
conda env create -f environment.yml
```

Activate it:

```bash
conda activate metaminer
```

> The first installation may require some time because the environment includes PyTorch, Transformers, Polars, Dash, Plotly, and other runtime dependencies.

### 3. Verify the installation

From inside the extracted bgMetaMiner directory, run:

```bash
python metaminer_CLI.py -h
```

You should see the bgMetaMiner command-line help page.

**Run bgMetaMiner from the extracted bgMetaMiner directory**, because the program uses bundled model and data files included with the distribution.

---

# Quick Start

The Linux command-line interface supports two main workflows:

1. processing an existing NCBI metadata file; and
2. downloading metadata from NCBI and then processing it.

To see all available options:

```bash
python metaminer_CLI.py -h
```

To see help for a specific command:

```bash
python metaminer_CLI.py normalize -h
```

or:

```bash
python metaminer_CLI.py download -h
```

---

## Normalize an existing NCBI metadata file

If you already have an NCBI genome metadata JSON or JSONL file:

```bash
python metaminer_CLI.py normalize \
    --input metadata.json \
    --outdir results
```

bgMetaMiner will extract and normalize the available genome metadata and write the generated files to the specified output directory.

---

## Download metadata for a taxon

bgMetaMiner can retrieve genome metadata directly using the NCBI Datasets command-line tool.

For example:

```bash
python metaminer_CLI.py download \
    --taxon "Escherichia coli" \
    --outdir results
```

A genus can also be supplied:

```bash
python metaminer_CLI.py download \
    --taxon Escherichia \
    --outdir results
```

NCBI Taxonomy IDs are also accepted:

```bash
python metaminer_CLI.py download \
    --taxon 562 \
    --outdir results
```

---

## Download metadata for an assembly accession

For a specific NCBI assembly:

```bash
python metaminer_CLI.py download \
    --accession GCA_000005825.2 \
    --outdir results
```

---

> ⚡ **Large datasets**
>
> For species with very large numbers of publicly available genomes, such as *Escherichia coli*, *Salmonella*, and other taxa with tens of thousands of assemblies, **16 GB RAM or more is strongly recommended**.

The amount of memory and processing time required will depend on the number of genomes being analysed and the completeness of their associated metadata.

---

# Interactive Metadata Exploration

bgMetaMiner includes an interactive dashboard for exploring normalized metadata.

Depending on the available metadata, users can examine characteristics such as:

* geographic distribution;
* assembly level;
* sequencing technology;
* submission year;
* genome assembly statistics;
* annotation source;
* isolation-source categories;
* NLP prediction confidence; and
* other assembly and BioSample attributes.

Filters can be combined to interactively investigate subsets of the downloaded genome collection.

See the [Dashboard Components](https://github.com/prekijpatel/MetaMiner/wiki/Components-of-Dashboard) for detailed dashboard documentation and examples.

---

This README is intended to provide installation instructions and a quick start.

For detailed documentation covering bgMetaMiner's workflow, metadata fields, normalization procedures, dashboard, output files, and usage examples, please visit:

### 📖 [bgMetaMiner Wiki](https://github.com/prekijpatel/bgMetaMiner/wiki)

---

# Feedback and Issues

bgMetaMiner is actively developed, and microbial metadata can contain many unusual or previously unseen formats.

If you encounter a bug, unexpected normalization, installation problem, or have a feature request, please open an issue.

# License

Please refer to the [`LICENSE`](https://github.com/prekijpatel/bgMetaMiner/blob/main/LICENSE) and accompanying repository files for the terms governing use and redistribution of bgMetaMiner.

---

## citation
If you find bgMetaMiner useful in your study, kindly cite. 
> Patel, J. K., & Elangovan, R. (2025). MetaMiner: Streamlined GUI Tool for Retrieving, Normalizing and Exploring Metadata. bioRxiv, 2025-08. https://doi.org/10.1101/2025.08.20.666107
