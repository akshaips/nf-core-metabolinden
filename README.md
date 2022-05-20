# ![nf-core/metabolinden](docs/images/nf-core-metabolinden_logo.png)

**Metabolomics quaLIty coNtrol anD paramEter optimizatioN**.

[![GitHub Actions CI Status](https://github.com/payamemami/metabolinden/workflows/nf-core%20CI/badge.svg)](https://github.com/payamemami/metabolinden/actions)
[![GitHub Actions Linting Status](https://github.com/payamemami/metabolinden/workflows/nf-core%20linting/badge.svg)](https://github.com/payamemami/metabolinden/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/payamemami/metabolinden.svg)](https://hub.docker.com/r/nfcore/metabolinden)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23metabolinden-4A154B?logo=slack)](https://nfcore.slack.com/channels/metabolinden)

## Introduction

**metabolinden** is a bioinformatics best-practise analysis pipeline for metabolomics data pre-processing and parameter tuning.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command (not usable now!):

    ```bash
    nextflow run payamemami/metabolinden -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

4. Start running your own analysis!

    ```bash
    nextflow run payamemami/metabolinden -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input '*.mzML' --identification_input 'database.tsv' --recalibration_masses 'lock_in_masses.csv'
    ```

See [usage docs](docs/) for all of the available options when running the pipeline.

## Pipeline Summary

The pipeline currently performs the following:

* Centroiding
* Recalibration
* Feature detection
* Alignment
* Grouping
* Identification based on internal standard
* Data exporting

## Documentation

The nf-core/metabolinden pipeline comes with documentation about the pipeline: [usage](docs/usage.md) and [output](docs/output.md).

## Credits

payamemami/metabolinden was originally written by Payam Emami.

We thank the following people for their extensive assistance in the development
of this pipeline:

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#metabolinden` channel](https://nfcore.slack.com/channels/metabolinden) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/metabolinden for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

In addition, references of tools and data used in this pipeline are as follows:

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->
