# Configuring the Definition File and Running the Evaluator Container for MPRALegNet

## Overview

This container includes:
    - Evaluator API script for genomic sequence evaluation.
    - Installed Python dependencies required by the Evaluator.

## Usage

We encourage using pre-built containers for this model that are hosted on Hugging Face: https://huggingface.co/Ni-os/MPRALegNet-container.

However, if you are building the container using the provided definition file, ensure you have the following directory structure on the host:

```bash
/path/to/MPRALegNet/src/evaluator_container_apptainer/
    ├── evaluator.def                                # Evaluator container definition file
    ├── evaluator.sif                                # Evaluator container image (not included on Github due to storage restrictions)
    ├── MPRALegNet_evaluator_API.py                  # Evaluator API script
    ├── data_loader.py                               # Data loading and validation
    ├── config.py                                    # Configuration file
    ├── evaluator_content_handler.py                 # API communication handlers
    ├── evaluator_metrics_calculator.py              # Correlation and specificity calculations
    ├── evaluator_data                
    │   ├── evaluator_input_sample_test.json         # Sample input JSON for evaluator
    │   ├── evaluator_message_gosai_5seqs.json       # Gosai 5 test sequences
    │   ├── evaluator_message_more_complex.json      # Complex evaluator input example
    │   └── evaluator_message_simple_test.json       # Simple evaluator input example
    └── predictions/                                 # Output directory for prediction results

```

### Build the container (SIF)

```bash
apptainer build legnet_evaluator.sif evaluator.def
```

### Bind these directories and run the container

The recommended way to run the evaluator is with the `--containall` flag for isolation, combined with the `-B` flag to bind the necessary data directories.

```bash
apptainer run --containall \
    -B absolute/path/to/evaluator_data:/evaluator_data \
    -B absolute/path/to/predictions:/predictions \
    legnet_evaluator.sif PREDICTOR_HOST PREDICTOR_PORT /predictions
```

## Details

- The container sends and receives data via a RestAPI and requires `/evaluator_data` and `/predictions` to be mounted at runtime.

## Purpose

- Facilitates genomic sequence evaluation by interacting with the Predictor container.
- Handles input/output through bound directories.

## Example Command

```bash
apptainer run --containall \
    -B absolute/path/to/evaluator_data:/evaluator_data \
    -B absolute/path/to/predictions:/predictions \
    legnet_evaluator.sif 172.16.47.244 5000 /predictions
```

## Arguments

1. `PREDICTOR_HOST`: IP address or hostname of the Predictor server.
2. `PREDICTOR_PORT`: Port number the Predictor is listening on.
3. `/predictions`: Path inside the container for storing output predictions.