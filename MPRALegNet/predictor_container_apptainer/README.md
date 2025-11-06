# Configuring Definition File and Running Predictor Container for MPRALegNet

## Overview

This container provides a complete MPRALegNet Predictor implementation for genomic sequence expression prediction. The container includes:

- **REST API Server**: Flask-based API for handling prediction requests
- **MPRALegNet Model**: Pre-trained deep learning model for gene expression prediction
- **Sequence Processing**: Comprehensive preprocessing pipeline for DNA sequences
- **Error Handling**: Robust validation and error reporting
- **Dependencies**: Complete Conda environment with all required packages

> **Note**: The Conda environment (`legnet_environment.yml`) excludes PyTorch due to container build conflicts. PyTorch is installed manually in the container definition file.

## Usage

### Using Pre-built Containers

We recommend using the pre-built containers hosted on Zenodo for immediate use: **[Zenodo Link]**

### Building from Source

If you need to build the container locally, ensure you have the following directory structure:

```bash
/path/to/MPRALegNet/src/predictor_container_apptainer/
    ├── predictor.def                          # Predictor container definition file
    ├── legnet_predictor.sif                   # Predictor container image (not included on Github due to storage restrictions)
    ├── legnet_environment.yml                 # Conda environment specification
    ├── script_and_utils/
    │   ├── MPRALegNet_predictor_API.py        # Main Flask API server
    │   ├── MPRALegNet_predict.py              # Model prediction functions
    │   ├── model.py                           # LegNet model architecture
    │   ├── api_preprocessing_utils.py         # Sequence preprocessing utilities
    │   ├── schema_validation.py               # Request validation functions
    │   ├── error_checking_functions.py        # Error checking utilities
    │   ├── predictor_content_handler.py       # Content negotiation handlers
    │   └── predictor_help_message.json        # Help messsage file
```

## Network Configuration

- The container operates as a REST API server and does not require mounted data directories

- All data is transmitted via HTTP requests from evaluators

- The HOST and PORT arguments determine the network interface and port for incoming connections

## API Features

- Supported Formats: JSON and MessagePack for request/response

- Endpoints:

    - GET /help - API documentation and metadata

    - GET /formats - Supported content formats

    - POST /predict - Main prediction endpoint

- Validation: Comprehensive request validation and error handling

## Example Usage
```bash
# Run on specific IP and port
apptainer run --containall legnet_predictor.sif 172.16.47.244 5000

# Run on localhost
apptainer run --containall legnet_predictor.sif 127.0.0.1 5000
```
## Arguments

- HOST: IP address or hostname that the predictor will bind to

    - Use 0.0.0.0 to accept connections from any network interface

    - Use 127.0.0.1 for localhost-only access

- PORT: Port number for the predictor service (typically 5000+)

## Model Specifications

- Input: DNA sequences (200-2000 base pairs)

- Output: Log-scale expression predictions

- Supported Cell Types: K562, HepG2, WTC11

- Species: Homo sapiens only

## Security Notes
- The `--containall` flag provides enhanced container isolation

- No host directory mounts are required for normal operation

- Network access is needed for initial model weight downloads from HuggingFace Hub