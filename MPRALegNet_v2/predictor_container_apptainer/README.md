# Configuring the Definition File and Running the Predictor Container for MPRALegNet

## Overview

This container provides a complete MPRALegNet Predictor implementation for regulatory sequence activity prediction. The container includes:

- **REST API Server**: Flask-based API for handling prediction requests
- **MPRALegNet Model**: Pre-trained deep learning model for regulatory sequence activity prediction
- **Sequence Processing**: Comprehensive preprocessing pipeline for DNA sequences
- **Error Handling**: Robust validation and error reporting
- **Dependencies**: Complete Conda environment with all required packages

> **Note**: The Conda environment (`legnet_environment.yml`) excludes PyTorch due to container build conflicts. PyTorch is installed manually in the container definition file.

- Pretrained model weights used in this framework were downloaded from https://zenodo.org/records/10558183

## Usage

### Using Pre-built Containers

For immediate use, we recommend using the pre-built containers hosted on Hugging Face: https://huggingface.co/Ni-os/MPRALegNet-container

### Building from Source

If you need to build the container locally, ensure you have the following directory structure:

```bash
/path/to/MPRALegNet/src/predictor_container_apptainer/
    ├── predictor.def                          # Predictor container definition file
    ├── legnet_predictor.sif                   # Predictor container image (not included on Github due to storage restrictions)
    ├── legnet_environment.yml                 # Conda environment specification
    ├── script_and_utils/
    |   |   ├── hepg2_best_model_test1_val2.ckpt # model weights
    |   |   ├── k562_best_model_test1_val2.ckpt
    |   |   └── wtc11_best_model_test1_val2.ckpt
    │   ├── MPRALegNet_predictor_API.py        # Main Flask API server
    │   ├── MPRALegNet_predict.py              # Model prediction functions
    │   ├── model.py                           # LegNet model architecture
    │   ├── api_preprocessing_utils.py         # Sequence preprocessing utilities
    │   ├── schema_validation.py               # Request validation functions
    │   ├── error_checking_functions.py        # Error checking utilities
    │   ├── predictor_content_handler.py       # Content negotiation handlers
    │   └── predictor_help_message.json        # Help message file
```

### Build the container (SIF)

```bash
apptainer build legnet_predictor.sif predictor.def
```

## Network Configuration

- The container operates as a REST API server and does not require mounted data directories

- All data is transmitted via HTTP requests from evaluators

- The HOST and PORT arguments determine the network interface and port for incoming connections

## API Features

- **Supported Formats**: JSON and MessagePack for request/response

- **Endpoints**:

    - GET /help - API documentation and metadata

    - GET /formats - Supported content formats

    - POST /predict - Main prediction endpoint

- **Validation**: Comprehensive request validation and error handling

## Example Usage
```bash
# Run on specific IP and port
apptainer run --containall legnet_predictor.sif 172.16.47.244 5000

# Run on localhost
apptainer run --containall legnet_predictor.sif 127.0.0.1 5000
```
## Arguments

- **HOST**: IP address or hostname that the predictor will bind to

    - Use 0.0.0.0 to accept connections from any network interface

    - Use 127.0.0.1 for localhost-only access

- **PORT**: Port number for the predictor service (typically 5000+)

## Model Specifications

- **Input**: DNA sequences (200-2000 base pairs, 200 by default).

- **Output**: Regulatory sequence activity predictions.

- **Scale**: MPRALegNet is trained to predict log2(RNA/DNA) values, thus its 'native' predictions are in log2-scale. The conversion, if requested by an evaluator, is performed in the following way: 2^prediction.

- **Supported Cell Types**: K562, HepG2, WTC11.

- **Species**: Homo sapiens only.

- **API Formats**: JSON and MessagePack support.

- **Preprocessing**: Automatic adapter handling, sequence padding, one-hot encoding

## Security Notes

- The `--containall` flag provides enhanced container isolation

- No host directory mounts are required for normal operation