# Adapting and Containerizing the MPRALegNet Model for Genomic Model Evaluation API

## **Overview**

This repository contains both the implementation of the MPRALegNet CNN model for predicting regulatory sequence activity from DNA sequences and a dummy evaluator to ensure the model works correctly and its performance can be evaluated within the Genomic API for Model Evaluation (GAME) ecosystem.

The framework provides:
- **MPRALegNet Predictor**: A containerized deep learning model that predicts regulatory sequence activity from DNA sequences
- **Test Evaluator**: A containerized testing framework that formally validates the technical predictor functionality and reports fake performance metrics.
- **REST API Integration**: Standardized API endpoints for seamless communication between components.

---

---

## **MPRALegNet-API  Structure**

The MPRALegNet-API is organized as follows:

```bash
MPRALegNet/
    ├── evaluator_container_apptainer                    # Test evaluator for validating predictor
    │   ├── evaluator.def                                # Evaluator container definition file
    │   ├── evaluator.sif                                # Evaluator container image (not included on Github due to storage restrictions)
    │   ├── MPRALegNet_evaluator_API.py                  # Evaluator API script
    │   ├── data_loader.py                               # Data loading and validation
    |   ├── config.py                                    # Configuration file
    |   ├── evaluator_content_handler.py                 # API communication handlers
    |   ├── evaluator_metrics_calculator.py              # Correlation and specificity calculations
    │   ├── evaluator_data                
    │   │   ├── evaluator_input_sample_test.json         # Sample input JSON for evaluator
    │   │   ├── evaluator_message_gosai_5seqs.json       # Gosai 5 test sequences
    │   │   ├── evaluator_message_more_complex.json      # Complex evaluator input example
    │   │   └── evaluator_message_simple_test.json       # Simple evaluator input example
    │   └── predictions/                                 # Output directory for prediction results
    |
    └── predictor_container_apptainer              # MPRALegNet model implementation
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

---

---

# Quick Start: Running MPRALegNet with the Demo Data

## Prerequisites

- Apptainer installed in the host OS

- Basic familiarity with command line operations

## Step 1: Download Containers and Data

Create a working directory and download the required containers and test data:

```bash
# Create main working directory
mkdir MPRALegNet
cd MPRALegNet

# Download predictor container from Hugging Face
wget -O legnet_predictor.sif "https://huggingface.co/Ni-os/MPRALegNet-container/resolve/main/legnet_predictor.sif?download=true"

# Create test evaluator directory
mkdir test_evaluator
cd test_evaluator

# Download evaluator container Hugging Face
wget -O legnet_evaluator.sif "https://huggingface.co/Ni-os/MPRALegNet-container/resolve/main/legnet_evaluator.sif?download=true"

# Download test dataset
wget -O evaluator_data.zip "https://huggingface.co/Ni-os/MPRALegNet-container/resolve/main/evaluator_data.zip?download=true"
unzip evaluator_data.zip

# Create output directory for predictions
mkdir predictions
```

## Step 2: Determine Network Configuration

Identify the IP address where the predictor will run:

```bash
hostname -I
```

Note: Choose an available port above 5000 (e.g., 5000, 5001, 8080)

## Step 3: Start the MPRALegNet Predictor

Launch the predictor container with your specific IP and port:

```bash
# Replace with your actual IP and chosen port
apptainer run legnet_predictor.sif 172.16.47.xxx 5000
```

### Expected Output:

```bash
MPRALegNet Predictor is running on http://172.16.47.xxx:5000
```

The predictor will now be listening for API requests on the specified address.

## Step 4: Run the Test Evaluator
 
In a new terminal window, navigate to your test directory and start the evaluator:

```bash
cd MPRALegNet/test_evaluator

apptainer run --containall \
    -B /absolute/path/to/evaluator_data:/evaluator_data \
    -B /absolute/path/to/predictions:/predictions \
    legnet_evaluator.sif 172.16.47.xxx 5000 /predictions
```

## Step 5: Verify Successful Execution

After the evaluator completes, check the predictions directory for results:

```bash
ls -la predictions/
```

### Expected Output Files:

`fake_correlation_summary_MPRALegNet_evaluator.tsv` - Performance metrics

`*_predictions_from_MPRALegNet.json` - Raw prediction results

---
---

# Component Details

## MPRALegNet Predictor Features

- **Model Architecture**: LegNet deep learning model optimized for genomic sequences

- **Supported Cell Types**: K562 (lymphoblasts), HepG2 (hepatocytes), WTC11 (stem cells)

- **Input**: DNA sequences (200-2000 base pairs)

- **Output**: Log-scale regulatory sequence activity predictions

- **API Formats**: JSON and MessagePack support

- **Preprocessing**: Automatic adapter handling, sequence padding, one-hot encoding

## Test Evaluator Features

- **Validation Suite**: Comprehensive API endpoint testing

- **Fake Metrics**: Pearson correlation, specificity calculations

- **Sample Datasets**: Multiple test scenarios with varying complexity

- **Standardized Output**: Consistent result formatting for comparison

## API Endpoints

**Predictor Endpoints**:

- `GET /help` - API documentation and model metadata

- `GET /formats` - Supported request/response formats

- `POST /predict` - Main prediction endpoint

**Evaluator Endpoints**:

- Standardized GAME API for model evaluation

- Fake performance metric calculation

---
---

# Advanced Usage

## Building Containers from Source

If you need to rebuild the containers:

```bash
# Build predictor container
cd predictor_container_apptainer
apptainer build legnet_predictor.sif predictor.def

# Build evaluator container  
cd ../evaluator_container_apptainer
apptainer build legnet_evaluator.sif evaluator.def
```

## GPU Acceleration

For improved performance, run the predictor with GPU support:

```bash
apptainer run --nv legnet_predictor.sif 172.16.47.xxx 5000
```

---
---

# Troubleshooting
## Common Issues

1. Connection Refused

    - Verify whether the Predictor is running before starting the Evaluator

    - Check firewall settings and port availability

2. Container Permission Errors

    - Ensure you have sufficient privileges to run containers

    - Use --containall flag for enhanced security

3. Missing Dependencies

    - Containers include all dependencies, but verify Apptainer installation

## Getting Help

- Check container logs for detailed error messages

- Verify all file paths are absolute in bind mounts

- Ensure sufficient disk space for temporary files

---
---

# Conclusion

🎉 Success! You have successfully completed a full genomic model evaluation cycle using MPRALegNet within the GAME framework. The framework demonstrates:

- Containerized deep learning model deployment

- Standardized API communication between components

- Automated fake performance evaluation and metric calculation

- Reproducible research workflows

For more advanced usage, refer to the individual component documentation and explore the sample datasets to understand the full capabilities of the MPRALegNet genomic prediction system.