
# Adapting and containerizing the CAGI5 'regulation saturation' challenge MPRA data (Kircher et. al (2019) MPRA dataset) for GAME

## Overview

This document outlines the structure of the API codebase for containerizing an MPRA Evaluator using the data from Kircher, M., Xiong, C., Martin, B. et al. Saturation mutagenesis of twenty disease-associated regulatory elements at single base-pair resolution. Nat Commun 10, 3583 (2019). <https://doi.org/10.1038/s41467-019-11526-w>


The Evaluator processes genomic variants from multiple human and mouse cell lines, including: 
- **HepG2 (human hepatocarcinoma)** with promoters (F9, LDLR, LDLR.2) and enhancer (SORT1, SORT1-flip, SORT1.2); 
- **HEL92.1.7 (human erythroleukemia)** with promoters (GP1BA, HBB, HBG1) and enhancer (BCL11A); 
- **HEK293T (human embryonic kidney)** with promoters (HNF4A, MSMB, TERT-HEK) and enhancer (MYCrs6983267);
- **K562 (human myeloid leukemia)** with promoter (PKLR-24h, PKLR-48h);
- **SF7996 (human glioblastoma cells)** with promoter (TERT-GAa, TERT-GBM, TERT-GSc); 
- **HeLa (human cervical adenocarcinoma)** with promoter (FOXE1); 
- **SK-MEL-28 (human melanoma)** with enhancer (IRF4); 
- **HaCaT (human keratinocytes)** with enhancer (IRF6); 
- **LNCaP (human prostate adenocarcinoma)** with enhancer (MYCrs11986220); 
- **Neuro-2a (mouse neuroblastoma cells)** with enhancers (RET, UC88); 
- **MIN6 (mouse insulinoma cells)** with enhancers (TCF7L2, ZFAND3); 
- and **NIH-3T3 (mouse embryonic fibroblasts)** with enhancer (ZRSh-13, ZRSh-13h2).

The data used in this Evaluator were downloaded from https://kircherlab.bihealth.org/satMutMPRA/ and can also be found on [GitHub](https://github.com/visze/sequence_cnn_models/tree/master/resources/satmut_examples). The Evaluator implements a REST API that communicates with Predictor to evaluate variant effects and calculate the correlation metrics between predicted and experimental measurements.

The sequences of promoter and enhancer regions are extracted with 500 kb flanks by default (resulting in a total reference genomic sequence of nearly 1 Mb). For evaluating models which consider longer genomic windows, the preprocessing script `create_evaluator_data.ipynb` can be rerun with the desired parameters.

---

### Using Pre-built Containers

For immediate use, we recommend using the pre-built containers hosted on Hugging Face: https://huggingface.co/datasets/Ni-os/Kircher-container

## 1. Evaluator Structure

The Kircher (2019) Evaluator is organized as follows:

```
/path/to/Kircher_2019_Evaluator/
    ├── evaluator_data/                 # Input data directory (must be mounted)
    ├── predictions/                    # Output directory (must be mounted)
    ├── config.py                       # Configuration file
    ├── data_loader.py                  # Data loading and validation
    ├── evaluator_content_handler.py    # API communication handlers
    ├── evaluator_metrics_calculator.py # Correlation calculations
    ├── kircher_evaluator.def           # This definition file
    ├── kircher_evaluator.py            # Main evaluator script
    └── README.md
```

---

## 2. Understanding the REST API

### 2.1. Evaluator REST API

- **Purpose**: Interfaces with the Predictor REST API to send genomic sequences and receive predictions. It handles input validation, HTTP communication, and output storage.
- **Core Script**: `kircher_evaluator.py`.
- **Key Features**
    1. **Dynamic Path Handling**:
        - Uses `os.path.exists` to determine if the script is running inside the container.
        - Adjusts paths to input (`EVALUATOR_DATA_DIR`) accordingly.
    2. **JSON Payload Creation**:
        - Reads .tsv files from `/evaluator_data` using pandas.
        - Extracts genomic sequences from the reference genome (hg38).
        - Creates a structured JSON payload with prediction tasks and sequences.
        - Validates unique sequence keys to prevent JSON parsing issues.
    3. **REST API Communication**:
        - Establishes an HTTP connection with the Predictor server.
        - Negotiates data formats (JSON/MessagePack) with the predictor.
        - Sends POST requests with genomic data and processes responses.
        - Handles successful (200 OK) and error (4xx/5xx) responses.
    4. **Data Processing**:
        - Estimates the log2 fold changes between the reference and alternative alleles.
        - Computes the Pearson correlation coefficient between the predicted and experimentally measured variant effects for each regulatory element, separately for each element
        - Supports both linear-scale and log-scale model outputs.
    5. **Output Management**:
        - Saves raw predictions and correlation results to mounted directories.
        - Generates comprehensive summary reports with timestamps.

---

### Structure of the API JSON Format

The Evaluator sends JSON payloads with the following structure:
```json
{
    "request": "predict",
    "readout": "point",
    "prediction_tasks": [
        {
            "name": "Kircher_F9",
            "type": "expression", 
            "cell_type": "human hepatocytes (HepG2)",
            "scale": "log",
            "species": "homo_sapiens"
        }
    ],
    "sequences": {
        "variant1_reference": "ATCG...",
        "variant1_alternate": "ATGG...",
        ...
    }
}
```

## 3. Configuring and Running the API Services

### 3.1 Configuring the Containers Using Definition Files

Containerizing the Evaluator involves creating **definition (.def) files** that specify the structure, dependencies, and environment of the container. This ensures a consistent and reproducible environment for running the Predictor and Evaluator services.

### **Purpose of the [Definition Files](https://apptainer.org/docs/user/latest/definition_files.html)**

The definition files provide a declarative way to define:

1. **Base Image**: The starting point of the container (e.g., ubuntu or python:3.13-slim).
2. **File Structure**: Which files to copy into the container during build and which directories to bind/mount at runtime.
3. **Dependencies**: System-level and Python libraries required for the API to function.
4. **Execution Environment**: Configurations such as environment variables, Python environments, and permissions.
5. **Entry point**: Specifies the script to run when the container starts.

### Why Containers for REST APIs?

1. **Reproducibility**: Containers guarantee the same environment across different systems.
2. **Modularity**: Encapsulates all dependencies, scripts, and configurations in a single container.
3. **Portability**: Allows seamless deployment on various platforms, whether it is local or HPC clusters.
4. **HTTP Standardization**: REST APIs provide standardized communication between containers.

### Container Configuration Details

### Evaluator API Container

1. **Base Image**:
    - `python:3.13-slim` is chosen for its lightweight nature and Python-specific optimizations.
    - This base image minimizes overhead while supporting Python dependencies like pandas, scipy, and requests.
2. **File Inclusions**:
    - **Core Script**: All Python modules are copied into `/Kircher_2019_Evaluator/` in the container.
    - **Dynamic Mounting**: Input (`/evaluator_data`) and output (`/predictions`) directories are mounted at runtime to allow flexibility in handling data without modifying the container.
3. **Sequence Extraction Parameters**:
    - Promoter and enhancer regions are extracted with 500 kb flanks by default, resulting in total genomic sequences of approximately 1 Mb.
    - For models requiring longer genomic contexts, use the preprocessing script `create_evaluator_data.ipynb` to regenerate data with custom window sizes.
4. **Dependencies**:
    - Installs required Python packages for genomic data processing and HTTP communication.
    - Includes **pyfastx** for reference genome sequence extraction.
5. **Runtime Needs**:
    - **Mounted Directories**: Binds host directories to `/evaluator_data` (input) and `/predictions` (output).
    - **HTTP Communication**: Connects to Predictor via REST API on specified host and port.

---

### 3.2. Building The Containers

Once the definition files for the Evaluator and Predictor APIs are configured, the next step is to build the container images using the `apptainer build` command.

### Steps to Build the Containers

1. **Set Up the Definition File**:
    - Ensure that the `.def` file includes all necessary files, directories, and dependencies.
    - Files are copied into the container using the `%files` section.
    - Environment-specific variables are listed in the  `%enviroment`  section.
    - Key configurations like `chmod` for full directory and file access are defined in the `%post` section. Dependencies are installed in this section.
    - Execution commands are defined in the `%runscript` section.
    - Metadata and help information are included in `%labels` and `%help` sections.
    - Additional information can be found here: <https://apptainer.org/docs/user/main/definition_files.html>

2. **Run the Build Command**: 

    ```bash
        cd Kircher_2019_Evaluator
        apptainer build kircher_evaluator.sif kircher_evaluator.def
    ```

    - The first argument specifies the output SIF file.
    - The second argument specifies the definition file to use.

3. **Verify Build Success**:

    After the build completes, confirm the presence of the SIF files:

    ```bash
        ls -l *.sif
    ```

    - The resulting `kircher_evaluator.sif` should now be available in the working directory.

### Data Directory Structure

The evaluator_data directory should contain the following structure:

```bash
evaluator_data/
    ├── K562/
    │   ├── PKLR-24h_sequences.tsv
    │   ├── PKLR-24h_values.tsv
    │   ├── PKLR-48h_sequences.tsv
    │   └── PKLR-48h_values.tsv
    ├── HepG2/
    │   ├── F9_sequences.tsv
    │   ├── F9_values.tsv
    │   ├── LDLR_sequences.tsv
    │   ├── LDLR_values.tsv
    │   └──   ...
    └── ...
```

### ⚠️ Permissions Management ⚠️

If you encounter permission errors during container build or execution:

- During container build phase, ensure that the directories and files have correct permissions using `chmod` in the `%post` section:

```bash
    chmod -R 755 /Kircher_2019_Evaluator
```

- This ensures all scripts and resources are accessible inside the container at runtime.

- Prevents issues with restricted permissions that could cause runtime errors.

---

### 3.3. Running the Evaluator API

1. **Start the Predictor First**: Ensure the Predictor container is running and accessible via REST API.

2. **Run the Evaluator Container**:

    ```bash
    apptainer run \
        -B /absolute/path/to/evaluator_data:/evaluator_data \
        -B /absolute/path/to/predictions:/predictions \
        kircher_evaluator.sif PREDICTOR_HOST PREDICTOR_PORT /predictions
    ```

3. **Key arguments**
    - `PREDICTOR_HOST`: IP address or hostname of the Predictor REST API server.
    - `PREDICTOR_PORT`: Port number for the Predictor API (typically 5000)
    - `/predictions`: Container path for output predictions (mounted from host)

    **Example**:

    ```bash
    apptainer run \
        -B /home/user/evaluator_data:/evaluator_data \
        -B /home/user/predictions:/predictions \
        kircher_evaluator.sif 172.17.0.2 5000 /predictions
    ```

4. **Expected Output**:

    The Evaluator will:

    - Process all cell types and corresponding elements;

    - Extract sequences from hg38 reference genome (200 nt long by default);

    - Send sequences to Predictor via REST API;

    - Calculate Pearson correlations between predictions and experimental values;

    - Generate output files in the predictions directory:

        - Raw predictions: `{cell_type}_cell_type_{element}_element_*_predictions.json`

        - Correlation summary: `correlation_summary_Kircher_{timestamp}.tsv`

5. **Verify Successful Execution**:

    Check the predictions directory for output files and review the correlation summary for evaluation results.

6. **Using the Evaluator with Predictors requiring other sequence lengths**

    To allow for compatibility with Predictors which cannot directly handle default sequence lengths of 200bps, the Evaluator can be executed by specifying the custom window size using `--sequence-length` parameter

    - `--sequence-length` or `-l`: The length of extracted sequences. Default is 200.

    **Example**:

    ```bash
    apptainer run \
        -B /home/user/evaluator_data:/evaluator_data \
        -B /home/user/predictions:/predictions \
        kircher_evaluator.sif 172.17.0.2 5000 /predictions --sequence-length 300
    ```


---
