"""Configuration Script for Evaluator Name, Input File, and Preferred Data Format"""

import os
import json
from datetime import datetime

# --- Core Evaluator Settings ---
# Evaluator name for predictions file and metrics CSV
evaluator_base = "Kircher"

# Dictionary mapping cell types to their corresponding genetic elements
# Used for organizing and processing data specific to each cell line
ELEMENTS_BY_CELL_TYPE = {
    "HepG2": ["F9", "LDLR", "LDLR.2", "SORT1", "SORT1-flip", "SORT1.2"],
    "HEL92.1.7": ["GP1BA", "HBB", "HBG1", "BCL11A"],
    "HEK293T": ["HNF4A", "MSMB", "MYCrs6983267", "TERT-HEK"],
    "K562": ["PKLR-24h", "PKLR-48h"],
    "SF7996": ["TERT-GAa", "TERT-GBM", "TERT-GSc"],
    "HeLa": ["FOXE1"],
    "SK-MEL-28": ["IRF4"],
    "HaCaT": ["IRF6"],
    "LNCaP": ["MYCrs11986220"],
    "Neuro-2a": ["RET", "UC88"],
    "MIN6": ["TCF7L2", "ZFAND3"],
    "NIH-3T3": ["ZRSh-13", "ZRSh-13h2"],
}

# Dictionary mapping cell type codes to their full biological names
# Contains standardized nomenclature for different cell lines used in the study
CELL_TYPES_NAMES = {
    "HepG2": "human hepatocytes (HepG2)",
    "HEL92.1.7": "human erythroleukemia cells (HEL92.1.7)",
    "HEK293T": "human embryonic kidney cells (HEK293T)",
    "K562": "lymphoblasts (K562)",
    "SF7996": "glioblastoma cells (SF7996)",
    "HeLa": "cervical adenocarcinoma cells (HeLa)",
    "SK-MEL-28": "melanoma cells (SK-MEL-28)",
    "HaCaT": "keratinocytes (HaCaT)",
    "LNCaP": "prostate adenocarcinoma cells (LNCaP)",
    "Neuro-2a": "mouse neuroblastoma cells (Neuro-2a)",
    "MIN6": "mouse insulinoma cells (MIN6)",
    "NIH-3T3": "mouse embryonic fibroblasts (NIH-3T3)",
}

CELL_TYPES_ORGANISM = {
    "HepG2": "homo_sapiens",
    "HEL92.1.7": "homo_sapiens",
    "HEK293T": "homo_sapiens",
    "K562": "homo_sapiens",
    "SF7996": "homo_sapiens",
    "HeLa": "homo_sapiens",
    "SK-MEL-28": "homo_sapiens",
    "HaCaT": "homo_sapiens",
    "LNCaP": "homo_sapiens",
    "Neuro-2a": "mus_musculus",
    "MIN6": "mus_musculus",
    "NIH-3T3": "mus_musculus",
}

# --- Directory Settings ---
# Get the absolute path of the script's directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Determine if running inside a container or not
if os.path.exists("/.singularity.d"):
    # Running inside the container
    EVALUATOR_DATA_DIR = "/evaluator_data"
    
    # Read build timestamp from Apptainer's auto-generated labels file to version the evaluator name
    try:
        with open('/.singularity.d/labels.json', 'r') as f:
            labels = json.load(f)
        raw_build_date = labels.get('org.label-schema.build-date', '')
        
        # Example format: "Friday_28_November_2025_18:6:29_PST"
        # Strip day-of-week and timezone, keep the core date+time
        parts = raw_build_date.split('_')
        # parts: ['Friday', '28', 'November', '2025', '18:6:29', 'PST']
        date_str = f"{parts[1]}_{parts[2]}_{parts[3]}_{parts[4]}"
        
        dt = datetime.strptime(date_str, "%d_%B_%Y_%H:%M:%S")
        build_timestamp = dt.strftime("%Y%m%d-%H%M%S")
        timezone_label = parts[5] if len(parts) > 5 else "UNK"
        EVALUATOR_NAME = f"{evaluator_base}_{build_timestamp}_{timezone_label}"
    
    except Exception as e:
        print(f"Warning: Could not parse build timestamp from labels.json: {e}")
        EVALUATOR_NAME = f"{evaluator_base}_unknown"
        
else:
    # Running outside the container
    EVALUATOR_DATA_DIR = os.path.join(SCRIPT_DIR, "evaluator_data")
    EVALUATOR_NAME = f"{evaluator_base}_dev"

# --- API Communication Settings ---
REQUEST_FORMAT = "application/json"
REQUEST_FORMAT = REQUEST_FORMAT.lower()

RESPONSE_FORMAT = "application/msgpack"
RESPONSE_FORMAT = RESPONSE_FORMAT.lower()

# HTTP request retry
MAX_RETRIES = 50
RETRY_INTERVAL = 30  # Seconds to wait between each retry attempt
