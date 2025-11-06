"""Configuration Script for Evaluator Name, Input File, and Preferred Data Format"""

import os

# --- Core Evaluator Settings ---
# Evaluator name for predictions file and metrics CSV
EVALUATOR_NAME = "Kircher"

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

# --- Directory Settings ---
# Get the absolute path of the script's directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Determine if running inside a container or not
if os.path.exists("/.singularity.d"):
    # Running inside the container
    EVALUATOR_DATA_DIR = "/evaluator_data"
else:
    # Running outside the container
    EVALUATOR_DATA_DIR = os.path.join(SCRIPT_DIR, "evaluator_data")

# --- API Communication Settings ---
REQUEST_FORMAT = "application/json"
REQUEST_FORMAT = REQUEST_FORMAT.lower()

RESPONSE_FORMAT = "application/msgpack"
RESPONSE_FORMAT = RESPONSE_FORMAT.lower()

# HTTP request retry
MAX_RETRIES = 50
RETRY_INTERVAL = 30  # Seconds to wait between each retry attempt
