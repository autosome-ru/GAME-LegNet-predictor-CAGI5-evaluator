# MPRALegNet_predict.py
#!/usr/bin/env python
# coding: utf-8
import os
import sys
import tqdm
import torch
import numpy as np

# Get the absolute path of the script directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Determine if running inside a container or not
if os.path.exists("/.singularity.d"):
    # Running inside the container - use fixed container path
    UTILS_DIR = "/predictor_container_apptainer"
else:
    # Running outside the container - use relative path to utils
    UTILS_DIR = os.path.join(SCRIPT_DIR, "..", "script_and_utils")

BASE_DIR = SCRIPT_DIR

# Add the directory containing `api_preprocessing_utils.py` to the Python path
sys.path.append(UTILS_DIR)

# Import the process_sequence function
from api_preprocessing_utils import process_sequence

# Import functions to download model from HuggingFace
from huggingface_hub import hf_hub_download
from safetensors.torch import load_file
import json
from model import LegNet


def get_device():
    """
    Automatically detects and returns the best available computing device.
    
    Checks for CUDA-enabled GPU availability and falls back to CPU if no GPU is available.
    This ensures optimal performance while maintaining compatibility.
    
    Returns:
        torch.device: CUDA device if available, otherwise CPU device
    """
    if torch.cuda.is_available():
        return torch.device("cuda")
    else:
        return torch.device("cpu")


def download_and_load_model(cell_type="k562", repo_id="Ni-os/MPRALegNet", device="cpu"):
    """
    Download and load pre-trained MPRALegNet model weights from HuggingFace Hub.
    
    This function handles the complete model loading pipeline:
    1. Downloads model configuration and weights from HuggingFace
    2. Initializes the LegNet architecture with correct parameters
    3. Loads cell-type specific pre-trained weights
    4. Sets model to evaluation mode for inference
    
    Args:
        cell_type (str): Target cell type for prediction. Must be one of:
                        "k562", "hepg2", or "wtc11"
        repo_id (str): HuggingFace repository ID containing the model
        device (torch.device): Device to load the model onto (CPU or CUDA)
    
    Returns:
        LegNet: Loaded and configured MPRALegNet model in evaluation mode
    """
    # Download main config file from HuggingFace Hub
    config_path = hf_hub_download(repo_id=repo_id, filename="config.json")

    # Load model configuration parameters
    with open(config_path, "r") as f:
        config = json.load(f)

    # Initialize LegNet model with architecture parameters from config
    model = LegNet(
        in_ch=config["in_ch"],
        stem_ch=config["stem_ch"],
        stem_ks=config["stem_ks"],
        ef_ks=config["ef_ks"],
        ef_block_sizes=config["ef_block_sizes"],
        pool_sizes=config["pool_sizes"],
        resize_factor=config["resize_factor"],
        activation=torch.nn.SiLU,
    ).to(device)

    # Map cell types to their corresponding weight files
    weight_files = {
        "hepg2": "weights/hepg2_best_model_test1_val2.safetensors",
        "k562": "weights/k562_best_model_test1_val2.safetensors",
        "wtc11": "weights/wtc11_best_model_test1_val2.safetensors",
    }

    # Download cell-type specific pre-trained weights
    weights_path = hf_hub_download(
        repo_id=repo_id, filename=weight_files[cell_type.lower()]
    )

    # Load weights into model architecture
    state_dict = load_file(weights_path)
    model.load_state_dict(state_dict)

    # Set model to evaluation mode for inference (disables dropout, batch norm in eval mode)
    model.eval()
    print(f"Model for {cell_type} loaded!")
    return model


# Hardcoded adapter sequences required by MPRALegNet model architecture
# These 15bp sequences are added to each input sequence for proper model recognition
upstream_adapter_seq = "AGGACCGGATCAACT"
downstream_adapter_seq = "CATTGCGTGAACCGA"


# Updated Prediction Function
def predict_legnet(sequences, cell_type, window_step=50):
    """
    Predict expression values using the MPRALegNet model for genomic sequences.
    
    This is the main prediction function that processes sequences through the full pipeline:
    1. Sequence preprocessing (adapter handling, padding, encoding)
    2. Model inference on both original and reverse complement sequences
    3. Aggregation of predictions from multiple windows for long sequences
    4. Averaging of predictions from both strands for robust results
    
    Args:
        sequences (dict): Dictionary mapping sequence IDs to DNA sequence strings
        cell_type (str): Target cell type for prediction ("k562", "hepg2", or "wtc11")
        window_step (int): Step size for sliding window processing of long sequences (>200bp)
    
    Returns:
        dict: Dictionary mapping sequence IDs to predicted expression values.
              For single sequences: single value list
              For multi-window sequences: averaged predictions across windows
    """
    # Detect and set computation device (GPU if available, otherwise CPU)
    device = get_device()

    print("Loading pre-trained model weights for Human Legnet")
    model_legnet = download_and_load_model(cell_type, device=device)

    predictions = {}

    # Process each sequence with progress tracking
    for seq_id, seq in tqdm.tqdm(
        sequences.items(), desc="Predictions in progress", unit="sequence"
    ):
        # Process sequence through full preprocessing pipeline
        # Returns list of (encoded_seq, reversed_encoded_seq) tuples
        # - For short sequences: single tuple in list
        # - For long sequences: multiple tuples for each window
        processed_list = process_sequence(
            seq, 
            target_length=200,      # Core sequence length without adapters
            seq_size=230,           # Final sequence length with adapters (200 + 15 + 15)
            upstream_adapter=upstream_adapter_seq,
            downstream_adapter=downstream_adapter_seq, 
            window_step=window_step
        )

        window_predictions = []

        # Process each window (or single sequence) and its reverse complement
        for encoded_seq, reversed_encoded_seq in processed_list:
            # Convert encoded sequence to tensor format expected by MPRALegNet
            # Reshape to (batch_size, channels, sequence_length) format
            seq_tensor = torch.tensor(
                np.array(encoded_seq)
                .reshape(1, len(encoded_seq), 4) # Reshape to (1, seq_len, 4 channels)
                .transpose(0, 2, 1),             # Transpose to (1, 4, seq_len)
                device=device,
                dtype=torch.float32,
            )

            # Predict with original sequence orientation
            pred = model_legnet(seq_tensor)

            # Convert reverse complement encoded sequence to tensor
            reversed_seq_tensor = torch.tensor(
                np.array(reversed_encoded_seq)
                .reshape(1, len(reversed_encoded_seq), 4)
                .transpose(0, 2, 1),
                device=device,
                dtype=torch.float32,
            )

            # Predict with reverse complement sequence
            reversed_pred = model_legnet(reversed_seq_tensor)

            # Average predictions from both strands for strand-independent prediction
            result = (pred + reversed_pred) / 2
            window_predictions.append(result.detach().cpu().flatten().tolist())

        # Aggregate predictions from all windows (for long sequences)
        if window_predictions:
            # Average predictions across all windows to get final sequence prediction
            avg_prediction = np.mean(window_predictions, axis=0).tolist()
            predictions[seq_id] = avg_prediction
        else:
            # Fallback for edge case with no predictions
            predictions[seq_id] = [0.0]

    return predictions
