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
if os.path.exists('/.singularity.d'):
    # Running inside the container
    BASE_DIR = '/script_and_utils'
    UTILS_DIR = '/script_and_utils'
else:
    # Running outside the container
    # Get the current working directory
    BASE_DIR = SCRIPT_DIR
    UTILS_DIR = SCRIPT_DIR

# Add the directory containing `error_checking_functions.py` to the Python path
sys.path.append(UTILS_DIR)

# Import the process_sequence function
from api_preprocessing_utils import process_sequence

# Import model
from model import LegNet

# Define the model directory
MODEL_DIR = os.path.join(BASE_DIR, 'weights')

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


def load_model(cell_type="k562", device="cpu"):
    """
    Load pre-trained MPRALegNet model weights.
    
    This function handles the complete model loading pipeline:
    1. Initializes the LegNet architecture with correct parameters
    2. Loads cell-type specific pre-trained weights
    3. Sets model to evaluation mode for inference
    
    Args:
        cell_type (str): Target cell type for prediction. Must be one of:
                        "k562", "hepg2", or "wtc11"
        device (torch.device): Device to load the model onto (CPU or CUDA)
    
    Returns:
        LegNet: Loaded and configured MPRALegNet model in evaluation mode
    """

    # Initialize LegNet model with architecture parameters
    stem_ks = 5 if cell_type.lower() == "wtc11" else 11
    ef_ks = 3 if cell_type.lower() == "wtc11" else 9

    model = LegNet(
        in_ch=4,
        stem_ch=64,
        stem_ks=stem_ks, 
        ef_ks=ef_ks,
        ef_block_sizes=[80, 96, 112, 128],
        pool_sizes=[2, 2, 2, 2] ,
        resize_factor=4,
        activation=torch.nn.SiLU,
    ).to(device)

    # Map cell types to their corresponding weight files
    weight_files = {
        "hepg2": "hepg2_best_model_test1_val2.ckpt",
        "k562": "k562_best_model_test1_val2.ckpt",
        "wtc11": "wtc11_best_model_test1_val2.ckpt",
    }

    # Download cell-type specific pre-trained weights
    weights_path = weight_files[cell_type.lower()]
    model_file_path = os.path.join(MODEL_DIR, weights_path)

    # Load weights into model architecture
    checkpoint = torch.load(model_file_path, map_location=device, weights_only=False)

    # Get the state_dict from checkpoint
    state_dict = checkpoint['state_dict']

    # Remove prefixes that PyTorch Lightning adds
    new_state_dict = {}
    for key, value in state_dict.items():
        # Remove common PyTorch Lightning prefixes
        if key.startswith('model.'):
            new_key = key[6:]  # Remove 'model.'
        elif key.startswith('net.'):
            new_key = key[4:]  # Remove 'net.'
        elif key.startswith('module.'):
            new_key = key[7:]  # Remove 'module.' (for DataParallel)
        else:
            new_key = key
        
        new_state_dict[new_key] = value
    

    # Load the modified state_dict into model
    model.load_state_dict(new_state_dict)

    # Set model to evaluation mode for inference (disables dropout, batch norm in eval mode)
    model.eval()
    
    print(f"Model for {cell_type} loaded!")
    return model


# Hardcoded adapter sequences required by MPRALegNet model architecture
# These 15bp sequences are added to each input sequence for proper model recognition
upstream_adapter_seq = "AGGACCGGATCAACT"
downstream_adapter_seq = "CATTGCGTGAACCGA"


# Updated Prediction Function
def predict_legnet(sequences, prediction_ranges, upstream_seq, downstream_seq, cell_type, scale, window_step=50):
    """
    Predict expression values using the MPRALegNet model for genomic sequences.
    
    This is the main prediction function that processes sequences through the full pipeline:
    1. Sequence preprocessing (adapter handling, padding, encoding)
    2. Model inference on both original and reverse complement sequences
    3. Aggregation of predictions from multiple windows for long sequences
    4. Averaging of predictions from both strands for robust results
    
    Args:
        sequences (dict): Dictionary mapping sequence IDs to DNA sequence strings
        prediction_ranges (dict): A collection of key-value pairs, where the keys should be identical to sequence ID keys 
                                    and values are arrays with the start and end region you want predictions for, within the provided sequence context. 
                                    Start and end are 0 indexed and inclusive, respectively (e.g. [0,1] is the first two bases).
        upstream_seq (str): Upstream flanking sequences to add to each sequence in sequences.
        downstream_seq (str): Downstream flanking sequences to add to each sequence in sequences.
        cell_type (str): Target cell type for prediction ("k562", "hepg2", or "wtc11")
        scale (str): Requested scale. 'log' or 'linear'
        window_step (int): Step size for sliding window processing of long sequences (>200bp)
    
    Returns:
        dict: Dictionary mapping sequence IDs to predicted expression values.
              For single sequences: single value list
              For multi-window sequences: averaged predictions across windows
    """
    # Detect and set computation device (GPU if available, otherwise CPU)
    device = get_device()

    print("Loading pre-trained model weights for MPRALegNet")
    model_legnet = load_model(cell_type, device=device)

    predictions = {}

    # Process each sequence with progress tracking
    for seq_id, seq in tqdm.tqdm(
        sequences.items(), desc="Predictions in progress", unit="sequence"
    ):
        
        seq_specific_prediction_ranges = None
        if prediction_ranges is not None:
                seq_specific_prediction_ranges = prediction_ranges.get(seq_id) or None    

        len_upstream_seq = len(upstream_seq) if upstream_seq is not None else 0
        len_downstream_seq = len(downstream_seq) if downstream_seq is not None else 0
        seq_size = len_upstream_seq + len(seq) + len_downstream_seq # Final sequence length with adapters
        
        
        # Process sequence through full preprocessing pipeline
        # Returns list of (encoded_seq, reversed_encoded_seq) tuples
        # - For short sequences: single tuple in list
        # - For long sequences: multiple tuples for each window
        processed_list = process_sequence(
            seq, 
            seq_specific_prediction_ranges,
            upstream_seq=upstream_seq,
            downstream_seq=downstream_seq,
            target_length=200,      # Core sequence length without adapters
            seq_size=seq_size,           
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
            
            # Convert to linear scale if requested
            if scale == 'linear':
                # Since model was trained on log2(RNA/DNA), convert back using 2^prediction
                avg_prediction = [2**x for x in avg_prediction]
            
            predictions[seq_id] = avg_prediction
        else:
            # Fallback for edge case with no predictions
            fallback_value = 0.0 if scale == 'log' else 1.0  # 2^0 = 1 in linear scale
            predictions[seq_id] = [fallback_value]

    return predictions
