# api_preprocessing_utils.py
import os
import sys

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

from error_checking_functions import BadRequestError

# Adapted from DREAM_RNN predictor
# https://github.com/de-Boer-Lab/GAME-DREAM-RNN-human-k562-predictor/blob/main/src/dream_rnn_script_and_utils/dream_rnn_preprocessing_utils.py#L82
# Assay-specific offset: in Agarwal et al. 2025's lentiMPRA constructs, on which
# DREAM-RNN model (and MPRALegNet model) is trained on, the regulatory probe ends 145 bp upstream
# of the EGFP TSS. 
# Used in Branch 2.1 of process_sequence to anchor the model's input window
# relative to the TSS when only the raw genomic sequence + prediction_ranges
# (no flanks) are provided.
MPRA_PROBE_TO_TSS_OFFSET = 145

def contains_adapters(seq, upstream_adapter, downstream_adapter):
    """
    Check if sequence already contains adapter sequences.

    Verifies whether the input sequence starts with the upstream adapter
    and ends with the downstream adapter, indicating it's already processed.

    Args:
        seq (str): The input sequence to check
        upstream_adapter (str): Upstream adapter sequence expected at start
        downstream_adapter (str): Downstream adapter sequence expected at end

    Returns:
        bool: True if both adapters are present at their respective positions
    """
    return seq.startswith(upstream_adapter) and seq.endswith(downstream_adapter)


def remove_adapters(seq, upstream_adapter, downstream_adapter):
    """
    Remove adapter sequences if present.

    If the sequence contains both upstream and downstream adapters,
    removes them to get the core sequence. Otherwise returns original sequence.

    Args:
        seq (str): The input sequence that may contain adapters
        upstream_adapter (str): Upstream adapter sequence to remove from start
        downstream_adapter (str): Downstream adapter sequence to remove from end

    Returns:
        str: Core sequence without adapter sequences
    """
    if contains_adapters(seq, upstream_adapter, downstream_adapter):
        return seq[len(upstream_adapter) : -len(downstream_adapter)]
    return seq


def pad_sequence(seq, target_length):
    """
    Pad a sequence with 'N' nucleotides to reach target length.

    For sequences shorter than the target length, adds 'N' padding at the beginning.
    For sequences longer than target length, truncates to the first target_length bases.
    MPRALegNet requires consistent input sequence length for model processing.

    Args:
        seq (str): The input sequence (without adapters) to pad
        target_length (int): Model-dependent desired length of the core sequence
                             (200 for MPRALegNet)

    Returns:
        str: Padded or truncated sequence of exactly target_length
    """
    seq_len = len(seq)

    if seq_len < target_length:
        # Add 'N' padding to the beginning of short sequences
        total_padding = target_length - seq_len
        pad = "N" * total_padding
        padded_seq = pad + seq
        return padded_seq

    elif seq_len == target_length:
        # Sequence is already correct length
        return seq

    else:
        # For sequences longer than target_length, take the first target_length bases
        truncated_seq = seq[:target_length]
        return truncated_seq

def _pad_with_flanks(seq, upstream_seq, downstream_seq, target_length):
    
    """
    Pads a probe (variable region) shorter than target_length with flank
    context. Takes the tail of the upstream flank and the head of the downstream flank
    so the probe sits in the same context it would in the actual construct.
    
    If the available flanks are not long enough to fill target_length, the
    remaining is filled with N-padding.
    
    Args:
        seq (str): The variable region sequence (probe)
        upstream_seq (str): The upstream flank sequence provided by the evaluator
        downstream_seq (str): The downstream flank sequence provided by the evaluator
        target_length (int): The desired length of the output sequence (200 bp)
        
    Returns:
        padded_seq (str): The sequence with flanks applied and padded to target_length if needed
    """
    
    if len(seq) >= target_length:
        return seq
    
    # Calculate how much padding is needed
    total_padding = target_length - len(seq)
    upstream_pad_len = total_padding // 2
    downstream_pad_len = total_padding - upstream_pad_len
    
    # Take the TAIL of the upstream_seq and the HEAD of the downstream_seq
    upstream_pad = upstream_seq[-upstream_pad_len:] if upstream_pad_len > 0 else ""
    downstream_pad = downstream_seq[:downstream_pad_len] if downstream_pad_len > 0 else ""
    
    padded_seq = upstream_pad + seq + downstream_pad
    
    # If flanked sequence is still shorter than target_length, pad with Ns
    if len(padded_seq) < target_length:
        padded_seq = pad_sequence(padded_seq, target_length)
    
    return padded_seq


def process_long_sequence(seq, target_length, window_step=50):
    """
    Process long sequence by extracting overlapping windows.

    For sequences longer than the MPRALegNet target length, creates multiple
    overlapping windows using a sliding window approach. This allows the model
    to process different regions of long sequences.

    Args:
        seq (str): The input sequence (without adapters)
        target_length (int): Desired length for each window (200 for MPRALegNet)
        window_step (int): Step size for sliding window extraction

    Returns:
        list: List of window sequences, each of target_length
    """
    seq_len = len(seq)
    windows = []

    # If sequence is already within target length, return as single window
    if seq_len <= target_length:
        return [seq]

    # Extract overlapping windows using sliding window approach
    for start in range(0, seq_len - target_length + 1, window_step):
        end = start + target_length
        window = seq[start:end]
        windows.append(window)

    # Ensure we capture the end of the sequence if the window step missed it
    if not windows or (seq_len - target_length) % window_step != 0:
        last_window = seq[-target_length:]
        if last_window not in windows:  # Avoid duplicates
            windows.append(last_window)

    return windows


def add_adapters(seq, upstream_adapter, downstream_adapter):
    """
    Add upstream and downstream adapter sequences to a given sequence.

    Appends MPRALegNet-specific adapter sequences to the core sequence.
    This is required for proper model recognition and processing.

    Args:
        seq (str): The input sequence after padding, if needed
        upstream_adapter (str): Upstream adapter sequence to prepend
        downstream_adapter (str): Downstream adapter sequence to append

    Returns:
        str: The sequence with adapters added at both ends
    """
    seq_with_adapters = upstream_adapter + seq + downstream_adapter
    return seq_with_adapters


# One-hot encode sequences
def one_hot_encode(seq):
    """
    Convert a DNA sequence to a one-hot encoded representation.

    Maps each nucleotide to a 4-element binary vector:
    A -> [1,0,0,0], G -> [0,1,0,0], C -> [0,0,1,0], T -> [0,0,0,1], N -> [0,0,0,0]

    Args:
        seq (str): DNA sequence consisting of 'A', 'T', 'G', 'C', and 'N'

    Returns:
        list: A list of 4-element lists representing one-hot encoded sequence
    """
    mapping = {
        "A": [1, 0, 0, 0],
        "G": [0, 1, 0, 0],
        "C": [0, 0, 1, 0],
        "T": [0, 0, 0, 1],
        "N": [0, 0, 0, 0],
    }  # 'N' represents ambiguous base or padding
    return [mapping.get(base, [0, 0, 0, 0]) for base in seq]


def reverse_complement(seq, mapping={"A": "T", "G": "C", "T": "A", "C": "G", "N": "N"}):
    """
    Generate the reverse complement of a DNA sequence.

    MPRALegNet uses both original and reverse complement sequences for
    strand-independent prediction, improving model robustness.

    Args:
        seq (str): Input DNA sequence
        mapping (dict): Nucleotide complement mapping dictionary

    Returns:
        str: Reverse complement of the input sequence
    """
    s = "".join(mapping[s] for s in reversed(seq))
    return s

# Validate if prediction ranges are within flanked sequence bounds (adapted from DREAM-RNN predictor)
# This could be different for different Predictors, so kept separate
# https://github.com/de-Boer-Lab/GAME-DREAM-RNN-human-k562-predictor/blob/main/src/dream_rnn_script_and_utils/dream_rnn_preprocessing_utils.py#L82
def _validate_prediction_range_bounds(prediction_ranges, seq_length):
    """
    Validate if prediction ranges are within the bounds of the flanked sequence length.
    Called only when prediction_ranges are actually used --
    CASE 2.1 (in process_sequence) no flanks provided + prediction_ranges provided.

    Args:
        prediction_ranges (list): List containing start and end indices of prediction ranges
        seq_length (int): Length of the flanked sequence

    Raises:
        BadRequestError: If prediction ranges are out of bounds.
    """
    if prediction_ranges is None:
        return

    pr_start, pr_end = prediction_ranges
    
    if pr_start < 0 or pr_end > seq_length:
        raise BadRequestError(
            f"prediction_ranges: [{pr_start}, {pr_end}] are out of bounds for flanked sequence length: {seq_length}."
        )

# Adapted from DREAM-RNN predictor
# https://github.com/de-Boer-Lab/GAME-DREAM-RNN-human-k562-predictor/blob/main/src/dream_rnn_script_and_utils/dream_rnn_preprocessing_utils.py#L82
def _crop_upstream_of_tss(seq, pr_start, target_length, offset):
    
    """
    Returns a window of the sequence upstream of the TSS, anchored by 
    the prediction range start and a known offset, for an MPRA-style assay.
    Ideal window: target_length bp ending `offset` bp before pr_start (the TSS).
    
    When the sequence is long enough, this returns exactly that window:
        [window_start ...... window_end] -- offset bp (145 bp) -- [pr_start = TSS]
        |------- target_length --------|
    
    When pr_start is too close to the start of the sequence to fit the
    full target_length window upstream of TSS - offset, we clamp the start
    to position 0 and pull a full target_length window from there. This
    pulls maximum real sequence context rather than inserting unnecessary
    Ns just to preserve the "biologically correct" offset to the TSS.
    
    This may return a window that is shorter than target_length if the 
    sequence itself is too short, but in that case we will pad with Ns 
    in process_sequence and log a warning.
    
    Args:
        seq (str): The input sequence
        pr_start (int): The start index of the prediction range (interpreted as TSS)
        target_length (int): The desired length of the output window
        offset (int): The known distance from the prediction range start to the TSS,
                      145 bp in this case.
                      
    Returns:
        window (str): The sequence window upstream of the TSS, ideally 
                      of length target_length
    """
    
    # Clamp to actual sequence bounds
    start = max(0, pr_start - offset - target_length)
    end = min(len(seq), start + target_length)
    window = seq[start:end]
    
    if len(window) < target_length:
        # If we don't have enough sequence upstream of the TSS, 
        # we can pad with Ns in process_sequence
        sys.stderr.write(
            f"LOG: TSS-anchored window pulled only {len(window)} bp of real "
            f"upstream sequence (wanted {target_length}); shortfall will be "
            f"left-padded with Ns. Predictions may be unreliable.\n"
        )
    
    return window

# Full preprocessing pipeline for a sequence
def process_sequence(
    seq, prediction_ranges, upstream_seq, downstream_seq, target_length, seq_size, upstream_adapter, downstream_adapter, window_step=50
):
    """ "
    Complete preprocessing pipeline for MPRALegNet sequence processing.

    Processes a sequence through the full pipeline: adapter removal, length
    normalization, adapter addition, and one-hot encoding. Handles both
    short sequences and long sequences (via sliding window). Always returns
    both original and reverse complement encodings.

    Args:
        seq (str): The input sequence from the API request
        prediction_ranges (list or None): List of prediction ranges to consider
        upstream_seq (str): Upstream flanking sequence provided by evaluator
        downstream_seq (str): Downstream flanking sequence provided by evaluator
        target_length (int): Length of core sequence before adapters (200 for MPRALegNet)
        seq_size (int): Final sequence size including adapters (230 for MPRALegNet)
        upstream_adapter (str): Upstream adapter sequence to prepend (HARDCODED 15bp for MPRALegNet)
        downstream_adapter (str): Downstream adapter sequence to prepend (HARDCODED 15bp for MPRALegNet)
        window_step (int): Step size for sliding window processing of long sequences

    Returns:
        list: List of tuples, each containing (encoded_seq, reversed_encoded_seq)
              For short sequences: single tuple in list
              For long sequences: multiple tuples for each window

    Logic:
        1. If constant flanks (upstream_seq / downstream_seq) provided:
            1.1. Check probe (variable region) length against target_length (200bp)
                    1.1.1. probe length <= target_length:
                        Pad with flanks to reach target_length (NOT Ns)
                        - Take tail of upstream_seq + probe + head of downstream_seq
                        - If target length not reached, pad with Ns
                        - Ignore prediction_ranges
                    1.1.2. probe length > target_length:
                        use sliding window, ignore prediction_ranges

        2. If no constant flanks:
            2.1. If probe length <= target length:
                    - Simply use the entire probe
                    - Pad with Ns, as necessary
            2.2. If probe length > target length and prediction_ranges provided (marks TSS/reporter location):
                    - MPRA-specific offset: probe ends 145bp upstream of TSS (MPRA_PROBE_TO_TSS_OFFSET)
                    - Region of interest:
                        start = max(0, PR_start - target_length - 145)
                        end = min(start + target_length, len(seq))
                        Pull as much real sequence context
                    - Pad with Ns if sequence still < target_length
            2.2. Else: use sliding window (probe length > target length and no pred_ranges)

        3. Add hardcoded MPRALegNet adapters (15bp upstream + 15bp downstream)
        4. One-hot encode the final sequence
    """

    # Bases may arrive in lower case (the request validation accepts them): upper-case the sequence and
    # the flanks once so that the adapters, the padding and the encoding all see upper case.
    seq = seq.upper()
    upstream_seq = upstream_seq.upper() if upstream_seq is not None else None
    downstream_seq = downstream_seq.upper() if downstream_seq is not None else None

    # First check if constant flanks are provided
    flanks_provided = bool(upstream_seq) or bool(downstream_seq)

    # Check and remove adapters if present (in case of pre-processed sequences)
    seq_without_adapters = remove_adapters(seq, upstream_adapter, downstream_adapter)

    if flanks_provided:
        # Branch 1: Flanks provided, ignore prediction_ranges
        #sys.stderr.write("LOG: Flanks provided. Applying flanks and ignoring prediction_ranges... \n")
        upstream_seq = upstream_seq if upstream_seq is not None else ""
        downstream_seq = downstream_seq if downstream_seq is not None else ""   

        seq_len = len(seq_without_adapters)

        if seq_len <= target_length:
            # Branch 1.1: Pad probe with flank context
            seq_with_flanks = _pad_with_flanks(
                seq_without_adapters, upstream_seq, downstream_seq, target_length
            )

            # Add hardcoded MPRALegNet adapters
            seq_with_adapters = add_adapters(seq_with_flanks, upstream_adapter, downstream_adapter)

            # Generate both original and reverse complement one-hot encodings
            encoded_seq = one_hot_encode(seq_with_adapters)
            reversed_seq = reverse_complement(seq_with_adapters)
            reversed_encoded_seq = one_hot_encode(reversed_seq)

            return [(encoded_seq, reversed_encoded_seq)]

        else:
            # Branch 1.2: probe is longer than target length - use sliding window
            windows = process_long_sequence(
                seq_without_adapters, target_length, window_step
            )

            results = []
            for window in windows:
                # Add adapters to each window
                window_with_adapters = add_adapters(
                    window, upstream_adapter, downstream_adapter
                )

                # Generate one-hot encodings for both strands for each window
                encoded_window = one_hot_encode(window_with_adapters)
                reversed_window = reverse_complement(window_with_adapters)
                reversed_encoded_window = one_hot_encode(reversed_window)

                results.append((encoded_window, reversed_encoded_window))

            return results
    else:
        # Branch 2: No flanks provided
        #sys.stderr.write("LOG: No flanks provided. Handling sequence based on length and prediction_ranges... \n")

        seq_len = len(seq_without_adapters)

        if seq_len <= target_length:
            # Branch 2.1: probe is shorter than target length, use as is and pad with Ns
            seq_without_adapters = pad_sequence(seq_without_adapters, target_length)

            # Add hardcoded MPRALegNet adapters
            seq_with_adapters = add_adapters(seq_without_adapters, upstream_adapter, downstream_adapter)

            # Generate both original and reverse complement one-hot encodings
            encoded_seq = one_hot_encode(seq_with_adapters)
            reversed_seq = reverse_complement(seq_with_adapters)
            reversed_encoded_seq = one_hot_encode(reversed_seq)

            return [(encoded_seq, reversed_encoded_seq)]
        
        elif prediction_ranges is not None:
            # Branch 2.2: probe is longer than target length, prediction_ranges provided
            # marks the TSS -- pull upstream window
            #sys.stderr.write("LOG: No flanks provided but prediction_ranges provided. Applying MPRA-specific TSS-anchored cropping... \n")
            _validate_prediction_range_bounds(prediction_ranges, len(seq))
            pr_start, _ = prediction_ranges
            seq_upstream_of_tss = _crop_upstream_of_tss(
                seq_without_adapters, pr_start, target_length, MPRA_PROBE_TO_TSS_OFFSET
            )
            if len(seq_upstream_of_tss) < target_length:
                seq_upstream_of_tss = pad_sequence(seq_upstream_of_tss, target_length)

            # Add hardcoded MPRALegNet adapters
            seq_with_adapters = add_adapters(seq_upstream_of_tss, upstream_adapter, downstream_adapter)

            # Generate both original and reverse complement one-hot encodings
            encoded_seq = one_hot_encode(seq_with_adapters)
            reversed_seq = reverse_complement(seq_with_adapters)
            reversed_encoded_seq = one_hot_encode(reversed_seq)

            return [(encoded_seq, reversed_encoded_seq)]
        else:
            # Branch 2.3: probe is longer than target length, no prediction_ranges provided -- crop to centre
            #sys.stderr.write("LOG: No flanks or prediction_ranges provided. Using sliding window... \n")
            
            windows = process_long_sequence(
                seq_without_adapters, target_length, window_step
            )

            results = []
            for window in windows:
                # Add adapters to each window
                window_with_adapters = add_adapters(
                    window, upstream_adapter, downstream_adapter
                )

                # Generate one-hot encodings for both strands for each window
                encoded_window = one_hot_encode(window_with_adapters)
                reversed_window = reverse_complement(window_with_adapters)
                reversed_encoded_window = one_hot_encode(reversed_window)

                results.append((encoded_window, reversed_encoded_window))

            return results

