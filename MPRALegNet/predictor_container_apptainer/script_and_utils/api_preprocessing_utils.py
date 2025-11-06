# api_preprocessing_utils.py


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


# Full preprocessing pipeline for a sequence
def process_sequence(
    seq, target_length, seq_size, upstream_adapter, downstream_adapter, window_step=50
):
    """ "
    Complete preprocessing pipeline for MPRALegNet sequence processing.

    Processes a sequence through the full pipeline: adapter removal, length
    normalization, adapter addition, and one-hot encoding. Handles both
    short sequences and long sequences (via sliding window). Always returns
    both original and reverse complement encodings.

    Args:
        seq (str): The input sequence from the API request
        target_length (int): Length of core sequence before adapters (200 for MPRALegNet)
        seq_size (int): Final sequence size including adapters (230 for MPRALegNet)
        upstream_adapter (str): Upstream adapter sequence (15bp for MPRALegNet)
        downstream_adapter (str): Downstream adapter sequence (15bp for MPRALegNet)
        window_step (int): Step size for sliding window processing of long sequences

    Returns:
        list: List of tuples, each containing (encoded_seq, reversed_encoded_seq)
              For short sequences: single tuple in list
              For long sequences: multiple tuples for each window
    """
    # Step 1: Check and remove adapters if present (in case of pre-processed sequences)
    seq_without_adapters = remove_adapters(seq, upstream_adapter, downstream_adapter)

    # Step 2: Handle different sequence lengths
    seq_len = len(seq_without_adapters)

    if seq_len <= target_length:
        # Case 1: Sequence is <= target_length
        padded_seq = pad_sequence(seq_without_adapters, target_length)
        seq_with_adapters = add_adapters(
            padded_seq, upstream_adapter, downstream_adapter
        )

        # Generate both original and reverse complement one-hot encodings
        encoded_seq = one_hot_encode(seq_with_adapters)
        reversed_seq = reverse_complement(seq_with_adapters)
        reversed_encoded_seq = one_hot_encode(reversed_seq)

        return [(encoded_seq, reversed_encoded_seq)]

    else:
        # Case 2: Sequence is longer than target_length - use sliding window
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
