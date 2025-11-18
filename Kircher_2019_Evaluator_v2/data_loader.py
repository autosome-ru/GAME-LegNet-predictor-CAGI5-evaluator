"""Handle Loading and Validating Evaluator Input/Request Data"""

import os
import json
import pandas as pd
import pyfastx
from collections import Counter
import functools

from config import CELL_TYPES_NAMES, EVALUATOR_DATA_DIR


class DuplicateKeysError(ValueError):
    """Raised when duplicate keys are found in a JSON object."""

    pass


# Internal helper function to detect duplicates during JSON parsing
def _detect_duplicates(pairs, duplicate_keys_state):
    """
    Detects duplicate keys during JSON parsing and counts occurrences of each key.

    This function intercepts the key-value pairs provided by `json.loads` and ensures that
    duplicate keys are flagged. It constructs the dictionary normally but counts how often
    each key appears, recording any keys that occur more than once.

    Args:
        pairs (list of tuple): A list of key-value pairs at the current level of the JSON.
        duplicate_keys_state (dict): The dictionary to update with any duplicates found.

    Returns:
        result_dict: A dictionary created from the key-value pairs.
    """

    # Use a local Counter to count occurrences of keys at this level
    local_counts = Counter()
    result_dict = {}
    for key, value in pairs:
        # Increment the count for each key
        local_counts[key] += 1
        # If the key is a duplicate, record it in the duplicate_keys dictionary
        if local_counts[key] > 1:
            duplicate_keys_state[key] = local_counts[key]
        # Add the key-value pair to the resulting dictionary
        result_dict[key] = value
    return result_dict


def _process_results(data, duplicate_keys):
    """
    Checks the duplicate_keys dictionary and prints a report.

    Args:
        data (dict): The dictionary of parsed data.
        duplicate_keys (dict): The dictionary of duplicates.

    Returns:
        data or None: The parsed data if no duplicates. None, if duplicates are found.
    """
    # Report duplicates if any were found
    if duplicate_keys:
        print("Duplicate keys found:")
        error_messages = [
            f"Key: '{key}', Count: {count}" for key, count in duplicate_keys.items()
        ]
        raise DuplicateKeysError(f"Duplicate keys found:\n" + "\n".join(error_messages))
    else:
        print("No duplicates found.")
        return data  # Return the parsed data if no duplicates.


# Function to check for duplicate keys in JSON object


def check_duplicates_from_string(json_string):
    """
    Parses a JSON string to detect and report any duplicate keys at the same level in the same object.
    This function ensures that no keys are silently overwritten in dictionaries.

    The function uses a helper to track the number of times each key appears during parsing,
    leveraging the `object_pairs_hook` parameter of `json.loads()` to intercept key-value pairs
    before they are processed into a dictionary. If duplicates are detected at any level, they
    are reported with their counts. Keys reused in separate objects within arrays (e.g. lists)
    are not considered duplicates.

    Args:
        json_string (str): The JSON content as a string to parse and check for duplicates.

    Raises:
        json.JSONDecodeError: If the string is not valid JSON.
        DuplicateKeysError: If duplicate keys are found in the JSON structure.

    Returns:
        dict: The parsed data if no errors or duplicates are found.
    """

    # Initialize a dictionary to track duplicate keys and their counts
    duplicate_keys = {}

    # Create a 1-argument hook callable by "freezing" the duplicate_keys dict
    # as the second argument to the helper.
    hook = functools.partial(_detect_duplicates, duplicate_keys_state=duplicate_keys)

    # Parse the JSON string using the helper to track duplicates
    data = json.loads(json_string, object_pairs_hook=hook)

    return _process_results(data, duplicate_keys)


# Function for check for duplicate keys if input file is in JSON format


def check_duplicates_from_json(json_file_path):
    """
    Parses a JSON file to detect and report any duplicate keys at the same level in the same object.
    This function ensures that no keys are silently overwritten in dictionaries.

    The function uses a helper to track the number of times each key appears during parsing,
    leveraging the `object_pairs_hook` parameter of `json.load()` to intercept key-value pairs
    before they are processed into a dictionary. If duplicates are detected at any level, they
    are reported with their counts and paths. Keys reused in separate objects within arrays
    (e.g. lists) are not considered duplicates.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        json.JSONDecodeError: If the file content is not valid JSON.
        DuplicateKeysError: If duplicate keys are found in the JSON structure.

    Returns:
        dict: The parsed data if no errors or duplicates are found.
    """

    # Initialize a dictionary to track duplicate keys and their counts
    duplicate_keys = {}

    # Create a 1-argument hook callable by "freezing" the duplicate_keys dict
    # as the second argument to the helper.
    hook = functools.partial(_detect_duplicates, duplicate_keys_state=duplicate_keys)

    # Open and parse the JSON file, using the helper to track duplicates
    with open(json_file_path, "r") as file:
        data = json.load(file, object_pairs_hook=hook)

    return _process_results(data, duplicate_keys)


def load_and_validate_data(cell_type, element, seq_length=200):
    """
    Loads genomic variant data and creates formatted JSON.

    This function serves as the main data preparation pipeline. It:
    1. Loads and filters variant data from TSV files
    2. Extracts reference and alternative sequences around variants
    3. Formats data into standardized JSON structure for prediction requests

    The output JSON follows the ML predictor API specification with prediction tasks
    and sequences for both reference and alternative alleles.

    Args:
        cell_type (str): Cell type to filter variants by (e.g., 'K562', 'HepG2')
        element (str): Genomic element type to filter by (e.g., 'F9', 'PLKR')
        seq_length (int, optional): Length of sequences to extract around variants.
                                   Defaults to 200bp.

    Raises:
        ValueError: If the file type is unsupported or data is malformed

    Returns:
        dict: Validated data dictionary formatted for predictor API, containing:
              - request type ("predict")
              - readout type ("point")
              - prediction tasks metadata
              - sequences for reference and alternative alleles
    """
    
    # Load variant data from TSV file
    file_path = os.path.join(EVALUATOR_DATA_DIR, cell_type, f"{element}_variants.tsv")

    sequence_dataFrame = pd.read_csv(file_path, sep="\t")
    sequence_dataFrame = sequence_dataFrame[[
        "Name",
        "Chromosome",
        "Position",
        "Ref",
        "Alt",
        "Element",
        "Cell_Type",
    ]]

    # Filter variants by specified cell type and genomic element
    filtered_df = sequence_dataFrame[
        (sequence_dataFrame["Cell_Type"] == cell_type)
        & (sequence_dataFrame["Element"] == element)
    ]

    if filtered_df.empty:
        print(
            f"Warning: No data found for cell_type '{cell_type}' and element '{element}'"
        )
        return None

    # Load extracted reference genome
    genome_file_path = os.path.join(EVALUATOR_DATA_DIR, cell_type, f"{element}_sequence.fasta.gz")

    try:
        genome = pyfastx.Fastx(genome_file_path, comment = True)
    except Exception as e:
        raise IOError(f"Error loading FASTA file {genome_file_path}: {str(e)}") from e

    for name, seq, comment in genome:
        # name = name + "_" + comment
        genome_seq = seq

    # Load Start and End porsitions of extracted genome sequence
    metadata_file_path = os.path.join(EVALUATOR_DATA_DIR, cell_type, f"{element}_metadata.json")

    with open(metadata_file_path, "r") as file:
        metadata = json.load(file)

    genome_start = metadata["start"]
    genome_end = metadata["end"]

    # Create dictionary mapping sequence names to their genomic sequences
    variant_dict = {}

    for _, row in filtered_df.iterrows():
        variant = row["Name"]

        # Extract reference sequence (wild-type allele)
        ref_sequence = get_sequence(
            genome = genome_seq,
            genome_start = genome_start,
            genome_end = genome_end,
            length=seq_length,
            pos=row.Position,
            ref=row.Ref,
            alt=row.Ref,  # Use ref again to get unmodified sequence
        )

        if ref_sequence is not None:
            variant_dict[f"{variant}_reference"] = ref_sequence

        # Extract alternative sequence (mutant allele with SNP variant)
        alt_sequence = get_sequence(
            genome = genome_seq,
            genome_start = genome_start,
            genome_end = genome_end,
            length=seq_length,
            pos=row.Position,
            ref=row.Ref,
            alt=row.Alt,
        )

        if alt_sequence is not None:
            variant_dict[f"{variant}_alternate"] = alt_sequence

    # Define prediction task metadata for the API
    prediction_task = [
        {
            "name": f"Kircher_{element}",
            "type": "expression",
            "cell_type": CELL_TYPES_NAMES[cell_type],
            "scale": "log",  # Model should return log-scale predictions
            "species": "homo_sapiens",
        }
    ]

    # Build the complete evaluator object for predictor API
    evaluator_dict = {
        "request": "predict",
        "readout": "point",  # Single value prediction per sequence
        "prediction_tasks": prediction_task,
        "sequences": variant_dict,
    }

    json_string = json.dumps(evaluator_dict, indent=4)

    # Validate JSON structure and check for duplicate keys
    data_dict = check_duplicates_from_string(json_string)

    return data_dict


def get_sequence(genome, genome_start, genome_end, length, pos, ref, alt):
    """
    Extract genomic sequence centered around a variant position with fixed length.

    This function extracts sequences from the reference genome and incorporates
    genetic variants (SNPs or deletions). Sequences are centered around the variant
    position and padded to the specified length. Handles both substitutions and
    deletions by modifying the reference sequence accordingly.

    Args:
        genome (str): Loaded extracted genome Fasta long (more 1 mb) sequence
        genome_start (str): Coordinate where long extracted genome sequence starts
        genome_end (str): Coordinate where long extrancted genome sequence ends
        length (int): Total length of sequence to extract (centered around variant)
        pos (int): Genomic position (1-based coordinate system)
        ref (str): Reference nucleotide at the position (single character)
        alt (str): Alternative nucleotide ('-' for deletion, single char for substitution)

    Returns:
        str or None: Extracted and modified genomic sequence in uppercase,
                    or None if reference nucleotide doesn't match or extraction fails

    Raises:
        ValueError: If reference nucleotide is not a single character or position is invalid
    """

    # Validate reference nucleotide format
    if not isinstance(ref, str) or len(ref) != 1:
        raise ValueError(f"Reference nucleotide should be single character, got {ref}")

    # Ensure start, end positions and variant coordinate are integers
    try:
        pos = int(pos)
        genome_start = int(genome_start) 
        genome_end = int(genome_end)
    except (ValueError, TypeError) as e:
        raise ValueError(f"genome start, genome end, variant position must be integers, got {genome_start}, {genome_end}, {pos}: {e}")

    # Verify reference nucleotide matches the genome at this position
    observed_ref = genome[pos - genome_start].upper()
    if observed_ref != ref.upper():
        return None

    # Calculate sequence coordinates centered around variant position
    half_len = length // 2
    start = pos - half_len - genome_start
    end = pos + half_len - genome_end

    # Adjust end coordinate for even-length sequences to maintain exact length
    if length % 2 == 0:
        end -= 1

    try:
        # Position of the variant within the extracted sequence
        ref_pos_in_seq = half_len

        if alt == "-":
            # Handle deletion: remove reference nucleotide from sequence
            # Extract sequence including one extra base to accommodate deletion
            seq = genome[start: end + 1]

            # Remove the reference nucleotide at variant position
            modified_seq = seq[:ref_pos_in_seq] + seq[ref_pos_in_seq + 1 :]
        else:
            # Handle substitution: replace reference with alternative nucleotide
            seq = genome[start: end]

            # Replace reference nucleotide with alternative
            modified_seq = seq[:ref_pos_in_seq] + alt + seq[ref_pos_in_seq + 1 :]

        return modified_seq.upper()

    except Exception as e:
        print(f"Error processing {pos}-{ref}>{alt}: {str(e)}")
        return None