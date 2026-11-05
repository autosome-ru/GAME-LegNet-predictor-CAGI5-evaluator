# schema_validation.py
"""Validation and Preprocessing of Payload"""
import numpy as np
from error_checking_functions import *

SUPPORTED_SCALES = ["linear", "log"]
DEFAULT_SCALE = "log"
# The species the model actually predicts for -- returned as `species_actual`.
# `species` in a request is the species of the cells the experiment was performed in.
MODEL_SPECIES = "homo_sapiens"

def unsupported_task_reason(task):
    """
    Checks whether MPRALegNet can perform a single prediction task, based on its species and type.

    A task MPRALegNet cannot perform does NOT fail the whole request: it is returned with
    null `*_actual` keys and an error object in `predictions` (see Predictor API Schema),
    while the other tasks are still predicted. (Any `cell_type` is accepted: K562, HepG2 and
    WTC11 each have their own model; any other cell type is predicted with the K562 model and
    reported as `cell_type_actual: "K562"`.)

    Args:
        task (dict): One object of `prediction_tasks` from the Evaluator request.

    Returns:
        str or None: A message describing why the task cannot be performed, or None if it can.
    """
    reasons = []

    # --- MODEL SPECIFIC: MPRALegNet only supports homo_sapiens ---
    if task.get('species', '').lower() != MODEL_SPECIES:
        reasons.append(
            f"This predictor only supports species: {MODEL_SPECIES}. Received '{task.get('species')}'."
        )

    # MPRALegNet predicts "expression" and the "expression_*" subtypes (e.g. expression_mRNA,
    # expression_pol1, expression_pol2, expression_pol3 -- see API spec).
    task_type = task.get('type', '')
    if not (task_type == "expression" or task_type.startswith("expression_")):
        reasons.append(
            f"prediction type '{task_type}' is not supported by MPRALegNet. "
            "Only 'expression' and 'expression_*' subtypes are accepted."
        )

    return " ".join(reasons) if reasons else None

def validate_request_payload(payload):
    """
    Performs all validation checks on the incoming request payload.
    Returns a dictionary of errors. If the dictionary is empty, validation passed.

    Args:
        payload (dict): The incoming request payload to validate

    Raises:
        BadRequestError: If any validation checks fail, with detailed error messages
        about what aspects of the payload are invalid

    Returns:
        None: Validation passes silently if no errors are found
    """
    errors = {"bad_prediction_request": []}

    # First confirm all mandatory keys are present
    errors = check_mandatory_keys(payload.keys(), errors)
    if any(errors.values()):
        flagged_errors = [msg for sublist in errors.values() for msg in sublist]
        raise BadRequestError(flagged_errors)

    # Check for mandatory keys inside each task object.
    errors = check_prediction_task_mandatory_keys(payload["prediction_tasks"], errors)
    if any(errors.values()):
        # Fail immediately if any task is missing keys, before we try to access them.
        flagged_errors = [msg for sublist in errors.values() for msg in sublist]
        raise BadRequestError(flagged_errors)

    # Perform all other validation checks
    errors = check_key_values_readout(payload["readout"], errors)
    errors = check_prediction_task_name(payload["prediction_tasks"], errors)
    errors = check_prediction_task_type(payload["prediction_tasks"], errors)
    errors = check_prediction_task_cell_type(payload["prediction_tasks"], errors)
    errors = check_prediction_task_species(payload["prediction_tasks"], errors)
    errors = check_prediction_task_scale(payload["prediction_tasks"], errors)

    if 'prediction_ranges' in payload:
        errors = check_seq_ids(payload['prediction_ranges'], payload['sequences'], errors)
        errors = check_prediction_ranges(payload['prediction_ranges'], errors)

    if 'upstream_seq' in payload:
        errors = check_key_values_upstream_flank(payload['upstream_seq'], errors)
    if 'downstream_seq' in payload:
        errors = check_key_values_downstream_flank(payload['downstream_seq'], errors)

    if any(errors.values()):
        flagged_errors = [msg for sublist in errors.values() for msg in sublist]
        raise BadRequestError(flagged_errors)


def preprocess_data(payload):
    """
    Handles model-specific data preprocessing and validation.

    This function performs model-specific checks and preprocessing on the payload,
    including sequence specifications and readout type compatibility.
    It ensures the data meets the specific requirements of the
    MPRALegNet model.

    Args:
        payload (dict): The validated request payload to preprocess

    Raises:
        PredictionFailedError: If the readout type is not supported or sequence specifications are invalid

    Returns:
        dict: The processed sequences dictionary ready for model inference
    """
    sequences = payload.get("sequences", {})

    if "upstream_seq" in payload or "downstream_seq" in payload:
        print(
            "INFO: Upstream/Downstream sequences found. They will be used to construct the biological context before cropping/padding."
        )

    # Check that the final sequences meet model specifications.
    # This includes character set.
    errors = {
        'bad_prediction_request': [],
        'prediction_request_failed': []
        }
    errors = check_seqs_specifications(sequences, errors)

    # Model-specific: Validate readout type compatibility
    # MPRaLegNet only supports "point" readout, not "track" or "interaction_matrix"
    readout_type = payload["readout"]
    if not readout_type in ["point"]:
        raise PredictionFailedError(
            f"MPRALegNet cannot process '{readout_type}' readout type."
        )

    # Task-level limitations (species, type) are handled per task by `unsupported_task_reason`.
    for task in payload["prediction_tasks"]:
        req_scale = task.get('scale')
        if req_scale and req_scale.lower() not in SUPPORTED_SCALES:
            errors['bad_prediction_request'].append(
                f"Unsupported scale: '{req_scale}'. Supported scales are: {SUPPORTED_SCALES}."
            )

    # Raise errors if any model-specific validation checks failed
    if any(errors.values()):
        flagged_errors = [msg for sublist in errors.values() for msg in sublist]
        raise PredictionFailedError(flagged_errors)

    return sequences
