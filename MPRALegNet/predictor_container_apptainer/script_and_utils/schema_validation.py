# schema_validation.py
"""Validation and Preprocessing of Payload"""

from error_checking_functions import *


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

    if "prediction_ranges" in payload:
        errors = check_seq_ids(
            payload["prediction_ranges"], payload["sequences"], errors
        )
        errors = check_prediction_ranges(
            payload["prediction_ranges"], payload["sequences"], errors
        )

    if "upstream_seq" in payload:
        errors = check_key_values_upstream_flank(payload["upstream_seq"], errors)
    if "downstream_seq" in payload:
        errors = check_key_values_downstream_flank(payload["downstream_seq"], errors)

    if any(errors.values()):
        flagged_errors = [msg for sublist in errors.values() for msg in sublist]
        raise BadRequestError(flagged_errors)


def preprocess_data(payload):
    """
    Handles model-specific data preprocessing and validation.

    This function performs model-specific checks and preprocessing on the payload,
    including sequence specifications, readout type compatibility, and species
    restrictions. It ensures the data meets the specific requirements of the
    MPRALegNet model.

    Args:
        payload (dict): The validated request payload to preprocess

    Raises:
        BadRequestError: If readout type is not supported
        PredictionFailedError: If sequence specifications or species are invalid

    Returns:
        dict: The processed sequences dictionary ready for model inference
    """
    sequences = payload.get("sequences", {})

    # Model-specific: MPRALegNet has hardcoded adapter sequences and ignores
    # flanking sequences from the request, but we inform the user about this
    if "upstream_seq" in payload or "downstream_seq" in payload:
        print(
            "INFO: MPRALegNet has adapter sequences hardcoded and does not use flanking sequences from the request."
        )

    # Check that the final sequences meet model specifications.
    # This includes character set.
    errors = {"prediction_request_failed": []}
    errors = check_seqs_specifications(sequences, errors)

    # Model-specific: Validate readout type compatibility
    # MPRaLegNet only supports "point" readout, not "track" or "interaction_matrix"
    readout_type = payload["readout"]
    if not readout_type in ["point"]:
        raise BadRequestError(
            f"MPRALegNet cannot process '{readout_type}' readout type."
        )

    # Model-specific: Ensure MPRALegNet Predictor only supports homo_sapiens
    # Check all prediction tasks for species compatibility
    for task in payload["prediction_tasks"]:
        if task.get("species").lower() != "homo_sapiens":
            errors["bad_prediction_request"].append(
                f"This predictor only supports species: homo_sapiens. Received '{task.get('species')}' for task '{task.get('name')}'."
            )

    # Raise errors if any model-specific validation checks failed
    if any(errors.values()):
        flagged_errors = [msg for sublist in errors.values() for msg in sublist]
        raise PredictionFailedError(flagged_errors)

    return sequences
