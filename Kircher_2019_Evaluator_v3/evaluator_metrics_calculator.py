"""Calculate the final evaluation metrics."""

# NOTE: Every evaluator will do this slightly differently depending on how the data is presented

from scipy.stats import pearsonr
import pandas as pd
import numpy as np
import sys
from datetime import datetime, timezone
from config import CELL_TYPES_NAMES, EVALUATOR_NAME


def calculate_and_save_metrics(
    predictions_content, measured_values_path, cell_type=None, element=None, sequence_ids=None
):
    """
    Calculates custom evaluation metrics by comparing model predictions with experimental measurements.

    This function serves as the main entry point for evaluation. It processes prediction tasks,
    validates the data structure, calculates Pearson correlations between predicted and measured
    variant effects, and returns formatted results.

    Args:
        predictions_content (dict): Dictionary containing model predictions with
                                  'prediction_tasks' structure. Each task should contain
                                  predictions for reference and alternative sequences.
        measured_values_path (str): Path to TSV file containing experimentally
                                   measured variant effects for validation.
        cell_type (str or list): Cell type(s) to filter variants by.
        element (str or list): Genomic element type(s) to filter by
                                        (e.g., 'promoter', 'enhancer').
        sequence_ids (keys, optional): IDs of the sequences sent to the Predictor, in the order sent. If given,
                                  a task is only scored if its predictions are for exactly these sequences.

    Returns:
        list: A list of dictionaries containing correlation results for each prediction task.
              A metric that cannot be calculated is recorded as "NaN".
    """

    print("----- Starting Evaluation Calculation -----")

    try:
        # Extract predictor name and prediction tasks from input
        predictor_name = predictions_content.get("predictor_name", "Unknown")
        all_tasks = predictions_content.get("prediction_tasks", [])

        # An error response (e.g. an HTTP error) has no prediction tasks: record "NaN" for this cell type and element
        if not all_tasks:
            print("WARNING: 'prediction_tasks' key missing or empty. The metric will be recorded as NaN.")
            all_tasks = [{"cell_type_requested": CELL_TYPES_NAMES[cell_type]}]

        # Pre-computation: Validate all tasks *once*
        # This map stores whether each task is valid and its metadata
        task_results = {}
        for task in all_tasks:
            if not isinstance(task, dict):
                continue

            task_name = task.get("cell_type_requested", "unknown")
            predictions = task.get("predictions", {})
            task_valid = True  # Assume valid until proven otherwise

            # Check if task contains errors or has no predictions
            if not isinstance(predictions, dict) or not predictions or "error" in predictions:
                print(
                    f"- Task (Cell: {task_name}): Found 'error' or no predictions. Marking as invalid."
                )
                task_valid = False
            # Predictions must be for exactly the sequences that were sent (not fewer, and no extra or different IDs)
            elif sequence_ids is not None and predictions.keys() != set(sequence_ids):
                print(
                    f"- Task (Cell: {task_name}): Predictions are not for exactly the {len(sequence_ids)} sequences that were sent. Marking as invalid."
                )
                task_valid = False

            # Store task metadata and predictions for later processing
            task_results[task_name] = {
                "is_valid": task_valid,
                "metadata": {k: v for k, v in task.items() if k != "predictions"},
                "predictions": predictions,
                "cell_type_actual": task.get("cell_type_actual", "unknown"),
                "cell_type_requested": task.get("cell_type_requested", "unknown"),
                "scale_prediction_actual": task.get(
                    "scale_prediction_actual", "linear"
                ),
            }

        # Calculate correlation metrics for all valid tasks
        correlation_results = _calculate_correlations(
            predictions_content=task_results,
            measured_values_path=measured_values_path,
            predictor_name=predictor_name,
            cell_type=cell_type,
            element=element,
            sequence_ids=sequence_ids,
        )

        return correlation_results

    except Exception as e:
        print(
            f"An unexpected error occurred during evaluation calculations: {e}",
            file=sys.stderr,
        )
        import traceback

        traceback.print_exc()


def _calculate_correlations(
    predictions_content, measured_values_path, predictor_name, cell_type, element, sequence_ids=None
):
    """
    Calculate Pearson correlations for each prediction task.

    Iterates through all prediction tasks, computes Pearson correlation between
    predicted and measured variant effects for valid tasks, and formats the results.

    Args:
        predictions_content (dict): Dictionary containing validated prediction tasks
        measured_values_path (str): Path to experimental measurements TSV file
        predictor_name (str): Name of the predictor/model being evaluated
        cell_type (str or list, optional): Cell type filter
        element (str or list, optional): Genomic element filter
        sequence_ids (keys, optional): IDs of the sequences sent to the Predictor, in the order sent

    Returns:
        list: Formatted correlation results for all tasks
    """

    all_task_correlation_results = []

    # Process each prediction task
    for actual_cell_type, results in predictions_content.items():
        # Get the scale of predictions (log or linear) for proper processing
        scale_prediction_actual = results.get("scale_prediction_actual", "linear")

        # Calculate Pearson correlation for valid tasks, return NaN for invalid ones
        pearson_r = (
            calculate_pearson_r(
                predictions_content=results,
                measured_values_path=measured_values_path,
                scale_prediction_actual=scale_prediction_actual,
                cell_type=cell_type,
                element=element,
                sequence_ids=sequence_ids,
            )
            if results["is_valid"]
            else float("nan")
        )

        # Print correlation results for monitoring
        print(f"----- Calculated Pearson r: ----- ")
        print(f"Predictor: {predictor_name}")
        print(f"Cell type requested: {results['cell_type_requested']}")
        print(f"Cell type ACTUAL: {results['cell_type_actual']}")
        print(f"Pearson r = {pearson_r} for {CELL_TYPES_NAMES[cell_type]}/{element}")
        print("=" * 50)

        # Create description and timestamp for results
        description = f"Kircher Variant Effects ({element}), Cell_Type_Requested - {CELL_TYPES_NAMES[cell_type]}, Cell_Type_Actual - {results["cell_type_actual"]}"
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")

        # Format results for output
        all_task_correlation_results.append(
            {
                "evaluator_name": EVALUATOR_NAME,
                "description": description,
                "predictor_name": predictor_name,
                "time_stamp": timestamp,
                "metric": "pearson_r",
                # A metric that could not be calculated is recorded as "NaN"
                "value": "NaN" if pearson_r is None or not np.isfinite(pearson_r) else str(pearson_r),
                "prediction_task(s)_data": str(results["metadata"]),
            }
        )

    return all_task_correlation_results


def calculate_pearson_r(
    predictions_content,
    measured_values_path,
    scale_prediction_actual,
    cell_type,
    element,
    sequence_ids=None,
):
    """
    Calculate Pearson correlation between predicted and measured variant effects.

    This function processes model predictions to compute variant effects between
    reference and alternative sequences, then correlates them with experimentally
    measured values. Handles both log-scale and linear-scale predictions.

    Args:
        predictions_content (dict): Dictionary containing model predictions with
                                  'predictions' key containing sequence predictions
        measured_values_path (str): Path to TSV file containing experimentally
                                   measured variant effects (log2 ref - log2 alt)
        scale_prediction_actual (str): Scale of predictions - 'log' for log-scale
                                      predictions, 'linear' or None for linear-scale
        cell_type (str or list, optional): Cell type(s) to filter variants by
        element (str or list, optional): Genomic element type(s) to filter by
        sequence_ids (keys, optional): IDs of the sequences sent to the Predictor, in the order sent

    Returns:
        float: Pearson correlation coefficient between predicted and measured variant effects,
               or None if calculation fails due to insufficient data or mismatches
    """

    # Validate input structure contains predictions
    if not predictions_content["predictions"]:
        raise ValueError("No prediction tasks found in predictions_content")

    try:
        # Extract predictions dictionary
        predictions_dict = predictions_content["predictions"]
    except (KeyError, IndexError) as e:
        raise ValueError(f"Invalid predictions structure: {e}")

    # Create DataFrame from Predictions
    # Convert dictionary to DataFrame with sequence IDs and predicted values
    predictions_df = pd.DataFrame(
        list(predictions_dict.items()), columns=["id_column", "Predicted_Value"]
    )

    # Put the predictions in the order the sequences were sent (reference, then alternate, for each variant),
    # whatever order the Predictor returned them in
    if sequence_ids is not None:
        predictions_df = predictions_df.set_index("id_column").loc[list(sequence_ids)].reset_index()

    # Extract single value if predictions are in list format
    predictions_df["Predicted_Value"] = predictions_df["Predicted_Value"].apply(
        lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x
    )

    # A prediction that is not a finite number invalidates the task, so that a Predictor cannot choose
    # which variants are scored by leaving some out
    predictions_df["Predicted_Value"] = pd.to_numeric(
        predictions_df["Predicted_Value"], errors="coerce"
    ).replace([np.inf, -np.inf], np.nan)
    if predictions_df["Predicted_Value"].isna().any():
        print("WARNING: Some predictions are not finite numbers. Skipping evaluation calculation for this task.")
        return None

    # Calculate variant effects based on prediction scale
    if scale_prediction_actual == "log":
        # For log-scale predictions: variant effect = log2(alternate) - log2(reference)
        # Assumes alternating reference/alternative pairs in the predictions
        value_predicted = (
            predictions_df["Predicted_Value"].iloc[1::2].values
            - predictions_df["Predicted_Value"].iloc[0::2].values
        )

    else:  # linear scale
        # For linear-scale predictions: variant effect = log2(alternate/reference)
        # Extract reference and alternative values from alternating positions
        ref_values = predictions_df["Predicted_Value"].iloc[0::2].values
        alt_values = predictions_df["Predicted_Value"].iloc[1::2].values

        # Avoid division by zero
        with np.errstate(divide="ignore", invalid="ignore"):
            value_predicted = np.log2(
                np.divide(alt_values, ref_values, where=ref_values != 0)
            )

        # Replace infinite values with NaN for clean processing
        value_predicted = np.where(
            np.isfinite(value_predicted), value_predicted, np.nan
        )

    # Create a results DataFrame with reference/alternative pairs and calculated variant effects
    delta_df = pd.DataFrame(
        {
            "ref_id": predictions_df["id_column"].iloc[0::2].values,
            "alt_id": predictions_df["id_column"].iloc[1::2].values,
            "Value_predicted": value_predicted,
        }
    )

    # Extract variant name from ID by removing reference/alternative suffixes
    delta_df["Name"] = delta_df["ref_id"].str.replace(
        r"_(?:reference|alternate)$", "", regex=True
    )
    # print(delta_df.shape)

    # Load experimentally measured data for validation
    try:
        print(f"Using measured data from: {measured_values_path}")
        measured_values = pd.read_csv(measured_values_path, sep="\t")
        number_of_sequences = len(measured_values)
        print(f"Found {number_of_sequences} sequences in measured data for comparison.")
    except Exception as e:
        print(
            f"FATAL: Could not load measured data from {measured_values_path} to get sequence count. {e}",
            file=sys.stderr,
        )
        return

    # Filter measured data by cell types if specified
    if cell_type is not None:
        # Convert single cell type to list for consistent filtering
        if not isinstance(cell_type, list):
            cell_type = [cell_type]

        measured_values = measured_values[measured_values.Cell_Type.isin(cell_type)]

    # Filter measured data by genomic elements if specified
    if element is not None:
        # Convert single element to list for consistent filtering
        if not isinstance(element, list):
            element = [element]

        measured_values = measured_values[measured_values.Element.isin(element)]

    # Merge predicted and measured data on variant names
    merged_data = measured_values.merge(delta_df, how="left", on="Name")
    merged_data = merged_data.dropna()

    print(f"Result {len(merged_data)} predicted sequences after processing")

    # Check if we have sufficient data for correlation calculation
    if len(merged_data) == 0:
        print("WARNING: No overlapping variants found after merging and cleaning")
        return None

    # Valid predictions (or measurements) that do not vary have no correlation: 0.0, not NaN
    if merged_data["Value_predicted"].std() == 0 or merged_data["Value"].std() == 0:
        print("WARNING: The predicted or the measured variant effects do not vary. Pearson r = 0.0")
        return 0.0

    # Calculate Pearson correlation between predicted and measured variant effects
    try:
        r, p_value = pearsonr(merged_data["Value"], merged_data["Value_predicted"])
        return 0.0 if np.isnan(r) else float(r)
    except Exception as e:
        print(f"Error calculating Pearson correlation: {e}")
        return None
