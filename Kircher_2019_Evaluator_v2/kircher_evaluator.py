"""RESTful Test Evaluator Utilizing Flask"""

import os
import sys
import argparse
import json
import pandas as pd
from datetime import datetime, timezone
from requests.exceptions import RequestException, HTTPError

import config
from data_loader import load_and_validate_data
from evaluator_content_handler import *
import evaluator_metrics_calculator


def run_evaluator(predictor_ip, predictor_port, output_dir, length = 200):
    """
    Preprocesses the data, sends request, receives response,
    saves the response, and returns file path.
    """

    # Validate output directory exists; create if it does not
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(
            f"Output directory '{output_dir}' did not exist. Created it successfully!"
        )

    all_correlations = []

    # Iterate through all cell type and element combinations defined in config
    for cell_type, elements in config.ELEMENTS_BY_CELL_TYPE.items():
        for element in elements:
            # Construct path to sequence data file for current cell type and element
            file_path = os.path.join(
                config.EVALUATOR_DATA_DIR, cell_type, f"{element}_variants.tsv"
            )

            # Skip if data file doesn't exist for this combination
            if not os.path.exists(file_path):
                print(
                    f"Skipping {config.CELL_TYPES_NAMES[cell_type]}/{element} - no data"
                )
                continue

            print("=" * 50)
            print(f"----- Parsing {element} -----")
            print(f"----- in Cell Type {config.CELL_TYPES_NAMES[cell_type]} -----")
            print(f"Sequences length = {length}")
            print("Extracting...")

            # Load and validate genomic data, returns formatted JSON for predictor API
            data_dict = load_and_validate_data(
                cell_type, element, seq_length=length
            )

            # Count total sequences to verify predictor returns correct number of predictions
            total_sequences = len(data_dict["sequences"])

            # Communicate with Predictor, send request, receive predictions
            predictor_url = f"http://{predictor_ip}:{predictor_port}"
            response_payload = None
            is_success_response = (
                False  # Track whether we got a successful 200 response
            )

            try:
                # Negotiate data formats with predictor (JSON vs MessagePack)
                req_fmt, resp_fmt = negotiate_formats(predictor_url)

                # Send prediction request to predictor - raises HTTPError for 4xx/5xx responses
                response = get_predictions(predictor_url, data_dict, req_fmt, resp_fmt)

                is_success_response = True
                print("Predictor returned 200 OK.")

                # Deserialize successful response based on negotiated format
                response_payload = deserialize_response(response, resp_fmt)

            except HTTPError as http_err:
                # Handle HTTP errors (4xx client errors, 5xx server errors)
                print(
                    f"Predictor returned HTTP {http_err.response.status_code}. Processing error payload..."
                )
                is_success_response = False
                try:
                    # Deserialize the error body (could be JSON or msgpack)
                    response_payload = deserialize_response(http_err.response, resp_fmt)
                except ValueError as decode_err:
                    # Handle cases where the error response itself is malformed
                    print(
                        f"Could not decode the error response body: {decode_err}",
                        file=sys.stderr,
                    )

                    # Create fallback error structure
                    response_payload = {
                        "predictor_name": "UnknownPredictor_ErrorResponse",
                        "error": [
                            {
                                "server_error": f"Failed to decode error response (Status {http_err.response.status_code}). Body: {http_err.response.text[:500]}..."
                            }
                        ],
                    }

            if response_payload is None:
                # Handle complete failure to get any response
                print(
                    "FATAL: No response payload received or processed.", file=sys.stderr
                )
                response_payload = {
                    "predictor_name": "UnknownPredictor_NoResponse",
                    "error": [
                        {
                            "evaluator_error": "No response payload could be processed after request."
                        }
                    ],
                }

            # Save predictions to file with descriptive filename
            predictor_name = response_payload.get(
                "predictor_name", "UnknownPredictor"
            ).replace(" ", "_")
            output_filename_base = f"{cell_type}_cell_type_{element}_element_{config.EVALUATOR_NAME}_predictions"
            output_filename = f"{output_filename_base}_from_{predictor_name}.json"
            saved_predictions_path = os.path.join(output_dir, output_filename)

            # Validate that predictor returned correct number of predictions
            for i, task in enumerate(
                response_payload.get("prediction_tasks", []), start=1
            ):
                preds = task.get("predictions", [])
                if len(preds) != total_sequences:
                    print(
                        f"Warning: Task {i} ('{task.get('name')}') has {len(preds)} predictions, but {total_sequences} sequences were sent to the Predictor."
                    )

            # Save raw predictions to JSON file
            try:
                with open(saved_predictions_path, "w", encoding="utf-8") as f:
                    json.dump(response_payload, f, ensure_ascii=False, indent=4)
                print(f"Raw predictions saved to {saved_predictions_path}")
            except IOError as e:
                print(
                    f"FATAL: Could not save predictions to {saved_predictions_path}. {e}",
                    file=sys.stderr,
                )
                return

            # Calculate evaluation metrics only for successful responses
            if is_success_response:
                # Path to measured experimental values for comparison
                file_path = os.path.join(
                    config.EVALUATOR_DATA_DIR, cell_type, f"{element}_variants.tsv"
                )

                # Calculate correlation between predictions and experimental measurements
                corellation_results = (
                    evaluator_metrics_calculator.calculate_and_save_metrics(
                        predictions_content=response_payload,
                        measured_values_path=file_path,
                        cell_type=cell_type,
                        element=element,
                    )
                )

                if corellation_results:
                    all_correlations.extend(corellation_results)

            else:
                # Skip metrics calculation for failed predictor responses
                print(
                    "Skipping metrics calculation because the Predictor did not return a 200 OK status."
                )

    # Generate final summary report if we have correlation results
    if all_correlations:
        evaluation_output = pd.DataFrame(all_correlations)

        # Create timestamped output filename
        correlation_summary_filepath = f"correlation_summary_{config.EVALUATOR_NAME}_{datetime.now(timezone.utc).strftime('%Y%m%d-%H%M%S')}.tsv"
        full_output_path = os.path.join(output_dir, correlation_summary_filepath)

        print(f"Saving ALL correlations to: {full_output_path}")
        evaluation_output.to_csv(full_output_path, sep="\t", index=False)
        print(
            f"Successfully saved {len(all_correlations)} correlation results to {full_output_path}"
        )

        # Print summary statistics
        print("\n=== Correlation Summary ===")
        print(f"Total calculations: {len(all_correlations)}")
        print(
            f"Average Pearson r: {evaluation_output['Value'].astype(float).mean():.4f}"
        )
        print(f"Min Pearson r: {evaluation_output['Value'].astype(float).min():.4f}")
        print(f"Max Pearson r: {evaluation_output['Value'].astype(float).max():.4f}")
    else:
        print("No correlation results were calculated!")


def parse_arguments():
    # Extract command line arguments
    parser = argparse.ArgumentParser(description='Evaluator for regulatory genomics predictions')
    
    parser.add_argument('predictor_ip', help='IP address of the predictor service')
    parser.add_argument('predictor_port', type=int, help='Port of the predictor service')
    parser.add_argument('output_dir', help='Directory for output results')
    
    # Add length flag
    parser.add_argument(
        '--sequence-length', 
        '-l', 
        type=int, 
        default=200,  # Default
        help='Length of sequence to extract around each variant (default: 200)'
    )
    
    return parser.parse_args()

if __name__ == "__main__":

    # Parse arguments
    args = parse_arguments()
    
    # Length validation
    if args.sequence_length <= 0:
        print("ERROR: Sequence length must be positive", file=sys.stderr)
        sys.exit(1)
    
    if args.sequence_length > 1000000:  # Maximum 1Mb by default. Use create_evaluator_data.ipynb to vary 
        print("WARNING: Very large sequence length requested")
    
    try:
        # Run the main evaluation pipeline
        run_evaluator(
            predictor_ip=args.predictor_ip,
            predictor_port=args.predictor_port,
            output_dir=args.output_dir,
            length=args.sequence_length,
        )
        print("Evaluation complete.")
        sys.exit(0)
        
    except (FileNotFoundError, ValueError) as e:
        # Handle data-related errors
        print(f"FATAL ERROR (Data): {e}", file=sys.stderr)
        sys.exit(1)
    except RequestException as e:
        # Handle network connection errors
        print(
            f"FATAL ERROR (Network): Could not connect to predictor at http://{args.predictor_ip}:{args.predictor_port}.",
            file=sys.stderr,
        )
        sys.exit(1)
    except Exception as e:
        # Handle any unexpected errors
        print(f"An unexpected fatal error occurred: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)