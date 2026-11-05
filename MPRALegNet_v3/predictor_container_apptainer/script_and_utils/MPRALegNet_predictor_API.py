# MPRALegNet_predictor_API.py

import os
import sys
import json
from flask import Flask

from config import legnet_DIR, HELP_FILE, PREDICTOR_NAME, SUPPORTED_REQUEST_FORMATS, SUPPORTED_RESPONSE_FORMATS
from error_checking_functions import *
from schema_validation import *
from predictor_content_handler import decode_request, encode_response

# Import the MPRALegNet prediction functions
from MPRALegNet_predict import predict_legnet

# Add DREAM_DIR to the Python path
if legnet_DIR not in sys.path:
    sys.path.append(legnet_DIR)

# --- Flask App and Central Error Handler ---
app = Flask(__name__)
# One of these works to maintain order when using jsonify()
app.config["JSON_SORT_KEYS"] = False
app.json.sort_keys = False


def create_error_response(error_key, messages, status_code):
    """
    Formats error response into a standarized JSON structure.

    Args:
        error_key (str): The category of the error (e.g. 'bad_prediction_request', 'prediction_request_failed').
        messages (list or str): A list of error message strings or a single message.
        status_code (int): Standard HTTP error status code based on the error.

    Returns:
        dict: A dictionary formatted for the standardized JSON error response.
    """
    if not isinstance(messages, list):
        messages = [str(messages)]
    error_payload = {"error": [{error_key: msg} for msg in messages]}
    print(error_payload)
    return error_payload, status_code


@app.errorhandler(APIError)
def handle_api_error(error):
    """This single handler catches all of our custom API errors."""
    # Get raw payload and status code
    payload, status_code = create_error_response(
        error.error_key, error.message, error.status_code
    )

    return encode_response(
        payload, status_code=status_code, isError=True, predictor_name=PREDICTOR_NAME
    )


@app.after_request
def after_request_callback(response):
    """This function runs after each request is processed."""
    print(f"\n--- Sending predictions back to Evaluator. ---")
    print(
        f"--- Request Complete. {PREDICTOR_NAME} Predictor is listening on http://{predictor_ip}:{predictor_port} ---\n"
    )
    return response


# --- API Endpoints ---
@app.route("/formats", methods=["GET"])
def formats_endpoint():
    """Provides the Predictor's supported formats"""
    supported_fmts = {
        "predictor_supported_request_formats": SUPPORTED_REQUEST_FORMATS,
        "predictor_supported_response_formats": SUPPORTED_RESPONSE_FORMATS,
    }
    try:
        return encode_response(
            supported_fmts,
            status_code=200,
            predictor_name=PREDICTOR_NAME,
            supported_response_formats=SUPPORTED_RESPONSE_FORMATS,
        )
    except Exception as e:
        raise ServerError(
            f"Error serializing supported format for /format endpoint: {e}"
        )


@app.route("/help", methods=["GET"])
def help_endpoint():
    """Provides the Predictor's help/metadata information."""
    try:
        with open(HELP_FILE, "r") as f:
            help_data = json.load(f)
        return encode_response(
            help_data,
            status_code=200,
            predictor_name=PREDICTOR_NAME,
            supported_response_formats=SUPPORTED_RESPONSE_FORMATS,
        )
    except Exception as e:
        raise ServerError(f"Error reading help file: {e}")


@app.route("/predict", methods=["POST"])
def predict():
    """
    The main endpoint for receiving sequences and returning predictions.
    
    This endpoint handles the complete prediction workflow:
    1. Decodes and validates incoming request
    2. Preprocesses sequence data for MPRALegNet model
    3. Runs model inference for each prediction task
    4. Formats and returns standardized prediction response
    
    The endpoint supports both JSON and MessagePack formats for request/response
    and handles various error conditions with appropriate HTTP status codes.
    
    Returns:
        Response: Flask response object containing predictions or error message
    """

    try:
        # Decode incoming request based on Content-Type header
        # Supports both JSON and MessagePack formats
        evaluator_request = decode_request(SUPPORTED_REQUEST_FORMATS)

        # Validate the payload using the imported function
        # These functions will raise an APIError on failure,
        # which will be caught automatically by @app.errorhandler
        validate_request_payload(evaluator_request)

        # Preprocess the data for MPRALegNet model compatibility
        # Handles sequence validation, adapter processing, and encoding
        sequences = preprocess_data(evaluator_request)
        readout_type = evaluator_request.get("readout")

        # ADDITION: Since we actually care about pred ranges, and upstream and downstream seqs
        prediction_ranges = evaluator_request.get("prediction_ranges")
        upstream_seq = evaluator_request.get("upstream_seq")
        downstream_seq = evaluator_request.get("downstream_seq")

        # Initialize response structure
        json_return = {"prediction_tasks": []}

        # Process each prediction task in the request
        for task in evaluator_request["prediction_tasks"]:
            # Determine cell type and scale for model selection
            # MPRALegNet supports K562, HepG2, and WTC11 cell types
            current_cell_type = "K562"  # Default fallback
            
            for possible_cell_type in ["K562", "HepG2", "WTC11"]:
                if possible_cell_type.lower() in task["cell_type"].lower():
                    current_cell_type = possible_cell_type
                    
            # Default output scale (MPRALegNet model was trained on log2(RNA/DNA))
            if task["scale"] in ["log", "linear"]:
                current_scale = task["scale"]
            else:
                current_scale = "log"

            # Run model inference based on readout type
            if readout_type == "point":
                # Call MPRALegNet prediction function with sequences and cell type
                # Uses sliding window with 50bp step for long sequences
                model_predictions = predict_legnet(sequences, prediction_ranges, upstream_seq, downstream_seq, current_cell_type, current_scale, 50)
            else:
                # MPRALegNet only supports "point" readout type
                raise PredictionFailedError(
                    f"Unsupported readout type: '{readout_type}'"
                )

            # Format response for this prediction task
            json_return["prediction_tasks"].append(
                {
                    "name": task["name"],                          # Task name from request
                    "type_requested": task["type"],                # Requested prediction type
                    "type_actual": "expression",                   # Actual prediction type delivered
                    "cell_type_requested": task["cell_type"],      # Requested cell type
                    "cell_type_actual": current_cell_type,         # Actual cell type used
                    "scale_prediction_requested": task["scale"],   # Requested output scale
                    "scale_prediction_actual": current_scale,      # Current output scale
                    "species_requested": task["species"],          # Requested species
                    "species_actual": "homo_sapiens",              # Actual species (MPRALegNet is human-specific)
                    "predictions": model_predictions,              # Model predictions dictionary
                }
            )

        # Encode and return successful response
        return encode_response(
            json_return,
            status_code=200,
            predictor_name=PREDICTOR_NAME,
            supported_response_formats=SUPPORTED_RESPONSE_FORMATS,
        )

    except Exception as e:
        # If it's already an APIError, re-raise it for the handler
        if isinstance(e, APIError):
            raise e
        # Otherwise, wrap the unknown error in a ServerError
        raise ServerError(f"An unexpected internal error occurred: {e}.")


if __name__ == "__main__":
    # Validate command line arguments
    if len(sys.argv) != 3:
        print(
            f"Invalid arguments! Arguments must have: <container image/python script> <ip_address> <port>"
        )
        sys.exit(1)

    # Extract IP and port from command line arguments
    predictor_ip = sys.argv[1]
    predictor_port = int(sys.argv[2])

    print(
        f"{PREDICTOR_NAME} Predictor is running on http://{predictor_ip}:{predictor_port}"
    )

    # Start Flask development server
    app.run(host=predictor_ip, port=predictor_port)
