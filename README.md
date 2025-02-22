# PtWAVE implementation summary

This repository provides Python-based code examples of the main algorithms implemented in PtWAVE (https://www.ptwave-ptbio.com) version 1.

*Notes: This code has been modified to explain our algorithms clearly. The actual PtWAVE server includes more complex data processing less relevant to the algorithms.

## Modeling of DNA Sequence Trace Chromatograms

1. EMSPs_Optimizer Class

This class is responsible for generating optimized sets of EMSPs. It uses a generator function select_EMSPs() to iteratively select and refine EMSP sets. Two modes are supported:

- Random Mode (random):
Randomly removes EMSPs to create different subsets.
Generates a coefficient matrix based on the selected EMSPs.

- Backstep Mode (backstep):
Removes a small percentage of EMSPs and evaluates whether the new set improves performance.
If the change is beneficial, the new EMSP set is accepted; otherwise, the previous set is retained.

2. Modeling_Chromatograms() Function

This function performs signal optimization by selecting an optimal set of EMSPs using different strategies:

- All Mode (all):
Uses all available EMSPs and generates a coefficient matrix.
Evaluates the model using Bayesian Information Criterion (BIC) and R² values.

- Random Mode (random):
Tries 10 different randomly selected EMSP sets.
Picks the one with the lowest BIC score.

- Backstep Mode (backstep):
Iteratively removes EMSPs, evaluates changes, and decides whether to accept or reject the new set.
Ensures that R² remains above 0.8 while minimizing BIC.

3. Model Evaluation

The function fit_evaluate_model() calculates performance metrics based on fitting result of non-negative linear modeling (NNLS) or non-negative LASSO regression:

Pearson's correlation coefficient (R²): Measures how well the predicted values fit the observed data.
Residual Sum of Squares (RSS): Determines the error in predictions.
Bayesian Information Criterion (BIC): Penalizes model complexity.

4. Final Adjustments

The optimized EMSP set is used to recalculate relative proportions (EMSP_rel) for each EMSP.
These proportions are adjusted based on the final R² value.

```python

def Modeling_Chromatograms(results, output_vec, decomposition_window, EMSPs, indel_max_size, fitting_algorithm, mode=None):
    """
    Optimize EMSPs based on the specified mode and find the optimal BIC value.
    """

    available_mutation_types = list(range(-indel_max_size * 2, indel_max_size + 1))
    eo = EMSPs_Optimizer(EMSPs, decomposition_window)

    if mode == "all":
        # Use all EMSPs for optimization
        coefficient_matrix = generate_coefficient_matrix(EMSPs, decomposition_window)
        bic_tuple = fit_evaluate_model(coefficient_matrix, output_vec, fitting_algorithm)

    elif mode == "random":
        # Try different random EMSP set and find the one with the lowest BIC
        EMSPs_selection = eo.select_EMSPs(available_mutation_types, "random")
        bic_tuple = min((fit_evaluate_model(EMSPs_selection.__next__()[1], output_vec, fitting_algorithm) for _ in range(10)), key=lambda x: x[0])

    elif mode == "backstep":
        # Reduce EMSPs step by step while optimizing
        min_bic_tuple = None
        for n in range(10):
            EMSPs_selection = eo.select_EMSPs(available_mutation_types, "backstep", random_seed=n)
            step_count = 0
            while step_count < 10:
                EMSPs, coefficient_matrix = EMSPs_selection.__next__()
                bic_tuple = fit_evaluate_model(coefficient_matrix, output_vec, fitting_algorithm)
                if step_count == 0 or (bic_tuple[0] < min_bic_tuple[0] and bic_tuple[1] > 0.8):
                    min_bic_tuple = bic_tuple
                    EMSPs_selection.send(True)
                    step_count += 1
                else:
                    EMSPs_selection.send(False)

        bic_tuple = min_bic_tuple

    # Save optimal results and calculate relative proportions
    r_squared, EMSPs, estimated_coefficients = bic_tuple[1], bic_tuple[2], bic_tuple[3]
    xtotal = estimated_coefficients.sum()

    for n, estimated_coefficient in enumerate(estimated_coefficients):
        EMSPs[n].EMSP_abs = estimated_coefficient
        EMSPs[n].EMSP_rel = estimated_coefficient / (1.0 * xtotal)

    EMSP_rel_array = round_percent([x.EMSP_rel for x in EMSPs], r_squared)
    
    for n, val in enumerate(EMSP_rel_array):
        EMSPs[n].EMSP_rel = val

    results.r_squared = np.round(r_squared, 2)
    results.bic = bic_tuple[0]

    return results, EMSPs

class EMSPs_Optimizer:
    """
    A class for optimizing EMSPs.
    """

    def __init__(self, EMSPs, decomposition_window):
        self.EMSPs = EMSPs
        self.decomposition_window = decomposition_window

    def select_EMSPs(self, available_mutation_types, mode, random_seed=1):
        """
        Generator function to create EMSP sets based on the given mode.
        """
        random.seed(random_seed)

        if mode == "random":
            while True:
                # Randomly remove some EMSPs and generate a coefficient matrix
                thinning_mutation_types = set(random.choices([v for v in available_mutation_types if v not in [-1, 0, 1]], k=len(available_mutation_types) - <Thinning Number>))
                selected_mutation_types = list(set(available_mutation_types) - thinning_mutation_types)
                self.sEMSPs = select_EMSPs(self.EMSPs, selected_mutation_types)
                self.coefficient_matrix =  generate_coefficient_matrix(self.sEMSPs, self.decomposition_window)
                yield self.sEMSPs, self.coefficient_matrix

        elif mode == "backstep":
            while True:
                # Remove a portion of EMSPs and evaluate if the change is acceptable
                thinning_mutation_types = random.sample([v for v in available_mutation_types if v not in [-1, 0, 1]], k=int(len(available_mutation_types) * 0.1))
                selected_mutation_types = [v for v in available_mutation_types if v not in thinning_mutation_types]
                self.sEMSPs = select_EMSPs(self.EMSPs, selected_mutation_types)
                self.coefficient_matrix = generate_coefficient_matrix(sEMSPs, self.decomposition_window)
                
                can_accept = yield self.sEMSPs, self.coefficient_matrix
                if can_accept:
                    available_mutation_types = selected_mutation_types

```

## Server implementation

PtWAVE realizes the front-end interface using Streamlit, and the back-end server is built with FastAPI and customized Python scripts. 

- Front-end

```python
import streamlit as st
import requests

def main():

    input_form, output_form = st.tabs(["Input parameters", "Output"])

    with input_form:

        <Setting input widgets of Streamlit>
        ...
        # Start analysis of input data and user parameters when the "Analyze" button is clicked.
        if st.button(label="Analyze"):
            on_button_click(input_data, user_parameters)
    with output_form:

        output_path = st.session_state.get("output", None)
        analysis_output = functions.load_json(f"{output_path}/output.json")

        <Writing output images using Streamlit>
        ...

def on_button_click(input_data, user_parameters):
    
    abi_control, abi_edited = functions.file_upload(abi_control_file, abi_edited_file, directory_path)

    # POST
    try:
        response = requests.post(API_URL,
            json={
                <Data Key>: <Input Data>
                }
            )
        response_data = response.json()

        <Checking response>
        ...
    
    except Exception as e:
        <Error handling>
        ...
```

- Back-end

```python
from fastapi import FastAPI

app = FastAPI()

@app.post("/analyze/")
async def analyze(params: InputParameters):
    data = await backend_analyze(params)
    return data

async def backend_analyze(params):
    try:
        ptwave(input_data, user_parameters)
        
        return {'status': 'success', 'output_directory': OUTPUT_DIRECTORY}

    except Exception as e:
        <Error handling>
        ...
```