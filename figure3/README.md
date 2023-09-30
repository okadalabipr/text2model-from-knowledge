# `figure3/`

This directory contains the data and code to reproduce Figure 3 in the paper (using LLMs to generate Text2Model files).

## About the files

### `biomass_models/`

This directory contains the biomass models used in the paper.

#### `original/` * `reconstructed/`

Includes the model files for the original model reported in [Kholodenko *et al*., 1999](https://doi.org/10.1074/jbc.274.42.30169), and the model reconstructed from the LLM's output, respectively.

Created using the Text2Model file which was manually created by the authors of the paper (`original.txt` & `reconstructed.txt`). The Text2Model file can be converted with the following code:

```python
from biomass import Text2Model

original = Text2Model("original.txt")
reconst = Text2Model("reconstructed.txt")

original.convert()
reconst.convert()
```

### `llm_prompt/`

This directory contains the LLM prompt used in the paper.

The original Text2Model output used in the paper was obtained through OpenAI's Playground interface on 2023-03-31.

The prompt used is documented in `original_prompt.txt`, and the output is documented in `original_completion.txt`. Furthermore, the parameter settings used in the Playground interface are documented in `openai_api_prompt.py` as a script that can be run to reproduce the results through OpenAI's API, however this was not the method used in the paper.

## Reproducing the results

The experimental data from Kholodenko *et al.*, 1999 that were used to fit the model parameters can be found in the `observable.py` file in both model directories:

```python
class Observable(DifferentialEquation):
    ...
    def set_data(self) -> None:
        self.experiments[self.obs_names.index("Total_phosphorylated_PLCg")] = {
            "EGF20nM": [0.0, 0.095, 0.055, 0.042, 0.038, 0.025],
            "EGF2nM": [0.0, 0.044, 0.048, 0.024, 0.013, 0.022]
        }
        ...

    def get_timepoint(self, obs_name: str) -> List[int]:
        if obs_name in self.obs_names:
            return [0, 15, 30, 45, 60, 120]
```

The code used to estimate the parameters and generate the plots can be found in the `/biomass_models/biomass_scipt.py` file as an executable script. Reproducing the results can take up to a couple of days depending on the available computational resources.
