# [`supplement`](supplement/)

This directory contains the supplementary material for the paper, which applies the methods in the paper to a JAK-STAT pathway model.

## Running the pre-requisite script

Before running the script wihtin this directory, run the pre-requisite script `pre-requisite.sh` to prepare the environment:

```bash
cd supplement
bash pre-requisite.sh
```

## Modified version of BioMASS

Since the experimental data used in this analysis contained varying time points for different experimental conditions, several modifications to the BioMASS package as well as the model files were made. The pre-requisite script will clone the original BioMASS repository and copy the modified files to the appropriate locations.

The modifications are documented below.

### [`biomass/temporal_dynamics.py`](./biomass_modification/temporal_dynamics.py)

The modifications in this file were made in mehods `_plot_experimental_data_with_error_bars` (lines 345, 362) and `_plot_experimental_data_without_error_bars` (lines 397, 414), indicated with comments showing the changes made.

### [`biomass_models/*/problem.py`](./biomass_models/KEGG_JAK-STAT/problem.py)

A single modification (in line 39) was made in the `problem.py` file for each model.

### [`biomass_models/*/observable.py`](./biomass_models/KEGG_JAK-STAT/observable.py)

To allow for the varying time points for each experimental condition, the `set_data` and `get_timepoint` methods were defined using dictionaries, unlike in the models used in the main text.

```python
def set_data(self) -> None:
    ...
    self.experiments[self.obs_names.index("IL13_cell")] = {
            "IL13_4": [1.71, 2.08, 2.37, 2.37, 2.88, 3.03, 3.10, 3.24, 3.61, 3.79, 3.17],
            "IL13_20": [4.59, 4.74, 5.61, 5.54, 6.64, 6.96, 6.49, 5.69, 6.82, 6.60, 5.69],
    }...

def get_timepoint(self, obs_name: str) -> List[int]:
    if obs_name in self.obs_names:
        if obs_name == "RacSurf":
            return {
                "IL13_4": [0, 10, 20, 30, 60, 90, 120],
                }...
```

## About the files

### [`biomass_models/`](./biomass_models/)

This directory contains the biomass models used in the supplement.

#### [`KEGG_JAK-STAT/`](./biomass_models/KEGG_JAK-STAT/)

Includes the model files for the KEGG JAK-STAT model. The Text2Model file used to generate the model can be found in [`KEGG_JAK-STAT.txt`](./biomass_models/KEGG_JAK-STAT.txt).

#### [`Raia2011_original`](./biomass_models/Raia2011_original/)

Includes the manually-reproduced model from a previous study ([Raia *et al*., 2011](https://doi.org/10.1158/0008-5472.CAN-10-2987)). The corresponding Text2Model file can be found in [`Raia2011_original.txt`](./biomass_models/Raia2011_original.txt).

Note that the original model was not fully reproduced due to the differences in how the original model and Text2Model represent each reaction. Furthermore, parameter estimation was not conducted for the original model.

#### [`Raia2011_reconst`](./biomass_models/Raia2011_reconst/)

Includes the model reconstructed from the LLM's output. The corresponding Text2Model file used to generate the model can be found in [`Raia2011_reconst.txt`](./biomass_models/Raia2011_reconst.txt).

### [`llm_prompt/`](./llm_prompt/)

This directory contains the LLM prompt used in the supplement.

The original Text2Model output used in the supplement was obtained through OpenAI's Playground interface on 2023-04-07.

The prompt used is documented in [`original_prompt.txt`](./llm_prompt/original_prompt.txt), and the output is documented in [`original_completion.txt`](./llm_prompt/original_completion.txt). Furthermore, the parameter settings used in the Playground interface are documented in [`openai_api_prompt.py`](./llm_prompt/openai_api_prompt.py) as a script that can be run to reproduce the results through OpenAI's API, however this was not the method used in the supplement.

### [`KEGG_script.py`](./KEGG_script.py)

This script can be used to generate the network visualizations and original Text2Model file for the KEGG JAK-STAT model. The script will create and output the results under `out/`.
