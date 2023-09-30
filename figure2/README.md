# [`figure2/`](figure2/)

This directory contains the data and code to generate Figure 2 in the paper (generating executable models from KEGG PATHWAYS).

## About the files

### [`biomass_models/`](./biomass_models/)

This directory contains the biomass models used in the paper.

#### [`KEGG_erbb_MCF7`](./biomass_models/KEGG_erbb_MCF7/) & [`KEGG_erbb_MDAMB231`](./biomass_models/KEGG_erbb_MDAMB231/)

Includes the model files for the MCF7 and MDAMB231 models, respectively. Both models were generated using the same Text2Model file ([`KEGG_erbb_text2model.txt`](./biomass_models/KEGG_erbb_text2model.txt)), but fitted to different experimental data.

### [`KEGG2Model/`](./KEGG2Model/)

This directory contains the code necessary to convert KEGG PATHWAYS to Text2Model files, as well as conducting co-occurrence analysis on the PubTator data.

#### Contents

- [`convert_network.py`](./KEGG2Model/convert_network.py): Converting KEGG PATHWAYS to Text2Model files.
- [`cooccurrence_analysis`](./KEGG2Model/cooccurrence_analysis.py): Conducting co-occurrence analysis on the PubTator data.
- [`KGML_parser.py`](./KEGG2Model/KGML_parser.py): Parsing KGML files.
- [`weight_visualizer.py`](./KEGG2Model/weight_visualizer.py): Visualizing weights of networks.

### [`pubtator/`](./pubtator/)

This directory contains the PubTator data used in the paper.

Although the script to reproduce the results can be found as [`prepare_data.py`](./pubtator/prepare_data.py) and  [`process_data.py`](./pubtator/process_data.py), this process can be data-intensive and time-consuming. Furthermore, the PubTator data that is used in the paper is the Jan. 2022 release, which seems to be no longer available. Therefore, the processed data is included as a compressed json file ([`pubtator_data.json.gz`](./pubtator/pubtator_data.json.gz)).

## Reproducing the results

### Network Visualizations

The network visualizations can be reproduced using the [`KEGG_script.py`](./KEGG_script.py) script. The script will create and output the results under `out/`.

### Simulation Results

The experimental data used to estimate the model parameters were obtained from a previous study ([Imoto *et al.*, 2020](https://doi.org/10.3390/cancers12102878)). The data can be found in the `observable.py` file in both model directories:

```python
class Observable(DifferentialEquation):
    ...
    def set_data(self) -> None:
        self.experiments[self.obs_names.index("Phosphorylated_AKT")] = {
            "EGF": [0.0, 0.242, 0.087, 0.088, 0.082, 0.045, 0.017, 0.043],
            "HRG": [0.0, 0.976, 1.0, 0.96, 0.876, 0.836, 0.77, 0.719],
        }
        ...

    def get_timepoint(self, obs_name: str) -> List[int]:
        if obs_name in self.obs_names:
            return [0, 5, 15, 30, 45, 60, 90, 120]
```

The code used to estimate the parameters and generate the plots can be found in the [`/biomass_models/biomass_scipt.py`](./biomass_models/biomass_script.py) file as an executable script. Reproducing the results can take up to a couple of days depending on the available computational resources.
