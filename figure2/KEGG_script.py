import os
import urllib.request

from KEGG2Model.convert_network import to_text2model
from KEGG2Model.coocurrence_analysis import CoocAnalyzer
from KEGG2Model.KGML_parser import KEGGpathway
from KEGG2Model.weight_visualizer import weightVis

# download KGML file of Human ErbB signaling pathway (hsa04012)
if not os.path.exists("hsa04012.xml"):
    url = "http://rest.kegg.jp/get/hsa04012/kgml"
    res = urllib.request.urlopen(url)
    with open("hsa04012.xml", "wb") as f:
        f.write(res.read())

outpath = "out"
if not os.path.exists(outpath):
    os.mkdir(outpath)

# process & visualize KGML file
p = KEGGpathway("hsa04012.xml")
p.visualize(os.path.join(outpath, "hsa04012.html"), show=False, informative=False)

# convert KGML file to Text2Model
graph, reactions = to_text2model(p)
# output Text2Model file
with open(os.path.join(outpath, "hsa04012_text2model.txt"), "w") as f:
    f.write(reactions)

# load pubtator data
pubtator = CoocAnalyzer(os.path.join("pubtator", "pubtator_data.json.gz"))
pubtator.load_pubtator_data()
pubtator.calculate_cooc()

# map text-mined information to the pathway and highlight components related to queries
vis = weightVis()
p.add_weights(pubtator)
queries = ["MCF-7", "MDA-MB-231"]
for query in queries:
    highlighted = p.highlight([query], method="AND", cooc_data=pubtator)
    vis.draw(
        highlighted,
        file_name=os.path.join(outpath, f"hsa04012_{query}.html"),
        scale=3,
        method="size",
        remove_lonely=True,
    )
