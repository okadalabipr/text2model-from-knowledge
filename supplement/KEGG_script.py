import os
import urllib.request

from KEGG2Model.convert_network import to_text2model
from KEGG2Model.KGML_parser import KEGGpathway

# download KGML file of Human JAK-STAT signaling pathway (hsa04630)
if not os.path.exists("hsa04630.xml"):
    print("Downloading KGML file of Human JAK-STAT signaling pathway (hsa04630)")
    url = "http://rest.kegg.jp/get/hsa04630/kgml"
    res = urllib.request.urlopen(url)
    with open("hsa04630.xml", "wb") as f:
        f.write(res.read())

outpath = "out"
if not os.path.exists(outpath):
    os.mkdir(outpath)

# process & visualize KGML file
print("Processing KGML file of Human JAK-STAT signaling pathway (hsa04630)")
p = KEGGpathway("hsa04630.xml")
print("Writing visualization results to hsa04630.html")
p.visualize(os.path.join(outpath, "hsa04630.html"), show=False, informative=False)

# convert KGML file to Text2Model
print("Converting KGML file to Text2Model")
graph, reactions = to_text2model(p)
# output Text2Model file
print("Writing Text2Model file to hsa04630_text2model.txt")
with open(os.path.join(outpath, "hsa04630_text2model.txt"), "w") as f:
    f.write(reactions)
