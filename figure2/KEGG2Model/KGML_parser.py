import os
import re
from time import sleep
import urllib.request
import xml.etree.ElementTree as ET
from typing import List, Literal, Optional, Tuple

import famplex
import networkx as nx
import numpy as np
from pyvis.network import Network

from .cooccurrence_analysis import CoocAnalyzer, _gilda_ground


class KEGGpathway:
    def __init__(self, kgml_path: str) -> None:
        # visualization options
        self.node_shapes = {
            "gene": "dot",
            "compound": "hexagon",
            "group": "square",
            "map": "ellipse",
            "ortholog": "diamond",
        }
        self.edge_colors = {
            "activation": "darkorange",
            "expression": "palevioletred",
            "repression": "seagreen",
            "inhibition": "navy",
            "binding/association": "firebrick",
            "phosphorylation": "gold",
            "dephosphorylation": "burlywood",
            "glycosylation": "lightsteelblue",
            "ubiquitination": "darkslateblue",
            "methylation": "palegoldenrod",
            "indirect effect": "grey",
            "": "grey",
            "state change": "black",
            "dissociation": "black",
            "missing interaction": "black",
            "compound": "limegreen",
        }

        if os.path.exists(kgml_path):
            root = ET.parse(kgml_path).getroot()
            self._parse_kgml(root)
        else:
            raise FileNotFoundError("Given KGML file was not found.")

    def _parse_kgml(self, root) -> None:
        self.kegg_entries = {}
        self.pathwayID = root.attrib["name"]
        self.title = root.attrib["title"]
        # parse nodes
        self.entries = self._parse_nodes(root)

        self._ground_entries()

        self.graph = nx.DiGraph()
        self.graph.add_nodes_from([value for value in self.entries.values()])

        # parse edges
        self.relations, self.edges = self._parse_edges(root)

        self.graph.add_edges_from(self.edges)

    def _parse_nodes(self, root) -> dict:
        entries = {}
        for entry in root.findall("./entry"):
            idx = entry.attrib["id"]
            kegg_IDs = entry.attrib["name"].split(" ")
            link = "https://rest.kegg.jp/get/" + "+".join(kegg_IDs)
            entry_type = entry.attrib["type"]
            graphics = entry.find("./graphics")
            x = int(graphics.attrib["x"])
            y = int(graphics.attrib["y"])
            if entry_type == "gene":
                components = kegg_IDs
                comp_dict = self._get_component_dicts(kegg_IDs)
                name = "_".join([comp["name"] for comp in comp_dict])
                node_name = "_".join(
                    [
                        comp["symbols"][0]
                        for comp in comp_dict
                        if "symbols" in comp and comp["symbols"]
                    ]  # discard entries wihout symbols from node name
                )
            elif entry_type == "compound":
                components = kegg_IDs
                comp_dict = self._get_component_dicts(kegg_IDs)
                name = "_".join([comp["name"] for comp in comp_dict])
                node_name = name
            elif entry_type == "group":
                components = [
                    component.attrib["id"] for component in entry.findall("./component")
                ]
                node_name = "_".join(
                    [entries[component][0] for component in components]
                )
                link = "https://rest.kegg.jp/get/" + "+".join(
                    [
                        comp_id
                        for component in components
                        for comp_id in entries[component][1]["kegg_IDs"]
                    ]
                )
                name = "+".join(
                    [entries[component][1]["name"] for component in components]
                )
                # adding "_group" key to indicate entries that belong to a group
                # using "_group" since "group" is already used by pyvis (vis.Network)
                for component in components:
                    entries[component][1]["_group"] = idx
            else:
                node_name = graphics.attrib["name"]
                name = node_name
                components = []
            entry_attrib = {
                "entry_ID": idx,
                "kegg_IDs": kegg_IDs,
                "name": name,
                "type": entry_type,
                "components": components,
                "level": x / 150,  # magic number; for visualization purposes
                "x": x,
                "y": y,
                "ref": link,
            }
            entries[idx] = (node_name, entry_attrib)

        return entries

    def _parse_edges(self, root) -> Tuple[list, list]:
        relations = []
        edges = []
        for relation in root.findall("./relation"):
            source_id = relation.attrib["entry1"]
            target_id = relation.attrib["entry2"]
            source = self.entries[source_id][0]
            target = self.entries[target_id][0]
            rel = relation.attrib["type"]
            relation_dict = {
                "relation": rel,
                "_source": source,
                "_target": target,
                "effect": 1,
                "indirect": False,
                "type": "",
            }

            subtypes = []
            for subtype in relation.findall("./subtype"):
                name = subtype.attrib["name"]
                value = subtype.attrib["value"]
                subtype_dict = {
                    "name": name,
                    "value": value,
                }
                if (name == "compound") or (name == "hidden compound"):
                    relation_dict["effect"] = 1
                    relation_dict["type"] = name
                    relation_dict["compound"] = value
                elif (name == "activation") or (name == "expression"):
                    relation_dict["effect"] = 1
                    relation_dict["type"] = name
                elif (name == "inhibition") or (name == "repression"):
                    relation_dict["effect"] = -1
                    relation_dict["type"] = name
                elif name == "indirect effect":
                    relation_dict["effect"] = 1
                    relation_dict["type"] = name
                    relation_dict["indirect"] = True
                elif (
                    (name == "state change")
                    or (name == "binding/association")
                    or (name == "phosphorylation")
                    or (name == "dephosphorylation")
                    or (name == "glycosylation")
                    or (name == "ubiquitination")
                    or (name == "methylation")
                ):
                    relation_dict["type"] = name
                elif (name == "dissociation") or (name == "missing interaction"):
                    relation_dict["effect"] = 1
                    relation_dict["type"] = name
                subtypes.append(subtype_dict)

            relation_dict["subtypes"] = subtypes
            relations.append((source_id, target_id, relation_dict))
            edges.append((source, target, relation_dict))
        return relations, edges

    def _get_component_dicts(self, kegg_IDs: List[str]) -> List[dict]:
        to_process = []
        for kegg_ID in kegg_IDs:
            kegg_ID_split = kegg_ID.rsplit(":")[1]
            if kegg_ID_split not in self.kegg_entries:
                to_process.append(kegg_ID)
        if to_process:
            new_entries = self._kegg_rest_get(to_process)
            for entry in new_entries:
                self.kegg_entries[entry["KEGGID"]] = entry
        return [self.kegg_entries[kegg_ID.rsplit(":")[1]] for kegg_ID in kegg_IDs]

    def _kegg_rest_get(
        self,
        kegg_IDs: List[str],
        batch_size: int = 10,
        retries: int = 3,
        delay: int = 3,
    ) -> List[dict]:
        output = []
        for i in range(0, len(kegg_IDs), batch_size):
            id_list = kegg_IDs[i : i + batch_size]
            url = "https://rest.kegg.jp/get/" + "+".join(id_list)
            for _ in range(retries):
                try:
                    response = urllib.request.urlopen(url)
                    break
                except urllib.error.URLERROR as e:
                    sleep(delay)
            else:
                raise e
            res = response.read().decode()

            lines = iter(res.splitlines())
            for line in lines:
                line = line.strip("\n")
                if line.startswith("ENTRY"):
                    entry_id = re.search(r"ENTRY\s+(\w+)\s", line).group(1)
                    entry_dict = {"KEGGID": entry_id}

                elif line.startswith("SYMBOL"):
                    symbol = re.sub(r"SYMBOL\s+", "", line).split(", ")
                    entry_dict["symbols"] = symbol

                elif line.startswith("NAME"):
                    name = re.search(r"NAME\s+(.+)", line).group(1).strip(";")
                    entry_dict["name"] = name

                elif line.startswith("DBLINKS"):
                    db_links = {}
                    db_link = re.search(r"DBLINKS\s+(.+):\s(\w+)", line)
                    db_name = db_link.group(1)
                    db_id = db_link.group(2)
                    db_links[db_name] = db_id

                    line = next(lines).strip("\n")
                    while re.match(r"\s+", line):
                        db_link = re.search(r"\s+(.+):\s+(\w+)", line)
                        db_name = db_link.group(1)
                        db_id = db_link.group(2)
                        db_links[db_name] = db_id
                        line = next(lines).strip("\n")
                    entry_dict["db_links"] = db_links

                if line.startswith("///"):
                    output.append(entry_dict)

        return output

    def visualize(
        self,
        out_path: str,
        show: bool = True,
        network: Optional[nx.DiGraph] = None,
        remove_lonely: bool = True,
        informative: bool = True,
        height: str = "1024px",
        width: str = "1024px",
        directed: bool = True,
        layout: bool = True,
        **kwargs,
    ) -> None:
        nt = Network(height, width, directed=directed, layout=layout, **kwargs)
        graph = self.graph.copy()
        if remove_lonely:
            graph.remove_nodes_from(list(nx.isolates(self.graph)))
        if not network:
            nt.from_nx(graph)
        else:
            nt.from_nx(network)
        for node in nt.nodes:
            if informative:
                node["shape"] = self.node_shapes[node["type"]]
            node["title"] = ""
            for key in ["entry_ID", "name", "type", "ref"]:
                if key in node:
                    node["title"] += f"{key.capitalize()}: {node[key]}\n"
            if node["type"] == "group":
                node["title"] += f"Group components: {node['components']}"
        for edge in nt.edges:
            edge["title"] = ""
            for key in ["from", "to", "relation", "type", "effect", "indirect"]:
                if key in edge:
                    edge["title"] += f"{key.capitalize()}: {edge[key]}\n"
            if informative:
                if "compound" in edge:
                    edge["color"] = self.edge_colors["compound"]
                    edge[
                        "title"
                    ] += f"Compound: {self.entries[edge['compound']][1]['name']}\n"
                else:
                    edge["color"] = self.edge_colors[edge["type"]]
                if ("indirect" in edge) and (edge["indirect"]):
                    edge["dashes"] = [5, 5]
                if ("effect" in edge) and (edge["effect"] == -1):
                    edge["arrows"] = {"to": {"enabled": True, "type": "bar"}}
        if show:
            nt.show(out_path)
        else:
            nt.save_graph(out_path)

    def visualize_kgml(
        self,
        out_path: str,
        height: str = "1024px",
        width: str = "1024px",
        directed: bool = True,
        layout: bool = True,
        **kwargs,
    ) -> None:
        nt = Network(height, width, directed=directed, layout=layout, **kwargs)
        for key, (label, node_dict) in self.entries.items():
            nt.add_node(key, label=label, **node_dict)
        for source_id, target_id, relation_dict in self.relations:
            nt.add_edge(source_id, target_id, **relation_dict)

        for node in nt.nodes:
            node["size"] = 10
            if node["type"] == "gene":
                node["shape"] = "dot"
            elif node["type"] == "compound":
                node["shape"] = "hexagon"
            elif node["type"] == "group":
                node["shape"] = "square"
            elif node["type"] == "map":
                node["shape"] = "ellipse"
            elif node["type"] == "ortholog":
                node["shape"] = "diamond"
            node["title"] = ""
            for key in ["entry_ID", "name", "type", "ref"]:
                if key in node:
                    node["title"] += f"{key.capitalize()}: {node[key]}\n"
            if node["type"] == "group":
                node["title"] += f"Group components: {node['components']}"
        for edge in nt.edges:
            edge["title"] = ""
            for key in ["from", "to", "relation", "type", "effect", "indirect"]:
                if key in edge:
                    edge["title"] += f"{key.capitalize()}: {edge[key]}\n"
            if "compound" in edge:
                edge["color"] = "limegreen"
                edge[
                    "title"
                ] += f"Compound: {self.entries[edge['compound']][1]['name']}\n"
            elif edge["type"] == "activation":
                edge["color"] = "darkorange"
            elif edge["type"] == "expression":
                edge["color"] = "palevioletred"
            elif edge["type"] == "repression":
                edge["color"] = "seagreen"
            elif edge["type"] == "inhibition":
                edge["color"] = "navy"
            elif edge["type"] == "binding/association":
                edge["color"] = "firebrick"
            elif edge["type"] == "phosphorylation":
                edge["color"] = "gold"
            elif edge["type"] == "dephosphorylation":
                edge["color"] = "burlywood"
            elif edge["type"] == "glycosylation":
                edge["color"] = "lightsteelblue"
            elif edge["type"] == "ubiquitination":
                edge["color"] = "darkslateblue"
            elif edge["type"] == "methylation":
                edge["color"] = "palegoldenrod"
            elif edge["type"] == "":
                edge["color"] = "grey"
            if ("indirect" in edge) and (edge["indirect"]):
                edge["dashes"] = [5, 5]
            if ("effect" in edge) and (edge["effect"] == -1):
                edge["arrows"] = {"to": {"enabled": True, "type": "bar"}}
        nt.show(out_path)

    def _ground_entries(self):
        for key, entry in self.entries.items():
            normalized_list = []
            for component in entry[1]["components"]:
                grounded = _gilda_ground(component, namespaces=["HGNC"])
                if grounded:
                    normalized_list.append(grounded)
            entry[1]["normalized"] = normalized_list
            # use famplex to find an entity for entries with more than two components
            if (len(entry[1]["components"]) > 1) and (len(normalized_list) > 0):
                initial = (
                    normalized_list[0].split("|")[1],
                    normalized_list[0].split("|")[0],
                )
                if famplex.in_famplex(*initial):
                    parent = famplex.parent_terms(*initial)[0]
                    flag = True
                for norm in normalized_list:
                    norm_term = (norm.split("|")[1], norm.split("|")[0])
                    if norm_term not in famplex.child_terms(*parent):
                        flag = False
                if flag:
                    entry[1]["parent"] = parent[1]
                    self.entries[key] = (parent[1], entry[1])

    def add_weights(self, cooc_data: CoocAnalyzer):
        # adding node weights (word count)
        entity_list = list(cooc_data.cv.get_feature_names_out())
        count_vec = np.squeeze(np.asarray(cooc_data.word_count.sum(0)))
        weight_sum = 0
        for node in self.graph.nodes.keys():
            components = self._fetch_components(node)
            weight = 0
            counter = 0
            for node_comp in components:
                if node_comp in entity_list:
                    idx = entity_list.index(node_comp)
                    weight += count_vec[idx]
                    counter += 1
            if weight != 0:
                # adding information to "weight" so the node size does not change
                self.graph.nodes[node]["weight"] = weight / counter
                self.graph.nodes[node]["exist"] = True
            else:
                self.graph.nodes[node]["weight"] = weight
                self.graph.nodes[node]["exist"] = False
            count = self.graph.nodes[node]["weight"]
            weight_sum += count
        # normalizing weights
        for node_values in self.graph.nodes.values():
            node_values["weight"] /= weight_sum
        # adding edge weights
        for source, target in self.graph.edges.keys():
            source_comps = self._fetch_components(source)
            target_comps = self._fetch_components(target)
            weight = 0
            counter = 0
            exist = False
            for s_comp in source_comps:
                for t_comp in target_comps:
                    wei, exi = cooc_data._fetch_weight_from_network(s_comp, t_comp)
                    weight += wei
                    exist = exist | exi
                    if exi:
                        counter += 1
            if exist:
                self.graph[source][target]["weight"] = (
                    weight / counter
                )  # averaging weights
            else:
                self.graph[source][target]["weight"] = weight
            self.graph[source][target]["exist"] = exist

    def _fetch_components(self, node) -> List[str]:
        # fetch names of the components of a given node
        components = []
        if self.graph.nodes[node]["type"] == "group":
            for sub in self.graph.nodes[node]["components"]:
                components.extend(self.entries[sub][1]["normalized"])
                if "parent" in self.entries[sub][1]:
                    components.append(
                        _gilda_ground(
                            self.entries[sub][1]["parent"], namespaces=["FPLX"]
                        )
                    )
        else:
            components.extend(self.graph.nodes[node]["normalized"])
            if "parent" in self.graph.nodes[node]:
                components.append(
                    _gilda_ground(self.graph.nodes[node]["parent"], namespaces=["FPLX"])
                )
        return components

    def highlight(
        self,
        query_list: List[str],
        method: Literal["AND", "OR"],
        cooc_data: CoocAnalyzer,
    ):
        pathway = self.graph.copy()
        # adding edge weights
        sub_counts, subG = cooc_data._analyze_sub_cooc(query_list, method)
        for source, target in pathway.edges.keys():
            source_comps = self._fetch_components(source)
            target_comps = self._fetch_components(target)
            weight = 0
            counter = 0
            exist = False
            for s_comp in source_comps:
                for t_comp in target_comps:
                    wei, exi = cooc_data._fetch_weight_from_network(
                        s_comp, t_comp, subG
                    )
                    weight += wei
                    exist = exist | exi
                    if exi:
                        counter += 1
            if exist:
                pathway[source][target]["weight"] = (
                    weight / counter
                )  # averaging weights
            else:
                pathway[source][target]["weight"] = weight
            pathway[source][target]["exist"] = exist
        weight_sum = nx.adjacency_matrix(pathway).sum()
        for source, target in pathway.edges.keys():
            pathway[source][target]["weight"] /= self.graph[source][target]["weight"]

        # adding node weights
        weight_sum = 0
        for node in pathway.nodes.keys():
            components = self._fetch_components(node)
            weight = 0
            counter = 0
            for node_comp in components:
                if node_comp in sub_counts:
                    weight += sub_counts[node_comp]
                    counter += 1
            if weight != 0:
                pathway.nodes[node]["weight"] = weight / counter
                pathway.nodes[node]["exist"] = True
            else:
                pathway.nodes[node]["weight"] = weight
                pathway.nodes[node]["weight"] = False
            count = pathway.nodes[node]["weight"]
            weight_sum += count
        # normalizing & highlighting sizes
        for node in pathway.nodes.keys():
            pathway.nodes[node]["weight"] /= weight_sum
            if self.graph.nodes[node]["exist"]:
                pathway.nodes[node]["weight"] /= self.graph.nodes[node]["weight"]

        return pathway
