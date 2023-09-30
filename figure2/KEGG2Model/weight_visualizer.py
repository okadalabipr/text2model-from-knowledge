from typing import Literal

import numpy as np
import networkx as nx
from matplotlib.colors import ListedColormap, rgb2hex
from pyvis.network import Network


class weightVis:
    def __init__(
        self,
        up_color: str = "#ff7f50",
        down_color: str = "#4169e1",
        mid_color: str = "#f5f5f5",
        none_color: str = "#808080",
        color_n: int = 2,
    ):
        self.transform = lambda weight: np.power(weight, color_n) / (
            np.power(weight, color_n) + 1
        )
        self.up_color = up_color
        self.down_color = down_color
        self.mid_color = mid_color
        self.none_color = none_color
        self.N = 256  # number of steps in the color map / 2
        N_2 = int(self.N / 2)
        mid_r, mid_g, mid_b = self.hex2rgb(self.mid_color)
        top = np.ones((N_2, 4))
        up_r, up_g, up_b = self.hex2rgb(self.up_color)
        top[:, 0] = np.linspace(mid_r / 255, up_r / 255, N_2)
        top[:, 1] = np.linspace(mid_g / 255, up_g / 255, N_2)
        top[:, 2] = np.linspace(mid_b / 255, up_b / 255, N_2)
        btm = np.ones((N_2, 4))
        btm_r, btm_g, btm_b = self.hex2rgb(self.down_color)
        btm[:, 0] = np.linspace(btm_r / 255, mid_r / 255, N_2)
        btm[:, 1] = np.linspace(btm_g / 255, mid_g / 255, N_2)
        btm[:, 2] = np.linspace(btm_b / 255, mid_b / 255, N_2)
        colors = np.vstack((btm, top))
        self.color_map = ListedColormap(colors)

    def draw(
        self,
        network: nx.DiGraph,
        file_name: str,
        method: Literal["size", "color"] = "size",
        scale: int = 10,
        remove_lonely: bool = True,
        directed: bool = True,
        layout: bool = True,
        show: bool = False,
        **kwargs,
    ):
        nt = Network("1024px", "1024px", directed=directed, layout=layout, **kwargs)
        graph = network.copy()
        if remove_lonely:
            graph.remove_nodes_from(list(nx.isolates(network)))
        nt.from_nx(graph, show_edge_weights=True)
        for node in nt.nodes:
            node["title"] = f"weight: {node['weight']}"
            if not node["exist"]:
                node["size"] = 0.7 * scale
                node["color"] = {
                    "border": self.none_color,
                    "background": "white",
                    "highlight": {"border": self.none_color, "background": "white"},
                    "hover": {"border": self.none_color, "background": "white"},
                }
                node["shapeProperties"] = {"borderDashes": [3, 3]}
                node["title"] = "Node not found"
            elif method == "color":
                node["color"] = rgb2hex(self.color_map(self.transform(node["weight"])))
            elif method == "size":
                node["size"] = node["weight"] * scale
                # node["size"] = scale + np.log2(node["weight"])
                node["color"] = (
                    self.down_color if node["weight"] <= 1 else self.up_color
                )
        for edge in nt.edges:
            edge["title"] = f"weight: {edge['weight']}"
            if not edge["exist"]:
                edge["color"] = self.none_color
                edge["width"] = 0.7
                edge["dashes"] = [3, 3]
                edge["title"] = "Edge not found"
            elif method == "color":
                edge["color"] = rgb2hex(self.color_map(self.transform(edge["weight"])))
            elif method == "size":
                edge["width"] = edge["weight"]
                edge["color"] = (
                    self.down_color if edge["weight"] <= 1 else self.up_color
                )
        # nt.show_buttons(filter_=["physics", "nodes"])
        if show:
            nt.show(file_name)
        else:
            nt.save_graph(file_name)

    def hex2rgb(self, hex: str):
        h = hex.lstrip("#")
        return tuple(int(h[i : i + 2], 16) for i in (0, 2, 4))
