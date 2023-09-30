from typing import Tuple

import networkx as nx

from .KGML_parser import KEGGpathway


def to_text2model(pathway: KEGGpathway) -> Tuple[nx.DiGraph, str]:
    output = pathway.graph.copy()
    _reassign_group_edges(pathway, output)
    output.remove_nodes_from(list(nx.isolates(output)))  # remove lonely nodes
    network = output.copy()
    remove_edges = set()
    for (source, target), edge_dict in network.edges.items():
        if edge_dict["relation"] == "GErel":
            # do nothing
            pass
        else:
            if edge_dict["type"] == "ubiquitination":
                new_node = "u_" + target
                if new_node not in output:
                    level = network.nodes[target]["level"] + 60 / 150
                    output.add_node(new_node, type="intermediate", level=level)
                    output.add_edge(target, new_node, type="transition")
            elif (edge_dict["type"] == "dephosphorylation") or (
                edge_dict["type"] == "inhibition"
            ):
                new_node = "a_" + target
                if new_node not in output:
                    level = network.nodes[target]["level"] + 60 / 150
                    output.add_node(new_node, type="intermediate", level=level)
                # move original edge target to new node
                output.add_edge(source, new_node, **edge_dict)
                # add new transition node
                output.add_edge(new_node, target, type="transition")
                # delete original edge
                remove_edges.add((source, target))
                if edge_dict["effect"] == -1:
                    # move all outgoing edges from the target node to new node
                    for child in network[target]:
                        output.add_edge(new_node, child, **network[target][child])
                        # remove copied edge from the original node
                        remove_edges.add((target, child))
                continue
            elif edge_dict["type"] == "dissociation":
                # ignore for now
                pass
            else:
                new_node = "a_" + target
                if new_node not in output:
                    level = network.nodes[target]["level"] + 60 / 150
                    output.add_node(new_node, type="intermediate", level=level)
                    output.add_edge(target, new_node, type="transition")
                if edge_dict["effect"] == 1:
                    for child in network[target]:
                        output.add_edge(new_node, child, **network[target][child])
                        remove_edges.add((target, child))
    output.remove_edges_from(remove_edges)
    reactions = _list_reactions(output)
    return output, reactions


def _reassign_group_edges(pathway: KEGGpathway, output: nx.DiGraph):
    # reassigns edges that are to/from a node that belongs to a group
    groups = [
        key for key, value in pathway.graph.nodes.items() if value["type"] == "group"
    ]
    remove_edges = set()
    for group in groups:
        components = [
            pathway.entries[id][0] for id in pathway.graph.nodes[group]["components"]
        ]
        for comp in components:
            for source in pathway.graph.predecessors(comp):
                output.add_edge(source, group, **pathway.graph[source][comp])
                remove_edges.add((source, comp))
            for target in pathway.graph.successors(comp):
                output.add_edge(group, target, **pathway.graph[comp][target])
                remove_edges.add((comp, target))
    output.remove_edges_from(remove_edges)


def _list_reactions(network):
    sorted_nodes = _sort_nodes_by_level(network)
    sorted_edges = sorted(
        nx.line_graph(network), key=lambda x: sorted_nodes.index(x[0])
    )
    reactions = ""
    for source, target in sorted_edges:
        edge_dict = network[source][target]
        if edge_dict["type"] == "transition":
            continue
        elif edge_dict["type"] == "bind":
            reactions += f"{source} binds {target} <-> {source}_{target}\n"
        elif edge_dict["relation"] == "GErel":
            if edge_dict["effect"] == 1:
                reactions += f"{source} transcribes {target}\n"
            elif edge_dict["effect"] == -1:
                reactions += f"{source} degrades {target}\n"
        elif edge_dict["relation"] == "PPrel":
            if edge_dict["type"] == "phosphorylation":
                product = "a_" + target
                reactions += f"{source} phosphorylates {target} -> {product}\n"
                new_reaction = f"{product} is dephosphorylated -> {target}\n"
                if new_reaction not in reactions:
                    reactions += new_reaction
            elif edge_dict["type"] == "dephosphorylation":
                product = target.strip("*")
                reactions += f"{source} dephosphorylates {target} -> {product}\n"
            elif edge_dict["type"] == "ubiquitination":
                product = "u_" + target
                reactions += f"{source} ubiquitinates {target} -> {product}\n"
                reactions += f"{product} is degraded"
            elif edge_dict["type"] == "dissociation":
                continue
            elif (edge_dict["effect"] == 1) or (
                edge_dict["type"] == "binding/association"
            ):
                product = "a_" + target
                reactions += f"{source} activates {target} -> {product}\n"
                new_reaction = f"{product} is deactivated -> {target}\n"
                if new_reaction not in reactions:
                    reactions += new_reaction
            elif edge_dict["effect"] == -1:
                product = target.strip("a_")
                reactions += f"{source} deactivates {target} -> {product}\n"
    return reactions


def _sort_nodes_by_level(network):
    # sorts nodes in given network by their "level" value
    sorted_nodes = sorted(
        [(key, value) for key, value in network.nodes.items()],
        key=lambda x: x[1]["level"],
    )
    return [item[0] for item in sorted_nodes]
