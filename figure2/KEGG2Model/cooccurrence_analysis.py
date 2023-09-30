import gzip
import json
import os
from typing import List, Literal, Optional, Tuple, Union

import gilda
import networkx as nx
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer


class CoocAnalyzer:
    def __init__(self, datapath: str) -> None:
        if os.path.exists(datapath):
            self.datapath = datapath
        else:
            raise FileNotFoundError("Given data directory was not found.")

    def load_pubtator_data(self) -> None:
        if self.datapath.endswith(".gz"):
            with gzip.open(self.datapath, "rb") as f:
                self.pubtator_data = json.load(f)
        else:
            with open(self.datapath, "r") as f:
                self.pubtator_data = json.load(f)

    def _dummy(self, tokens: List[str]) -> List[str]:
        return tokens

    def calculate_cooc(self, set_diag_0: bool = True) -> None:
        self.cv = CountVectorizer(preprocessor=self._dummy, tokenizer=self._dummy)
        self.word_count = self.cv.fit_transform(
            [
                entities
                for data in self.pubtator_data
                for entities in data["entities"]
                if len(entities) > 1
            ]
        )
        self.cooc_matrix = self.word_count.T * self.word_count
        if set_diag_0 is True:
            self.cooc_matrix.setdiag(0)
        self.network = nx.from_numpy_array(self.cooc_matrix.toarray())
        labels = {}
        for idx, label in enumerate(self.cv.get_feature_names_out()):
            labels[idx] = label
        self.network = nx.relabel_nodes(self.network, labels)

    def _analyze_sub_cooc(self, query_list: List[str], method: Literal["AND", "OR"]):
        # returns results of co-occurrence analysis on a subset of articles
        # that have been filtered by the query
        grounded_list = []
        for query in query_list:
            if len(query.split("|")) == 3:
                grounded = query
            else:
                grounded = _gilda_ground(query)
                if grounded is None:
                    print(f"Ignoring {query}: Invalid search query.")
                    continue
            grounded_list.append(grounded)
        if len(grounded_list) == 0:
            print("No valid search query.")
            return

        query_set = set(grounded_list)

        def operand(query_set: set, article_set: set):
            if method == "AND":
                return query_set.issubset(article_set)
            elif method == "OR":
                return not query_set.isdisjoint(article_set)

        # obtain a list of articles relevant to the query
        articles = []
        for idx, data in enumerate(self.pubtator_data):
            ent_set = set([])
            for entities in data["entities"]:
                ent_set |= set(entities)
            if operand(query_set, ent_set):
                articles.append(idx)
        # calculate the sub-co-occurrence network with the article list
        filtered_data = [
            entities
            for article in articles
            for entities in self.pubtator_data[article]["entities"]
        ]
        cv = CountVectorizer(preprocessor=self._dummy, tokenizer=self._dummy)
        word_count = cv.fit_transform(filtered_data)
        entity_list = list(cv.get_feature_names_out())
        count_vec = np.squeeze(np.asarray(word_count.sum(0)))
        entity_counts = dict(zip(entity_list, count_vec))
        sub_cooc_mtx = word_count.T * word_count
        sub_cooc_mtx.setdiag(0)
        sub_network = nx.from_numpy_array(sub_cooc_mtx.toarray())
        labels = {}
        for idx, label in enumerate(entity_list):
            labels[idx] = label
        return entity_counts, nx.relabel_nodes(sub_network, labels)

    def _fetch_weight_from_network(
        self, source: str, target: str, network: Union[nx.DiGraph, None] = None
    ) -> Tuple[int, bool]:
        if not network:
            network = self.network
        try:
            weight = network[source][target]["weight"]
            exist = True
        except KeyError:
            weight = 1
            exist = False
        return weight, exist


def _gilda_ground(query: str, **kwargs) -> Optional[str]:
    if len(query.split("|")) == 3:
        return query

    matches = gilda.ground(query, **kwargs)
    if matches:
        return "|".join(
            [matches[0].term.entry_name, matches[0].term.db, matches[0].term.id]
        )
    else:
        return None
