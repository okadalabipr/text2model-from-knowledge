import glob
import json
import os
import pickle
import sys
import tarfile
import xml.etree.cElementTree as ET

import gilda
import numpy as np
import spacy


def _gilda_ground(query: str, **kwargs) -> str or None:
    matches = gilda.ground(query, **kwargs)
    if matches:
        return "|".join(
            [matches[0].term.entry_name, matches[0].term.db, matches[0].term.id]
        )
    else:
        # print("No matches were found.")
        return None


if __name__ == "__main__":
    datapath = "./data"  # location of pubtator tar files

    # load spacy model
    nlp = spacy.load("en_core_sci_sm")

    # get list of tar files to process
    tar_list = os.listdir(datapath)

    # process each tar file
    for tar in tar_list:
        # get list of xml files in tar file
        with tarfile.open(os.path.join(datapath, tar), "r") as f:
            xml_list = f.getnames()
        tarpath = os.path.join(datapath, tar)
        for xml_name in xml_list:
            basename = os.path.basename(xml_name)
            outpath = "./pubtator_processed"

            # os.environ["PYSTOW_HOME"] = "./pystow_data"

            with tarfile.open(tarpath, "r") as f:
                with f.extractfile(xml_name) as biocxml:
                    try:
                        root = ET.parse(biocxml).getroot()
                    except ET.ParseError:
                        print("Bad XML file")
                        sys.exit(1)

            ids = []
            titles = []
            journals = []
            sentence_list = []
            entity_list = []
            span_list = []
            for document in root.findall("./document"):
                id = document.find("./id").text
                # check if pmcid if the id is a pmcid
                pmcid = document.find("./passage/infon[@key='article-id_pmc']")
                if pmcid:
                    pmcid = pmcid.text
                if pmcid == id:
                    id = "PMC" + id
                title = document.find("./passage/text").text
                journal = document.find("./passage/infon[@key='journal']")
                if journal is not None:
                    journal = journal.text.split(";")[0]
                else:
                    journal = "n/a"

                sentences = []
                entities = []
                spans = []
                for passage in document.findall("./passage"):
                    locations = []
                    annotations = []
                    lengths = []
                    for annotation in passage.findall("./annotation"):
                        loc = int(annotation.find("./location").attrib["offset"])
                        length = int(annotation.find("./location").attrib["length"])
                        name = annotation.find("./text").text
                        locations.append(loc)
                        annotations.append(name)
                        lengths.append(length)
                    locations = np.array(locations)
                    if (
                        len(annotations) > 0
                    ):  # checking if there are any annotations in the passage
                        passage_offset = int(passage.find("./offset").text)
                        text = passage.find("./text").text
                        doc = nlp(text)
                        for sent in doc.sents:
                            offset = sent.start_char + passage_offset
                            length = len(sent.text)
                            anno_inside = list(
                                *np.where(
                                    (locations >= offset)
                                    & (locations < offset + length)
                                )
                            )
                            if (
                                len(anno_inside) > 0
                            ):  # checking if there are any annotations in the sentence
                                normalized_names = []
                                in_spans = []
                                for idx in anno_inside:
                                    normalized = _gilda_ground(
                                        annotations[idx], context=sent.text
                                    )
                                    if normalized:
                                        normalized_names.append(normalized)
                                        start = locations[idx] - offset
                                        end = start + lengths[idx]
                                        in_spans.append((start, end))
                                if (
                                    len(normalized_names) > 0
                                ):  # checking if the normalization was successful
                                    sentences.append(sent.text)
                                    entities.append(normalized_names)
                                    spans.append(in_spans)
                ids.append(id)
                titles.append(title)
                journals.append(journal)
                sentence_list.append(sentences)
                entity_list.append(entities)
                span_list.append(spans)

            if not os.path.exists(outpath):
                os.mkdir(outpath)

            if len(ids) != 0:
                # only save pmid and entities for co-occurrence analysis
                with open(os.path.join(outpath, basename + "_id.pkl"), "wb") as output:
                    pickle.dump(ids, output)
                # with open(
                #     os.path.join(outpath, basename + "_title.pkl"), "wb"
                # ) as output:
                #     pickle.dump(titles, output)
                # with open(
                #     os.path.join(outpath, basename + "_journal.pkl"), "wb"
                # ) as output:
                #     pickle.dump(journals, output)
                # with open(
                #     os.path.join(outpath, basename + "_sentence.pkl"), "wb"
                # ) as output:
                #     pickle.dump(sentence_list, output)
                with open(
                    os.path.join(outpath, basename + "_entity.pkl"), "wb"
                ) as output:
                    pickle.dump(entity_list, output)
                # with open(
                #     os.path.join(outpath, basename + "_span.pkl"), "wb"
                # ) as output:
                #     pickle.dump(span_list, output)

    # combine pickled data to generate pubtator_data.json
    pmid_list = sorted(glob.glob(os.path.join(outpath, "*id.pkl")))
    # title_list = sorted(glob.glob(os.path.join(outpath, "*title.pkl")))
    # journal_list = sorted(glob.glob(os.path.join(outpath, "*journal.pkl")))
    # sentence_list = sorted(glob.glob(os.path.join(outpath, "*sentence.pkl")))
    # span_list = sorted(glob.glob(os.path.join(outpath, "*span.pkl")))
    entity_list = sorted(glob.glob(os.path.join(outpath, "*entity.pkl")))

    pmids_all = []
    for pmid in pmid_list:
        with open(pmid, "rb") as f:
            pmids = pickle.load(f)
        pmids_all.extend(pmids)

    # titles_all = []
    # for file in title_list:
    #     with open(file, "rb") as f:
    #         titles = pickle.load(f)
    #     titles_all.extend(titles)

    # journals_all = []
    # for file in journal_list:
    #     with open(file, "rb") as f:
    #         journals = pickle.load(f)
    #     journals_all.extend(journals)

    # sentences_all = []
    # for file in sentence_list:
    #     with open(file, "rb") as f:
    #         sentences = pickle.load(f)
    #     sentences_all.extend(sentences)

    # spans_all = []
    # for file in span_list:
    #     with open(file, "rb") as f:
    #         spans = pickle.load(f)
    #     spans_all.extend(spans)

    entities_all = []
    for file in entity_list:
        with open(file, "rb") as f:
            entities = pickle.load(f)
        entities_all.extend(entities)

    pubtator_data = []
    for pmid, entity in zip(pmids_all, entities_all):
        if len(entity) >= 2:  # set minimum number of annotated sentences in an article
            data_dict = {
                "pmid": pmid,
                # "title": title,
                # "journal": journal,
                # "sentences": sentence,
                # "spans": spans,
                "entities": entity,
            }
            pubtator_data.append(data_dict)

    # dump the data into a json file
    with open("pubtator_data.json", "w") as f:
        json.dump(pubtator_data, f)
