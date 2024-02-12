from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "CNTF_CTF1_CLCF1_IL27_IL6_IL11_IL13_IL31_LIF_OSM",
    "Rec",
    "IL13_Rec",
    "a_Rec",
    "JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2",
    "a_JAK1_JAK2_CNTFR_IL31RA_IL4R_IL6R_IL6ST_IL11RA_IL13RA1_IL13RA2_LIFR_OSMR_TYK2",
    "STAT",
    "a_STAT",
    "PTPN6",
    "SOCSmRNA",
    "SOCS4_SOCS7_SOCS1_SOCS2_SOCS3_SOCS6_SOCS5",
]

NUM: int = len(NAMES)

Species = make_dataclass(
    cls_name="Species",
    fields=[(name, int) for name in NAMES],
    namespace={"NAMES": NAMES, "NUM": NUM},
    frozen=True,
)

name2idx: Dict[str, int] = {k: v for v, k in enumerate(NAMES)}

V = Species(**name2idx)

del name2idx
