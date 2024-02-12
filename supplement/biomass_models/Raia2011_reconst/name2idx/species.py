from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "IL13",
    "Rec",
    "IL13_Rec",
    "JAK2",
    "pJAK2",
    "SOCS3",
    "p_IL13_Rec",
    "STAT5",
    "pSTAT5",
    "CD274mRNA",
    "SOCS3mRNA",
    "DecoyR",
    "IL13_DecoyR",
    "SHP1",
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
