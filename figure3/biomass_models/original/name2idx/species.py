from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "EGF",
    "EGFR",
    "Ra",
    "R2",
    "RP",
    "PLCg",
    "RPL",
    "RPLP",
    "PLCgP",
    "Grb2",
    "RG",
    "SOS",
    "RGS",
    "GS",
    "Shc",
    "RSh",
    "RShP",
    "ShP",
    "RShG",
    "ShG",
    "RShGS",
    "ShGS",
    "PLCgP_I",
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
