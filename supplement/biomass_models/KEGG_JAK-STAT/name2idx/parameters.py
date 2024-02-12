from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "kf27",
    "kr27",
    "kf28",
    "kr28",
    "V29",
    "K29",
    "kf30",
    "V31",
    "K31",
    "V32",
    "K32",
    "V33",
    "K33",
    "V34",
    "K34",
    "V36",
    "K36",
    "n36",
    "kf37",
    "kf38",
    "V39",
    "K39",
]

NUM: int = len(NAMES)

Parameters = make_dataclass(
    cls_name="Parameters",
    fields=[(name, int) for name in NAMES],
    namespace={"NAMES": NAMES, "NUM": NUM},
    frozen=True,
)

name2idx: Dict[str, int] = {k: v for v, k in enumerate(NAMES)}

C = Parameters(**name2idx)

del name2idx
