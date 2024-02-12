from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "kf1",
    "kr1",
    "V2",
    "K2",
    "V7",
    "K7",
    "V8",
    "K8",
    "V9",
    "K9",
    "n9",
    "V10",
    "K10",
    "n10",
    "kf11",
    "kf12",
    "kr12",
    "kf13",
    "kf14",
    "kf15",
    "kr15",
    "V16",
    "K16",
    "V17",
    "K17",
    "kf18",
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
