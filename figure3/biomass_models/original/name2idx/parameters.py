from dataclasses import make_dataclass
from typing import Dict, List

NAMES: List[str] = [
    "kf1",
    "kr1",
    "kf2",
    "kr2",
    "kf3",
    "kr3",
    "V4",
    "K4",
    "kf5",
    "kr5",
    "kf6",
    "kr6",
    "kf7",
    "kr7",
    "V8",
    "K8",
    "kf9",
    "kr9",
    "kf10",
    "kr10",
    "kf11",
    "kr11",
    "kf12",
    "kr12",
    "kf13",
    "kr13",
    "kf14",
    "kr14",
    "kf15",
    "kr15",
    "V16",
    "K16",
    "kf17",
    "kr17",
    "kf18",
    "kr18",
    "kf19",
    "kr19",
    "kf20",
    "kr20",
    "kf21",
    "kr21",
    "kf22",
    "kr22",
    "kf23",
    "kr23",
    "kf24",
    "kr24",
    "kf25",
    "kr25",
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
