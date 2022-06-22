import re
import textwrap
from dataclasses import dataclass


@dataclass
class SeqRecord:
    id: str  # noqa
    name: str
    seq: str
    haplotype: int = -1
    node: int = -1
    population: str = None

    def __post_init__(self):
        # Make sure no line breaks or spaces in sequence
        self.seq = re.sub(r"[^a-zA-Z]", "", self.seq)
        print(type(self.seq))

    @property
    def seqid(self):
        return f"{self.population}-{self.id}-{self.haplotype}"

    @property
    def description(self):
        return (
            f"id:{self.id}, node:{self.node}, name:{self.name}, "
            f"haplotype:{self.haplotype}, "
            f"population:{self.population}"
        )

    def as_array(self):
        return list(self.seq)

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        s = f">{self.id} {self.description}\n"
        s += "\n".join(
            textwrap.wrap(
                self.seq,
            )
        )
        return s

    def __repr__(self):
        return str(self)
