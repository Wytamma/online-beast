from pathlib import Path
from Bio.Align import MultipleSeqAlignment
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class BeastXML:
    """Class for editing BEAST XML files."""

    file_name: Path
    traits: list
    xml: ElementTree  # what is?

    def __init__(
        self,
        file_name: Path,
        traits: list = [],
    ):
        self.file_name = file_name
        self.traits = traits
        self.xml = self._load_xml()
        # should check to see if there are traits in the xml
        # that are not in self.traits.

    def _load_xml(self):
        return ET.parse(self.file_name)

    @property
    def alignment(self) -> MultipleSeqAlignment:
        msa = MultipleSeqAlignment([])
        data = self.xml.find("data")
        for sequence_el in data:
            msa.append(
                SeqRecord(
                    Seq(sequence_el.get("value")),
                    id=sequence_el.get("taxon"),
                    description="",
                )
            )
        return msa

    def get_sequence_ids(self) -> list:
        return [s.id for s in self.alignment]

    def _add_trait_data(self, sequence_id, traitname, deliminator, group):
        trait_el = self.xml.find(f".//*[@traitname='{traitname}']")
        if not trait_el:
            raise ValueError(f"Could not find trait with traitname '{traitname}'")

        trait = sequence_id.split(deliminator)[group]
        trait_el.set("value", f"{trait_el.get('value')},{sequence_id}={trait}")

    def add_sequence(self, record: Seq):
        for trait in self.traits:
            self._add_trait_data(record.id, **trait)
        data = self.xml.find("data")
        sequence_el = ET.Element(
            "sequence",
            {
                "id": f"seq_{record.id}",
                "taxon": record.id,
                "totalcount": "4",
                "value": record.seq,
            },
        )
        data.append(sequence_el)

    def write(self, out_file=None) -> None:
        if not out_file:
            out_file = self.file_name
        self.xml.write(out_file)
