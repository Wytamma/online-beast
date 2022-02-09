from pathlib import Path
from Bio.Align import MultipleSeqAlignment
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from datetime import datetime


class BeastXML:
    """Class for editing BEAST XML files."""

    file_name: Path
    traits: list
    xml: ElementTree  # what is?
    date_trait: bool
    date_format: str
    date_deliminator: str

    def __init__(
        self,
        file_name: Path,
        traits: list = [],
        date_trait: bool = True,
        date_format: str = "%Y-%m-%d",
        date_deliminator: str = "_",
    ):
        self.file_name = file_name
        self.traits = traits
        self.xml = self._load_xml()
        # should check to see if there are traits in the xml
        # that are not in self.traits.
        if date_trait:
            try:
                self._get_first_element_by_attribute("traitname", "date")
            except ValueError:
                raise ValueError(
                    f"Could not find date trait (use --no-date-trait flag?)."
                )
        self.date_trait = date_trait
        self.date_format = date_format
        self.date_deliminator = date_deliminator

    def _load_xml(self):
        return ET.parse(self.file_name)

    def _get_first_element_by_attribute(self, attr, value):
        el = self.xml.find(f".//*[@{attr}='{value}']")
        if not el:
            raise ValueError(f"Could not find trait with {attr}='{value}'")
        return el

    def _add_trait(self, sequence_id, traitname, deliminator, group):
        trait_el = self._get_first_element_by_attribute("traitname", traitname)
        trait = sequence_id.split(deliminator)[group]
        trait_el.set("value", f"{trait_el.get('value')},{sequence_id}={trait}")

    def _add_date_trait(self, sequence_id):
        date_el = self._get_first_element_by_attribute("traitname", "date")
        date = None
        for potential_date in sequence_id.split(self.date_deliminator):
            try:
                date = datetime.strptime(potential_date, self.date_format).strftime(
                    self.date_format
                )
            except ValueError:
                pass
        if not date:
            raise ValueError(
                f"Could not parse date trait with format '{self.date_format}' and deliminator '{self.date_deliminator}'"
            )
        date_el.set("value", f"{date_el.get('value')},{sequence_id}={date}")

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

    def add_sequence(self, record: Seq):
        if self.date_trait:
            self._add_date_trait(record.id)
        for trait in self.traits:
            self._add_trait(record.id, **trait)
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
