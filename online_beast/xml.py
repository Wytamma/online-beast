from pathlib import Path
from typing import List
from Bio.Align import MultipleSeqAlignment

from lxml import etree as ET
from xml.etree.ElementTree import ElementTree
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from datetime import datetime


class BeastXML:
    """Class for editing BEAST XML files."""

    file_name: Path
    traits: List[dict]
    xml: ElementTree
    date_trait: bool
    date_format: str
    date_delimiter: str

    def __init__(
        self,
        file_name: Path,
        traits: List[dict] = [],
        date_trait: bool = True,
        date_format: str = "%Y-%m-%d",
        date_delimiter: str = "_",
    ):
        self.file_name = file_name
        self.xml = self._load_xml()
        self.traits = traits
        for trait in traits:
            try:
                self._get_first_element_by_attribute("traitname", trait)
            except ValueError:
                raise ValueError(f"Could not find '{trait}' trait in xml.")
        # should check to see if there are traits in the xml
        # that are not in self.traits?
        if date_trait:
            try:
                self._get_first_element_by_attribute("traitname", "date")
            except ValueError:
                raise ValueError(f"Could not find 'date' trait in xml..")
        self.date_trait = date_trait
        self.date_format = date_format
        self.date_delimiter = date_delimiter

    def _load_xml(self):
        return ET.parse(str(self.file_name))

    def _get_first_element_by_attribute(self, attr, value):
        el = self.xml.find(f".//*[@{attr}='{value}']")
        if el is None:
            raise ValueError(f"Could not find trait with {attr}='{value}'")
        return el

    def _add_trait(self, sequence_id, traitname):
        trait_el = self._get_first_element_by_attribute("traitname", traitname)
        trait = self.get_trait_data(sequence_id, traitname)
        trait_el.set("value", f"{trait_el.get('value')},{sequence_id}={trait}")

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
        msa.sort()
        return msa

    def get_sequence_ids(self) -> list:
        return [s.id for s in self.alignment]

    def get_trait_data(self, sequence_id, traitname):
        if traitname == "date" and self.date_trait:
            date = None
            for potential_date in sequence_id.split(self.date_delimiter):
                try:
                    date = datetime.strptime(potential_date, self.date_format).strftime(
                        self.date_format
                    )
                except ValueError:
                    pass
            if not date:
                raise ValueError(
                    f"Could not parse date trait with date format '{self.date_format}' and date delimiter '{self.date_delimiter}'"
                )
            return date

        trait = next(t for t in self.traits if t["traitname"] == traitname)

        return sequence_id.split(trait["delimiter"])[trait["group"]]

    def add_sequence(self, record: SeqRecord):
        if record.id in self.get_sequence_ids():
            raise ValueError(f"New sequence id must be unique ({record.id})")
        try:
            self.alignment.append(record)
        except ValueError as e:
            raise e
        if self.date_trait:
            self._add_trait(record.id, "date")
        for trait in self.traits:
            self._add_trait(record.id, trait["traitname"])
        data = self.xml.find("data")
        sequence_el = ET.Element(
            "sequence",
            {
                "id": f"seq_{record.id}",
                "spec": "Sequence",
                "taxon": record.id,
                "totalcount": "4",
                "value": str(record.seq),
            },
        )
        sequence_el.tail = "\n"
        data.append(sequence_el)
        data[:] = sorted(data, key=lambda x: x.get("taxon"))

    def write(self, out_file=None) -> None:
        if not out_file:
            out_file = self.file_name
        ET.indent(self.xml, space="    ", level=0)
        with open(out_file, "w") as f:
            xml_string = ET.tostring(self.xml, pretty_print=True, encoding=str)
            f.write(xml_string)
