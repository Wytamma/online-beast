import filecmp
from lxml import etree as ET
import pytest
from online_beast import BeastXML
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


ebola_file_name = "data/ebola.xml"
GTR_file_name = "data/testGTR.xml"
bx_ebola = BeastXML(ebola_file_name, date_format="%d/%m/%Y")
bx_GTR = BeastXML(GTR_file_name, date_trait=False)


def test_init():
    assert bx_ebola.file_name == ebola_file_name
    assert ET.tostring(bx_ebola.xml) == ET.tostring(ET.parse(str(ebola_file_name)))
    assert bx_ebola.traits == []
    assert bx_ebola.date_trait == True
    assert bx_ebola.date_format == "%d/%m/%Y"
    assert bx_ebola.date_delimiter == "_"


def test_no_date_trait():
    assert bx_GTR.file_name == GTR_file_name
    assert bx_GTR.date_trait == False
    assert bx_GTR.date_format == "%Y-%m-%d"
    assert bx_GTR.date_delimiter == "_"


def test_get_trait():
    el = bx_ebola._get_first_element_by_attribute("traitname", "date")
    assert el.get("traitname") == "date"


def test_get_alignment():
    assert type(bx_ebola.alignment) == MultipleSeqAlignment
    assert type(bx_GTR.alignment) == MultipleSeqAlignment
    assert bx_ebola.alignment[0].id == "KM034550_25/05/2014"


def test_get_sequence_ids():
    assert len(bx_ebola.get_sequence_ids()) == 13


def test_get_trait_data():
    assert (
        bx_ebola.get_trait_data(sequence_id="KM034550_25/05/2014", traitname="date")
        == "25/05/2014"
    )


def test_add_sequence():
    id = "KM034550a_25/05/2014"
    sequence = bx_ebola.alignment[0].seq

    with pytest.raises(ValueError):
        bx_ebola.add_sequence(SeqRecord(Seq("ATG"), id=id))

    with pytest.raises(ValueError):
        bx_ebola.add_sequence(SeqRecord(sequence, id="KM034550_25/05/2014"))

    bx_ebola.add_sequence(SeqRecord(sequence, id=id))
    assert bx_ebola.alignment[-1].id == id
    assert bx_ebola.alignment[-1].seq == sequence


def test_add_sequence():
    bx_ebola.write("ebola.xml")
    assert filecmp.cmp("data/online-result-ebola.xml", "ebola.xml") == True
