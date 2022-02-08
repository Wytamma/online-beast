from io import StringIO
from pathlib import Path
from Bio import Phylo, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import xml.etree.ElementTree as ET
import typer
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from datetime import datetime

app = typer.Typer()


def get_MSA_from_xml(xml_file: Path) -> MultipleSeqAlignment:
    align = MultipleSeqAlignment([])
    data = ET.parse(xml_file).find("data")
    for sequence_el in data:
        align.append(
            SeqRecord(
                Seq(sequence_el.get("value")),
                id=sequence_el.get("taxon"),
                description="",
            )
        )
    return align


def get_tree_from_state_file(state_file: Path):
    with open(state_file) as f:
        lines = f.readlines()

    tree = lines[1]

    return tree.split(">")[1].split("</")[0]


def find_closest_sequence(MSA: MultipleSeqAlignment, new_sequence):
    aligner = PairwiseAligner()
    max_score = None
    seq_id = None
    with typer.progressbar(MSA) as progress:
        for i, sequence in enumerate(progress):
            score = sum(
                xi != yi for xi, yi in zip(str(sequence.seq), str(new_sequence))
            )
            if not max_score:
                max_score = score
                seq_id = i
            elif score < max_score:
                max_score = score
                seq_id = i
    if seq_id == "None":
        raise Exception("No Seq found?")
    return seq_id, max_score


def add_node_to_tree(tree, nearest_seq_id, name=None):
    clade = next(c for c in tree.get_terminals() if c.name == str(nearest_seq_id))
    clade.branch_length = clade.branch_length / 2
    clade.split(branch_length=clade.branch_length)
    clade.confidence = 1
    clade.clades[0].name = str(nearest_seq_id)
    if name:
        clade.clades[1].name = name
    clade.name = None
    return clade.clades[1]


def add_new_tree_to_state_file(tree, state_file, output):
    writer = Phylo.NewickIO.Writer([tree])

    newick_tree = next(writer.to_strings(format_branch_length="%1.17f"))

    with open(state_file) as f:
        state_file_lines = f.readlines()

    old_tree_line = state_file_lines[1]

    opening_tag = old_tree_line.split(">")[0]
    closting_tag = old_tree_line.split("</")[-1]

    state_file_lines[1] = f"{opening_tag}>{newick_tree}</{closting_tag}"
    if output:
        state_file = Path(f"{output}.state")
    with open(state_file, "w") as f:
        f.writelines(state_file_lines)
    return state_file


def add_new_sequence_to_xml(
    xml_file, sequence, seq_id, xml_output, dateFormat, deliminator
):
    xml_tree = ET.parse(xml_file)
    data = xml_tree.find("data")
    sequence_el = ET.Element(
        "sequence",
        {"id": f"seq_{seq_id}", "taxon": seq_id, "totalcount": "4", "value": sequence},
    )
    data.append(sequence_el)
    trait = xml_tree.find(".//*[@traitname='date']")
    if trait:
        date = None
        for potential_date_trait in seq_id.split(deliminator):
            try:
                date = datetime.strptime(potential_date_trait, dateFormat).strftime(
                    dateFormat
                )
            except:
                pass
        if date:
            typer.echo(f"Adding date data: {date}")
            trait.set("value", f"{trait.get('value')},{seq_id}={date}")
    if xml_output:
        xml_file = xml_output
    xml_tree.write(xml_file)
    return xml_file


def get_sequences_to_add(MSA, new_seq_MSA_fasta):
    records = SeqIO.parse(new_seq_MSA_fasta, "fasta")
    MSA_Ids = [s.id for s in MSA]
    return [record for record in records if record.id not in MSA_Ids]


@app.command()
def main(
    xml_file: Path,
    fasta_file: Path,
    state_file: Path = None,
    output: Path = None,
    dateFormat: str = "%d/%m/%Y",
    deliminator: str = "_",
):
    MSA = get_MSA_from_xml(xml_file)
    sequences_to_add = get_sequences_to_add(MSA, fasta_file)
    for sequence in sequences_to_add:
        typer.echo(f"Adding new sequence: {sequence.id}")
        if len(sequence) != MSA.get_alignment_length():
            raise ValueError("Sequences must all be the same length")
        closest_seq_id, max_score = find_closest_sequence(MSA, sequence.seq)
        if not state_file:
            state_file = Path(f"{xml_file}.state")
        newick_tree = get_tree_from_state_file(state_file)
        tree = Phylo.read(StringIO(newick_tree), "newick")
        new_clade = add_node_to_tree(tree, closest_seq_id, name=sequence.id)
        Phylo.draw_ascii(tree)
        new_clade.name = str(len(tree.get_terminals()) - 1)
        state_file = add_new_tree_to_state_file(tree, state_file, output=output)
        xml_file = add_new_sequence_to_xml(
            xml_file, sequence.seq, sequence.id, output, dateFormat, deliminator
        )
