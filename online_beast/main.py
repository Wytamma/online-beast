from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import typer
from typing import List, Optional, Tuple
from datetime import datetime

from .xml import BeastXML
from .state import StateTree

app = typer.Typer()


def find_closest_sequence(
    sequences_from_xml: MultipleSeqAlignment, new_sequence_record: SeqRecord
):
    max_score = None
    seq_id = None
    with typer.progressbar(sequences_from_xml) as progress:
        for i, sequence in enumerate(progress):

            if sequence.id == new_sequence_record.id:
                continue
            score = sum(
                xi != yi
                for xi, yi in zip(str(sequence.seq), str(new_sequence_record.seq))
            )
            if max_score == None:
                max_score = score
                seq_id = i
            elif score < max_score:
                max_score = score
                seq_id = i
    if seq_id == None:
        raise Exception("No Seq found?")
    return seq_id, max_score


def get_sequences_to_add(fasta_file, list_of_seq_ids: list):
    records = SeqIO.parse(fasta_file, "fasta")
    return [record for record in records if record.id not in list_of_seq_ids]


def decimal_year(dt: datetime):
    year_part = dt - datetime(year=dt.year, month=1, day=1)
    year_length = datetime(year=dt.year + 1, month=1, day=1) - datetime(
        year=dt.year, month=1, day=1
    )
    return dt.year + year_part / year_length


@app.command()
def main(
    xml_file: Path,
    fasta_file: Path,
    state_file: Path = None,
    output: Path = None,
    template: bool = typer.Option(
        False,
        help="XML file is a beast template. New sequences must be append to the end of the original fasta file.",
    ),
    date_trait: bool = True,
    date_format: str = "%Y-%m-%d",
    date_delimiter: str = "_",
    trait: Optional[List[str]] = typer.Option(
        None,
        help="Trait information 'traitname delimiter group' string seperated by spaces",
    ),
):
    if not state_file:
        state_file = Path(f"{xml_file}.state")
    state_tree = StateTree(state_file)
    traits = [
        {
            "traitname": t.split(" ")[0],
            "delimiter": t.split(" ")[1],
            "group": int(t.split(" ")[2]),
        }
        for t in trait
    ]

    beast_xml = BeastXML(
        xml_file,
        traits,
        date_trait=date_trait,
        date_format=date_format,
        date_delimiter=date_delimiter,
    )
    if template:
        date_trait = False
        records = list(SeqIO.parse(fasta_file, "fasta"))
        number_of_tips = state_tree.tree.count_terminals()
        alignment = MultipleSeqAlignment(records[:-number_of_tips])
        sequences_to_add = records[number_of_tips:]
    else:
        sequences_to_add = get_sequences_to_add(
            fasta_file, beast_xml.get_sequence_ids()
        )

    if not sequences_to_add:
        typer.echo("No new sequences found in the fasta file.")
        raise typer.Exit(code=1)

    typer.echo(f"Adding {len(sequences_to_add)} new sequences")
    for sequence in sequences_to_add:
        typer.echo(f"Adding new sequence: {sequence.id}")
        if not template:
            alignment = beast_xml.alignment
        if len(sequence) != alignment.get_alignment_length():
            raise ValueError("Sequences must all be the same length")

        closest_tree_node_id, max_score = find_closest_sequence(alignment, sequence)
        sampling_time_delta = 0
        if date_trait:
            closest_tree_node_date = beast_xml.get_trait_data(
                alignment[closest_tree_node_id].id, traitname="date"
            )
            closest_tree_node_dt = datetime.strptime(
                closest_tree_node_date, beast_xml.date_format
            )
            new_tree_node_date = beast_xml.get_trait_data(sequence.id, traitname="date")
            new_tree_node_dt = datetime.strptime(
                new_tree_node_date, beast_xml.date_format
            )
            sampling_time_delta = decimal_year(new_tree_node_dt) - decimal_year(
                closest_tree_node_dt
            )

        ids = [s.id for s in alignment]
        ids.append(sequence.id)
        index = sorted(ids).index(sequence.id)

        new_clade = state_tree.graft(
            closest_tree_node_id, index=index, sampling_time_delta=sampling_time_delta
        )
        # name = new_clade.name
        # new_clade.name = f"{sequence.id}"
        # state_tree.draw()
        # new_clade.name = name
        if not template:
            beast_xml.add_sequence(sequence)
        else:
            alignment.append(sequence)

    if not template:
        beast_xml.write(out_file=output)
    # else:
    #     for c in state_tree.tree.get_terminals():
    #         c.name = str(int(c.name) + 1)
    state_tree.draw()
    state_tree.write(out_file=output)
