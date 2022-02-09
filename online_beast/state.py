from pathlib import Path
from Bio import Phylo
from io import StringIO


class StateTree:
    """Class for editing BEAST state files."""

    file_name: str
    tree: str

    def __init__(self, file_name: Path):
        self.file_name = file_name
        self.tree = self._get_tree()

    def _get_tree(self):
        with open(self.file_name) as f:
            lines = f.readlines()
        tree_line = lines[1]
        newick_tree = tree_line.split(">")[1].split("</")[0]
        return Phylo.read(StringIO(newick_tree), "newick")

    @property
    def newick(self):
        writer = Phylo.NewickIO.Writer([self.tree])
        return next(writer.to_strings(format_branch_length="%1.17f"))

    def graft(self, branch_id, graft_point: float = 0.5, branch_length: float = None):
        clade = next(c for c in self.tree.get_terminals() if c.name == str(branch_id))
        if not clade:
            raise ValueError("Could not find branch with id", branch_id)
        total_branch_length = clade.branch_length
        clade.branch_length = total_branch_length * graft_point
        clade.split(branch_length=total_branch_length * (1 - graft_point))
        clade.confidence = 1
        clade.name = None
        original_branch, new_branch = clade.clades
        original_branch.name = str(branch_id)
        new_branch.name = str(len(self.tree.get_terminals()) - 1)
        if branch_length:
            new_branch.branch_length = branch_length
        return new_branch

    def draw(self):
        Phylo.draw_ascii(self.tree)

    def write(self, out_file=None) -> None:
        with open(self.file_name) as f:
            state_file_lines = f.readlines()

        old_tree_line = state_file_lines[1]
        opening_tag = old_tree_line.split(">")[0]
        closting_tag = old_tree_line.split("</")[-1]
        state_file_lines[1] = f"{opening_tag}>{self.newick}</{closting_tag}"

        if not out_file:
            out_file = self.file_name
        else:
            out_file = Path(f"{out_file}.state")

        with open(out_file, "w") as f:
            f.writelines(state_file_lines)
