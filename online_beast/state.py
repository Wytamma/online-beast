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

    def get_parent(self, child_clade):
        node_path = self.tree.get_path(child_clade)
        return node_path[-2]

    def graft(
        self, branch_id, sampling_time_delta: float = 0, graft_point: float = 0.5
    ):
        clade = next(c for c in self.tree.get_terminals() if c.name == str(branch_id))
        if not clade:
            raise ValueError("Could not find branch with id", branch_id)
        total_branch_length = clade.branch_length

        running_total = total_branch_length
        internal = False
        if total_branch_length + sampling_time_delta < 0:
            while total_branch_length + sampling_time_delta < 0:
                try:
                    clade = self.get_parent(clade)
                except IndexError:
                    # root
                    break
                internal = True
                total_branch_length += clade.branch_length

            # raise ValueError("negative branch bad :( - add less samples?")
            # need to write while loop to find parent node that works

        if sampling_time_delta <= 0:
            parent_branch_length = (
                total_branch_length - abs(sampling_time_delta)
            ) * graft_point
            clade.branch_length = parent_branch_length
            branch_length = total_branch_length - parent_branch_length
        else:
            parent_branch_length = (total_branch_length) * graft_point
            clade.branch_length = parent_branch_length
            branch_length = total_branch_length - parent_branch_length

        if internal:
            descendants = 1
        else:
            descendants = 2
        clade.split(n=descendants, branch_length=branch_length)
        clade.name = None

        if not internal:
            original_branch = clade.clades[0]
            original_branch.name = str(branch_id)

        new_branch = clade.clades[-1]
        new_branch.name = str(len(self.tree.get_terminals()) - 1)  # zero indexed
        if sampling_time_delta:
            if new_branch.branch_length + sampling_time_delta < 0:
                new_branch.name = "new_branch"
                self.draw()
                print(self.tree)
                raise ValueError("negative branch bad :(")
            new_branch.branch_length = new_branch.branch_length + sampling_time_delta
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
