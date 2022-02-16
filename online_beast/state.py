from pathlib import Path
from Bio import Phylo
from io import StringIO
from Bio.Phylo.BaseTree import TreeMixin


class StateTree:
    """Class for editing BEAST state files."""

    file_name: str
    tree: TreeMixin

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
        self,
        branch_id,
        index: int,
        sampling_time_delta: float = 0,
        graft_point: float = 0.5,
    ):
        try:
            clade = next(
                c for c in self.tree.get_terminals() if c.name == str(branch_id)
            )
        except StopIteration:
            raise ValueError("Could not find branch with id", branch_id)
        total_branch_length = clade.branch_length

        internal = False

        starting_total_branch_length = self.tree.total_branch_length()

        while total_branch_length + sampling_time_delta < 0:
            try:
                clade = self.get_parent(clade)
            except IndexError:
                # root
                ValueError("Sample date before root...")
            internal = True
            total_branch_length += clade.branch_length

        if sampling_time_delta < 0:
            parent_branch_length = (
                total_branch_length - abs(sampling_time_delta)
            ) * graft_point
            # parent_branch_length = (clade.branch_length) * graft_point
            og_parent = clade.branch_length
            clade.branch_length = parent_branch_length
            branch_length = total_branch_length - parent_branch_length
        else:
            parent_branch_length = (total_branch_length) * graft_point
            clade.branch_length = parent_branch_length
            branch_length = total_branch_length - parent_branch_length

        for leaf in self.tree.get_terminals():
            if int(leaf.name) < index:
                continue
            leaf.name = str(int(leaf.name) + 1)

        if internal:
            descendants = 1
        else:
            descendants = 2

        clade.split(n=descendants, branch_length=branch_length)
        clade.name = None

        if internal:
            children = list(clade.clades)
            clade.clades = []
            clade.split(branch_length=0)
            clade.clades[1] = children[-1]
            clade.clades[0].branch_length = og_parent - clade.branch_length
            clade.clades[0].clades = children[:-1]
            clade.clades[0].confidence = None
            clade.clades[0].name = None
        else:
            original_branch = clade.clades[0]
            if branch_id >= index:
                branch_id += 1
            original_branch.name = str(branch_id)

        new_branch = clade.clades[1]
        if sampling_time_delta:
            if new_branch.branch_length + sampling_time_delta < 0:
                raise ValueError("negative branch bad :(")
            new_branch.branch_length = new_branch.branch_length + sampling_time_delta

        error = 0.00001  # for floats
        if (
            abs(
                (self.tree.total_branch_length() - new_branch.branch_length)
                - starting_total_branch_length
            )
            > error
        ):
            print(
                self.tree.total_branch_length(),
                new_branch.branch_length,
                (self.tree.total_branch_length() - new_branch.branch_length),
                starting_total_branch_length,
            )
            raise ValueError("Tree length is wrong!")

        if index != None:
            new_branch.name = str(index)
        else:
            new_branch.name = str(len(self.tree.get_terminals()) - 1)  # zero indexed
        return new_branch

    def draw(self):
        Phylo.draw_ascii(self.tree)

    def write(self, out_file=None) -> None:
        with open(self.file_name) as f:
            state_file_lines = f.readlines()

        for c in self.tree.find_clades():
            if c.confidence:
                c.confidence = None

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
