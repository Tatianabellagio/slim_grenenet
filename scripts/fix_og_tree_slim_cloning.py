import pyslim
import tskit

baseline_tree = tskit.load('og_tree_offset.trees')

tables = baseline_tree.dump_tables()

tables.sites.clear()
for s in baseline_tree.sites():
    tables.sites.append(s.replace(ancestral_state=""))

tables.mutations.clear()
mm = pyslim.default_slim_metadata('mutation_list_entry')
for k, m in enumerate(baseline_tree.mutations()):
    tables.mutations.append(
        m.replace(derived_state=str(k), metadata={'mutation_list': [mm]}))

og_tree_offset_for_slim = tables.tree_sequence()

og_tree_offset_for_slim.dump('og_tree_offset_for_slim.trees')
