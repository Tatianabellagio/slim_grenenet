# Paths to the input text files
tree_file="data/og_tree_offset.trees"
## output
output_file="data/og_tree_offset_after_cloning.trees"

echo $tree_file
echo $output_file

slim \
    -d "tree='$tree_file'" \
    -d "output='$output_file'" \
    scripts/tree_after_cloning.slim


