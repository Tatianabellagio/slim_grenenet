# Paths to the input text files
tree_file="data/og_tree_offset.trees"
# Path to the SLiM script
slim_script="scripts/tree_after_cloning.slim"
## output
output="data/og_tree_offset_after_cloning.trees"
# Construct the SLiM command
slim_command="slim -d \"tree='$tree_file'\" -d \"output='$output'\" \"$slim_script\""
    
# Echo the command
echo "Running: $slim_command"

# Uncomment the following line to execute the SLiM command
eval "$slim_command"

echo "made it here?"