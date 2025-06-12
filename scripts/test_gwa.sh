
# Count the number of columns in the file
num_columns=$(awk -F ' ' 'NR==1 {print NF; exit}' "../data/greneNet_final_v1.1.recode.fam")

# Check if the number of columns is not 6
if [ "$num_columns" -ne 6 ]; then
    echo "fam file wrong dimensions"
    exit 1
else
    # Check if all values in the 6th column are 0
    if awk -F ' ' '{print $6}' "../data/greneNet_final_v1.1.recode.fam" | grep -qE '^[^0]|0[^.].*$'; then
        echo "run_gemma"
    else
        echo "dont run gemma"
    fi
fi


    