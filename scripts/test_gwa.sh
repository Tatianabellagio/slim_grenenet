base_directory=$(pwd)
echo $base_directory qE '^[^0]|0[^.].*$'
if grep -qE '^[^0]|0[^.0]|0\.[^0]|0\.0*[1-9]' "data/pheno.txt"; then
    echo 'hola'
else echo 'chau'

fi


    