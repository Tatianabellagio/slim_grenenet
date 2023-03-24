import requests

## defining download
## ecotypes to download (coma separated)
## example 9998,9999

with open('ecotypes_grenenet_1001g.txt') as f:
    lines = f.readlines()

ecotypes = str(lines).replace('[\'', '').replace('\']', '').replace(' ', '')


## regions of genome to donwload 
## example Chr1:1000..1010,Chr2:2000..2010

regions = 'Chr1:0..1000'

regions = 'Chr1:1000..552138'
#ecotypes = '9998,9999'

## example 
# https://tools.1001genomes.org/api/v1/vcfsubset/strains/9998,9999/regions/Chr1:1000..1010,Chr2:2000..2010/type/fullgenome/format/vcf

url = f'https://tools.1001genomes.org/api/v1/vcfsubset/strains/{ecotypes}/regions/{regions}/type/fullgenome/format/vcf'

req = requests.get(url)
open('test.vcf', 'wb').write(req.content)