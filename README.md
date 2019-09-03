# gdsc
GDSC RDF converter

## Check docker-compose.yaml
Edit docker-compose.yaml and build docker container using docker-compose .
```
docker-compose.yaml
--------------------------------
#      - [directory downloading gdsc data]:/work/data/raw         for mount
#      - [RDF file output directory]:/work/output                 for mount
--------------------------------
```

## Create docker container
```
docker-compose up -d

docker-compose exec gdsc_cnvtr bash
```

Check script files in container .
```
/work/
    ┣ scripts/
    ┃   ┣ __init__.py
    ┃   ┣ download_files.py
    ┃   ┣ gdsc_to_rdf.py
    ┃   ┣ omics_to_rdf.py
    ┃   ┗ preprocessing.py
    :
```


## Get Dataset
Run as follows. Select GDSC version (1 or 2).
```
$ mkdir output
$ cs scripts/
$ python3 download_files.py -v 2
```

Check raw files in container .
```
/work/
    ┣ data/
    ┃   ┗ raw/
    ┃       ┣ PANCANCER_IC.csv
    ┃       ┣ PANCANCER_ANOVA.csv
    ┃       ┗ cellosaurus.txt
    :
```


## Remove lines breaking rules
Sometimes the format of PANCANCER_ANOVA.csv is corrupted. As see bellow, the number of columns exceed that of the header if it contains multiple drug name.
```
drug_name,drug_id,drug_target,target_pathway,feature_name, ...
Navitoclax,1011,BCL2, BCL-XL, BCL-W,Apoptosis regulation,ABCB1_mut,13,736, ...
```

You can remove the lines by following command.
```
$ python3 remove_lines.py -v 2
```

## Create RDF files
Run as follow. (It takes about 5 minites / 10 core)
```
$ python3 preprocessing.py -v 2
$ python3 gdsc_to_rdf.py -v 2
```

Reference owl are here.

- https://github.com/med2rdf/med2rdf-ontology/blob/master/med2rdf-0.06.ttl
- http://www.bioassayontology.org/bao/bao_complete_merged.owl
- http://semanticscience.org/ontology/sio.owl

