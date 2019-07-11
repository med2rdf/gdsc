# gdsc
GDSC RDF converter

## Get Dataset  
Download data from the following web site .
* [GDSC](https://www.cancerrxgene.org/downloads)　["PANCANCER_IC"](https://www.cancerrxgene.org/translation/drug/download#ic50) ["PANCANCER_ANOVA"](https://www.cancerrxgene.org/translation/drug/download#anova)

```
- PANCANCER_IC_Wed Oct 10 03_15_39 2018.csv
- PANCANCER_ANOVA_Wed Oct 10 03_46_43 2018.csv
```

Correct wrong data about the format of csv.
For example, "Other, kinase" -> "Other kinase".

* [Cellosaurus](https://web.expasy.org/cellosaurus/)
```
- cellosaurus.txt
```

ftp://ftp.expasy.org/databases/cellosaurus
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

## Create RDF files
Check script and raw files in container .
```
/work/
    ┣ data/
    ┃   ┗ raw/
    ┃       ┣ PANCANCER_IC_Wed Oct 10 03_15_39 2018.csv
    ┃       ┣ PANCANCER_ANOVA_Wed Oct 10 03_46_43 2018.csv
    ┃       ┗ cellosaurus.txt
    ┃
    ┣ scripts/
    ┃   ┣ __init__.py
    :   ┣ gdsc_to_rdf.py
        ┣ omics_to_rdf.py           
        ┗ preprocessing.py
```

Run as follow. (It takes about 5 minites / 10 core)
```
cd /work/scripts
python3 preprocessing.py
python3 gdsc_to_rdf.py
```

Reference owl are here.

- https://github.com/med2rdf/med2rdf-ontology/blob/master/med2rdf-0.06.ttl
- http://www.bioassayontology.org/bao/bao_complete_merged.owl
- http://semanticscience.org/ontology/sio.owl

