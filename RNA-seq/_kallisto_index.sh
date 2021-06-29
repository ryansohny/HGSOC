#!/bin/bash

kallisto=<kallisto exec> # ver 0.46.1

$kallisto index \
-i GRCh37_rna.idx \
GRCh37_latest_rna.fna.gz # ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_rna.fna.gz
