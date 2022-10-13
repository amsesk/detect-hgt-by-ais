#!/bin/bash

GFF=$1
cat $GFF | cut -f1,9 | sed -r 's/(ID[=][0-9_]*).*/\1/' | grep -v "^#" > ${GFF}.idmap
