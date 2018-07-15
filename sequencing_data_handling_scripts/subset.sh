#!/bin/bash

#val=$( echo " $3 - $2 " | bc )

gunzip -cd $1 | sed -n $2,$3p |gzip > ~/REMIX/temp/$4.contaminant.fq.gz

