#!/bin/bash


bash scripts/alignrest.sh REMIX13.I pilonrefs/np.pilon.fna REMIX re_bams   
bash scripts/alignrest.sh REMIX14.I pilonrefs/gg.pilon.fna REMIX re_bams
bash scripts/alignrest.sh REMIX15.I pilonrefs/np.pilon.fna REMIX re_bams
bash scripts/alignrest.sh REMIX16.I pilonrefs/gg.pilon.fna REMIX re_bams

bash scripts/alignrest.sh REMIX13.U pilonrefs/np.pilon.fna REMIX re_bams   
bash scripts/alignrest.sh REMIX14.U pilonrefs/gg.pilon.fna REMIX re_bams
bash scripts/alignrest.sh REMIX15.U pilonrefs/np.pilon.fna REMIX re_bams
bash scripts/alignrest.sh REMIX16.U pilonrefs/gg.pilon.fna REMIX re_bams
