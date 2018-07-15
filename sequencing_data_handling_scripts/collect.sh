#bash scripts/collect.sh REMIX13.UP01 np

base=$1
gen=$2

paste -d"\t" \
<(cut -f 1,2,3,4,5 depths/UP/$base.I1.$gen.pilon.depth.bed) \
<(cut -f5 depths/UP/$base.I2.$gen.pilon.depth.bed) \
<(cut -f5 depths/UP/$base.I3.$gen.pilon.depth.bed) \
<(cut -f5 depths/UP/$base.I4.$gen.pilon.depth.bed) \
> depths/all/$base.I.$gen.depth.bed 

