unu slice -a 0 -p 0 -i $1 -o x.nrrd
unu slice -a 0 -p 1 -i $1 -o y.nrrd
unu slice -a 0 -p 2 -i $1 -o z.nrrd
unu join -i x.nrrd y.nrrd -a 0 -incr | unu sselect -i - -s z.nrrd -a 1 -th 0.5 -o - x | unu jhisto -b 400 400 | unu quantize -b 8 -max 1% -o $2
xdg-open $2
