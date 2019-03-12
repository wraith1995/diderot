rm *.nrrd
cd seg
make
cd ..
cd sub
make
cd ..
cd find
make
cd ..
python control.py
python render.py
unu quantize -b 16 -i out_0.nrrd -o out.png
xdg-open out.png
