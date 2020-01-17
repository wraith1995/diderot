cd ..
make local-install
cd test
../bin/diderotc --log --dump-pt --dump-ast --dump-simple --dump-high --dump-mid --dump-low --dump-tree fun.diderot
#../bin/diderotc --exec --log --dump-pt --dump-ast --dump-simple --dump-high --dump-mid --dump-low --dump-tree justTypes.diderot
