#!/bin/sh
#
# wrapper script for MkFrags.mkMakefile function
#

PROG=mkmk

if [ $# != 1 ] ; then
  echo "usage: $PROG.sh <dir>"
  exit 1
fi

DIR=$1

SRC=/home/teocollin/gitcode/diderot/src/compiler/gen/fragments/sources.cm

if test "smlnj" = "smlnj" ; then
/home/teocollin/installs/smlnj/bin/sml $SRC 2> /dev/null 1>/dev/null <<XXXX
MkFrags.mkMakefile "$DIR";
XXXX
exit $?
elif test "smlnj" = "mlton" ; then
  HERE=$(pwd)
  cd /home/teocollin/gitcode/diderot/src/compiler/gen/fragments
  make -s $PROG || exit $1
  cd $HERE
  /home/teocollin/gitcode/diderot/src/compiler/gen/fragments/$PROG $DIR || exit $1
  exit 0
else
  echo "unknown SML implementation"
  exit 1
fi
