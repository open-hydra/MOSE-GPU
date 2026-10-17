#!/bin/bash -
#===============================================================================
#
#          FILE: test.sh
#
#         USAGE: ./test.sh
#
#   DESCRIPTION: bash srcipt to run TecIO tests
#===============================================================================
function print_usage {
  echo "Bash srcipt to run TecIO tests"
  echo "Usage:"
  echo "   ./test.sh clean    => clean test folders"
  echo "   ./test.sh run      => run all the tests"
  exit 1
}

runit()
{
  echo "running:" $1 >> $resultsFile
  if test -x ./test* ; then
    ./test* >> $resultsFile
    result=$?
    if [[ $result == "0" ]]; then
      echo -e "${GREEN}success${NC} - $dir" 
    else 
      echo -e "${RED}fail${NC}      - $dir" 
    fi
  else
    echo -e "${RED}not run${NC}   - $dir" 
  fi
}

[ $# == 0 ] && print_usage

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

resultsFile=$(pwd)/testlog.txt
rm -rf $resultsFile > /dev/null 2>&1

dirs=`find . -maxdepth 1 -type d | grep -v '^\.$' | sed -e 's/^\.//' -e 's/^\///'`

if [[ "$1" == "run" ]]; then

  for dir in $dirs
  do
    [[ "$dir" == "flushpartitioned" ]] && continue
    [[ "$dir" == "gridsolution" ]] && continue
    [[ "$dir" == "ij_ordered" ]] && continue
    [[ "$dir" == "ijkpartitioned" ]] && continue
    [[ "$dir" == "multiplefiles" ]] && continue
    [[ "$dir" == "rewriteszl" ]] && continue
    cd $dir                      > /dev/null 2>&1
    ln -sf ../base.make makefile
    make                         > /dev/null 2>&1
    result=$?
    echo "building:" $dir >> $resultsFile
    if test $result -eq 0 ; then
      echo "Passed" >> $resultsFile
      runit
    else
      echo "Failed" >> $resultsFile
    fi
    cd ..                        > /dev/null 2>&1
  done

fi

if [[ "$1" == "clean" ]]; then
  
  for dir in $dirs; do
    cd $dir                      > /dev/null 2>&1
    rm -r build
    rm makefile *.plt test*      > /dev/null 2>&1
    cd ..                        > /dev/null 2>&1
  done
fi
