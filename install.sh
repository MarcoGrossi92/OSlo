#!/bin/bash -
#===============================================================================
#
#          FILE: intstall.sh
#
#         USAGE: run "./install.sh [options]" from OSlo master directory
#
#   DESCRIPTION: A utility script that builds OSlo project
#===============================================================================

# DEBUGGING
set -e
set -C # noclobber

# INTERNAL VARIABLES AND INITIALIZATIONS
readonly PROJECT="OSlo"
readonly DIR=$(pwd)
readonly PROGRAM=`basename "$0"`

function usage () {
    echo "Install script of $PROJECT"
    echo "Usage:"
    echo
    echo "$PROGRAM --help|-?"
    echo "    Print this usage output and exit"
    echo
    echo "$PROGRAM --build  |-b"
    echo "    Build the whole project via CMake"
    echo
    echo "$PROGRAM --compile|-c <type>"
    echo "    Compile with build <type> (DEBUG, RELEASE, TESTING)"
    echo
}

function build_project () {
  rm -rf bin build && mkdir -p build
  cd build
  cmake .. -DUSE_OPENMP=ON -DCMAKE_BUILD_TYPE=RELEASE
  make
}

function compile () {
  mkdir -p build
  cd build
  cmake .. -DUSE_OPENMP=ON -DCMAKE_BUILD_TYPE=$TYPE
  make
}

BUILD=0
TYPE=0

# RETURN VALUES/EXIT STATUS CODES
readonly E_BAD_OPTION=254

# PROCESS COMMAND-LINE ARGUMENTS
if [ $# -eq 0 ]; then
  usage
  exit 0
fi

while test $# -gt 0; do
  if [ x"$1" == x"--" ]; then
    # detect argument termination
    shift
    break
  fi
  case $1 in

    --build | -b )
      shift
      BUILD=1
      ;;

    --compile | -c )
      shift
      TYPE="$1"
      ;;

    --setvars | -s )
      shift
      SETVARS=1
      ;;

    -? | --help )
      usage
      exit
      ;;

    -* )
      echo "Unrecognized option: $1" >&2
      usage
      exit $E_BAD_OPTION
      ;;

    * )
      break
      ;;
  esac
done

if [ "$BUILD" != "0" ]; then
  build_project
elif [ "$TYPE" != "0" ]; then
  compile
else
  usage
fi
