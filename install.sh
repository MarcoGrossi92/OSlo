#!/bin/bash

set -e  # Exit on any command failure
set -u  # Treat unset variables as an error

PROGRAM=$(basename "$0")
readonly DIR=$(pwd)
BUILD_DIR="$DIR/build"
VERBOSE=false

function usage() {
    cat <<EOF

Install script for OSlo

Usage:
  $PROGRAM [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]

Global Options:
  -h       , --help         Show this help message and exit
  -v       , --verbose      Enable verbose output

Commands:
  build                     Perform a full build
    --compiler=<name>       Set compilers suit (intel,gnu)
    --use-openmp            Use OpenMP
    --use-mpi               Use MPI
    --use-sundials          Use Sundials

  compile                   Compile the program using the CMakePresets file

EOF
    exit 1
}


log() {
    if [ "$VERBOSE" = true ]; then
        # Bold and dim gray (ANSI escape: bold + color 90)
        echo -e "\033[1;90m$1\033[0m"
    fi
}

error() {
    # Bold red + [ERROR] tag, output to stderr
    echo -e "\033[1;31m[ERROR] $1\033[0m" >&2
}

task() {
    # Bold yellow + ==> tag, output to stdout
    echo -e "\033[1;38;5;186m==> $1\033[0m"
}


# Create default CMakePresets.json if it doesn't exist
function write_presets() {
  FC=$(grep '^CMAKE_Fortran_COMPILER:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  CC=$(grep '^CMAKE_C_COMPILER:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)

  cat <<EOF > CMakePresets.json
{
  "version": 3,
  "cmakeMinimumRequired": {
    "major": 3,
    "minor": 23
  },
  "configurePresets": [
    {
      "name": "default",
      "description": "Default preset",
      "binaryDir": "\${sourceDir}/build",
      "cacheVariables": {
        "CMAKE_BUILD_TYPE": "${BUILD_TYPE}",
        "CMAKE_Fortran_COMPILER": "${FC}",
        "CMAKE_C_COMPILER": "${CC}",
        "USE_OPENMP": "${USE_OPENMP}",
        "USE_MPI": "${USE_MPI}",
        "USE_SUNDIALS": "${USE_SUNDIALS}"
      }
    }
  ]
}
EOF
  log "CMakePresets.json created with default settings."
}


# Default global values
COMMAND=""
COMPILERS=""
BUILD_TYPE="RELEASE"
USE_OPENMP="false"
USE_MPI="false"
USE_SUNDIALS="false"

# Define allowed options for each command using regular arrays
CMD=("build" "compile")
CMD_OPTIONS_build=("--use-openmp" "--use-mpi" "--compilers" "--use-sundials")

# Parse options with getopts
while getopts "hv:-:" opt; do
    case "$opt" in
        -)
            case "$OPTARG" in
                verbose) VERBOSE=true ;;
                help) usage ;;
                *) error "Unknown global option '--$OPTARG'"; usage ;;
            esac
            ;;
        h) usage ;;
        v) VERBOSE=true ;;
        *) error "Unknown global option '-$opt'"; usage ;;
    esac
done
shift $((OPTIND -1))

# Ensure a command was provided
if [[ $# -eq 0 ]]; then
  error "No command provided!"
  usage
fi

COMMAND="$1"
# Check if the command is valid
if [[ ! " ${CMD[@]} " =~ " ${COMMAND} " ]]; then
  error "Unknown command '$COMMAND'"
  usage
fi
shift

# Parse command-specific options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --compilers=*)
            [[ "$COMMAND" == "build" ]] || { error " --compilers is only valid for 'build' command"; exit 1; }
            COMPILERS="${1#*=}"
            ;;
        --use-openmp)
            [[ "$COMMAND" == "build" ]] || { error " --use-openmp is only valid for 'build' command"; exit 1; }
            USE_OPENMP="true"
            ;;
        --use-mpi)
            [[ "$COMMAND" == "build" ]] || { error " --use-mpi is only valid for 'build' command"; exit 1; }
            USE_MPI="true"
            ;;
        --use-sundials)
            [[ "$COMMAND" == "build" ]] || { error " --use-sundials is only valid for 'build' command"; exit 1; }
            USE_SUNDIALS="true"
            ;;
        *)
            eval "opts=(\"\${CMD_OPTIONS_${COMMAND}[@]}\")"
            error "Unknown option '$1' for command '$COMMAND'. Valid options: ${opts[@]}"
            exit 1
            ;;
    esac
    shift
done


# Execute the selected command
case "$COMMAND" in
    build)
        task "Building project"
        log "Use OpenMP: $USE_OPENMP"
        log "Use MPI: $USE_MPI"
        log "Use Sundials: $USE_SUNDIALS"
        rm -rf $BUILD_DIR
        if [[ $COMPILERS == "intel" ]]; then 
            log "Using Intel compilers"
            export FC="ifx"
            export CC="icx"
        elif [[ $COMPILERS == "gnu" ]]; then 
            log "Using GNU compilers"
            export FC="gfortran"
            export CC="gcc"
        fi
        cmake -B $BUILD_DIR -DUSE_OPENMP=$USE_OPENMP -DUSE_MPI=$USE_MPI -DUSE_SUNDIALS=$USE_SUNDIALS -DCMAKE_BUILD_TYPE=$BUILD_TYPE || exit 1
        cmake --build $BUILD_DIR || exit 1

        task "Write CMakePresets.json"
        write_presets
        ;;
    compile)
        task "Compiling project using CMakePresets"
        cmake --preset default || exit 1
        cmake --build $BUILD_DIR || exit 1
        ;;
    *)
        error "Unknown command '$COMMAND'"
        usage
        ;;
esac