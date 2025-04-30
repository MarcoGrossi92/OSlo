#!/bin/bash

set -e  # Exit on any command failure
set -u  # Treat unset variables as an error

PROGRAM=$(basename "$0")
readonly DIR=$(pwd)
VERBOSE=false
CMD_OPTIONS=()  # Replace with a regular array

function usage() {
    cat <<EOF

Install script for OSlo

Usage:
  $PROGRAM [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]

Global Options:
  -v       , --verbose         Enable verbose output

Commands:
  build                        Perform the full build

  compile                      Compile the program
    --build-type=<build>       Set build type (release, debug, testing, default: release)

EOF
    exit 1
}


function log() {
    if [ "$VERBOSE" = true ]; then
        echo "$1"
    fi
}


# Default global values
COMMAND=""
BUILD_TYPE=""

# Define allowed options for each command using regular arrays
CMD_OPTIONS_COMPILE=("--build-type")

# Parse options with getopts
while getopts "v:-:" opt; do
    case "$opt" in
        -)
            case "$OPTARG" in
                verbose) VERBOSE=true ;;
                *) echo "Error: Unknown global option '--$OPTARG'"; usage ;;
            esac
            ;;
        v) VERBOSE=true ;;
        *) echo "Error: Unknown global option '-$opt'"; usage ;;
    esac
done
shift $((OPTIND -1))

# Ensure a command was provided
if [[ $# -eq 0 ]]; then
    echo "Error: No command provided!"
    usage
fi

COMMAND="$1"
shift

# Parse command-specific options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --build-type=*)
            [[ "$COMMAND" == "compile" ]] || { echo "Error: --build-type is only valid for 'compile' command"; exit 1; }
            BUILD_TYPE="${1#*=}"
            ;;
        *)
            echo "Error: Unknown option '$1' for command '$COMMAND'. Valid options: ${CMD_OPTIONS[$COMMAND]}"
            exit 1
            ;;
    esac
    shift
done


# Execute the selected command
case "$COMMAND" in
    build)
        log "Building project"
        rm -rf bin build && mkdir -p build
        cd $DIR/build
        cmake ..
        make
        ;;
    compile)
        if [[ -z "$BUILD_TYPE" ]]; then
            echo "Error: --build-type is required for 'compile' command!"
            exit 1
        fi
        log "Compiling with build type: $BUILD_TYPE"
        cd $DIR/build
        cmake .. -DCMAKE_BUILD_TYPE=$BUILD_TYPE
        make
        ;;
    *)
        echo "Error: Unknown command '$COMMAND'"
        usage
        ;;
esac