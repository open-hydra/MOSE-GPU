#!/bin/bash

set -e  # Exit on any command failure
set -u  # Treat unset variables as an error

PROGRAM=$(basename "$0")
readonly DIR=$(pwd)
BUILD_DIR="$DIR/build"
VERBOSE=false
project=FLINT

function usage() {
    cat <<EOF

Install script for $project

Usage:
  $PROGRAM [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]

Global Options:
  -h       , --help         Show this help message and exit
  -v       , --verbose      Enable verbose output

Commands:
  build                     Perform a full build
    --master=<name>         Set master (None, hydra)
    --compilers=<name>      Set compilers suit (intel,gnu)
    --use-sundials          Use Sundials
    --use-cantera           Use Cantera
    --use-tecio             Use TecIO

  compile                   Compile the program using the CMakePresets file

  update                    Download git submodules
    --remote                Use the latest remote commit

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
        "MASTER": "${MASTER_TYPE}",
        "CMAKE_Fortran_COMPILER": "${FC}",
        "CMAKE_C_COMPILER": "${CC}",
        "USE_CANTERA": "${USE_CANTERA}",
        "USE_SUNDIALS": "${USE_SUNDIALS}",
        "USE_TECIO": "${USE_TECIO}"
      }
    }
  ]
}
EOF
}


# Default global values
COMMAND=""
COMPILERS=""
MASTER_TYPE=""
BUILD_TYPE="RELEASE"
USE_SUNDIALS="false"
USE_CANTERA="false"
USE_TECIO="false"
REMOTE=false

# Define allowed options for each command using regular arrays
CMD=("build" "compile" "update")
CMD_OPTIONS_build=("--master --compilers --use-sundials --use-cantera --use-tecio")
CMD_OPTIONS_update=("--remote")

# Parse global options
while getopts "hv-:" opt; do
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
        ?) error "Unknown global option '-$OPTARG'"; usage ;;
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
        --master=*)
            [[ "$COMMAND" == "build" ]] || { error " --master is only valid for 'build' command"; exit 1; }
            if [[ ! "$1" =~ ^--master=(None|hydra)$ ]]; then
                error "Invalid value for --master. Valid values are 'None' or 'hydra'."
                exit 1
            fi
            MASTER_TYPE="${1#*=}"
            ;;
        --compilers=*)
            [[ "$COMMAND" == "build" ]] || { error " --compilers is only valid for 'build' command"; exit 1; }
            COMPILERS="${1#*=}"
            ;;
        --use-cantera)
            [[ "$COMMAND" == "build" ]] || { error " --use-cantera is only valid for 'build' command"; exit 1; }
            USE_CANTERA="true"
            ;;
        --use-sundials)
            [[ "$COMMAND" == "build" ]] || { error " --use-sundials is only valid for 'build' command"; exit 1; }
            USE_SUNDIALS="true"
            ;;
        --use-tecio)
            [[ "$COMMAND" == "build" ]] || { error " --use-tecio is only valid for 'build' command"; exit 1; }
            USE_TECIO="true"
            ;;
        --remote)
            [[ "$COMMAND" == "update" ]] || { error " --remote is only valid for 'update' command"; exit 1; }
            REMOTE=true
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
        task "Building $project"
        if [[ -z "$MASTER_TYPE" ]]; then
            error " --master is required for the 'build' command!"
            exit 1
        fi

        task "Cloning submodules"
        if [[ $MASTER_TYPE == "None" ]]; then
          git submodule update --init lib/OSlo lib/ORION
        fi

        if [[ $COMPILERS == "intel" ]]; then 
            export FC="ifx"
            export CXX="icpx"
            export CC="icx"
        elif [[ $COMPILERS == "gnu" ]]; then 
            export FC="gfortran"
            export CXX="g++"
            export CC="gcc"
        fi
        log "Build dir: $BUILD_DIR"
        log "Build type: $BUILD_TYPE"
        log "Master: $MASTER_TYPE"
        log "Use Cantera: $USE_CANTERA"
        log "Use Sundials: $USE_SUNDIALS"
        log "Use TecIO: $USE_TECIO"
        if [[ -z "${FC+x}" || -z "${CC+x}" || -z "${CXX+x}" ]]; then
          log "Compilers not set. CMake will decide."
        else
          log "Compilers: FC=$FC, CC=$CC", CXX=$CXX
        fi
        rm -rf $BUILD_DIR
        cmake -B $BUILD_DIR -DMASTER=$MASTER_TYPE -DUSE_CANTERA=$USE_CANTERA -DUSE_SUNDIALS=$USE_SUNDIALS -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DUSE_OPENMP=OFF -DUSE_MPI=OFF -DUSE_TECIO=$USE_TECIO || exit 1
        cmake --build $BUILD_DIR || exit 1
        log "[OK] Compilation successful"

        task "Write CMakePresets.json"
        write_presets
        log "[OK] CMakePresets.json created"
        ;;
    compile)
        task "Compiling $project using CMakePresets"
        cmake --preset default || exit 1
        cmake --build $BUILD_DIR || exit 1
        log "[OK] Compilation successful"
        ;;
    update)
        task "Updating git submodules"
        if [[ "$REMOTE" == "true" ]]; then
            log "Updating submodules to latest remote commit"
            git submodule update --init --remote
        else
            log "Updating submodules to current commit"
            git submodule update --init
        fi
        log "[OK] Submodules updated"
        ;;
    *)
        error "Unknown command '$COMMAND'"
        usage
        ;;
esac