#!/bin/bash

set -e  # Exit on any command failure
set -u  # Treat unset variables as an error

PROGRAM=$(basename "$0")
readonly DIR=$(pwd)
BUILD_DIR="$DIR/build"
VERBOSE=false
project=ORION

function usage() {
    cat <<EOF

Install script for $project

Usage:
  $PROGRAM [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]

Global Options:
  -h       , --help         Show this help message and exit
  -v       , --verbose      Enable verbose output

Commands:
  build                     Perform the full build
    --compiler=<name>       Set compilers suit (intel,gnu)
    --use-tecio             Use TecIO

  compile                   Compile the program using the CMakePresets file

  setvars                   Set $project paths in environment variables

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


function download_extra () {
  [[ -f doxygen ]] && return
  if [[ $1 == linux ]]; then
    url=https://www.doxygen.nl/files/doxygen-1.10.0.linux.bin.tar.gz
    curl -O $url
    tar zxf *.tgz
  fi
}


function define_path () {
  rm -f .setvars.sh
  echo 'export ORIONDIR='$DIR >> .setvars.sh
  if [[ $SHELL == *"zsh"* ]]; then
    echo 'ORION () { '$DIR'/bin/app/converter $@; }' >> .setvars.sh
    RCFILE=$HOME/.zshrc
  elif [[ $SHELL == *"bash"* ]]; then
    echo 'function ORION () { '$DIR'/bin/app/converter $@; }' >> .setvars.sh
    RCFILE=$HOME/.bashrc
  fi
  log "RC file: $RCFILE"
  echo 'export -f ORION' >> .setvars.sh
  grep -v "ORION" $RCFILE > tmpfile && mv tmpfile $RCFILE
  echo 'source '$DIR'/.setvars.sh' >> $RCFILE
  source $RCFILE --force
}


# Create default CMakePresets.json if it doesn't exist
function write_presets() {
  FC=$(grep '^CMAKE_Fortran_COMPILER:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)
  CXX=$(grep '^CMAKE_CXX_COMPILER:FILEPATH=' "$BUILD_DIR/CMakeCache.txt" | cut -d= -f2-)

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
        "CMAKE_CXX_COMPILER": "${CXX}",
        "USE_TECIO": "${USE_TECIO}"
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
USE_TECIO=false

# Define allowed options for each command using regular arrays
CMD=("build" "compile" "setvars")
CMD_OPTIONS_build=("--compilers --use-tecio")

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
        --compilers=*)
            [[ "$COMMAND" == "build" ]] || { error " --compilers is only valid for 'build' command"; exit 1; }
            COMPILERS="${1#*=}"
            ;;      
        --use-tecio)
            [[ "$COMMAND" == "build" ]] || { error " --use-tecio is only valid for 'build' command"; exit 1; }
            USE_TECIO=true
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

        # download Doxygen
        #./doxygen .Doxyfile

        task "Configuring and building $project"
        rm -rf $BUILD_DIR
        if [[ $COMPILERS == "intel" ]]; then 
            log "Using Intel compilers"
            export FC="ifx"
            export CXX="icpx"
        elif [[ $COMPILERS == "gnu" ]]; then 
            log "Using GNU compilers"
            export FC="gfortran"
            export CXX="g++"
        fi
        log "Build dir: $BUILD_DIR"
        log "Use TecIO: $USE_TECIO"
        log "Build type: $BUILD_TYPE"
        if [[ -z "${FC+x}" || -z "${CXX+x}" ]]; then
          log "Compilers not set. CMake will decide."
        else
          log "Compilers: FC=$FC, CXX=$CXX"
        fi
        cmake -B $BUILD_DIR -DUSE_TECIO=$USE_TECIO -DCMAKE_BUILD_TYPE=$BUILD_TYPE || exit 1
        cmake --build $BUILD_DIR || exit 1
        log "[OK] Compilation successful"

        task "Write CMakePresets.json"
        write_presets
        log "[OK] CMakePresets.json created"

        task "Defining environment variables"
        define_path
        log "[OK] Environment variables defined"
        ;;
    compile)
        task "Compiling $project using CMakePresets"
        cmake --preset default || exit 1
        cmake --build $BUILD_DIR || exit 1
        log "[OK] Compilation successful"
        ;;
    setvars)
        task "Setting $project environment variables"
        define_path
        log "[OK] Environment variables defined"
        ;;
    *)
        error "Unknown command '$COMMAND'"
        usage
        ;;
esac