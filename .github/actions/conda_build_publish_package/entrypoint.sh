#!/bin/bash

set -ex
set -o pipefail

go_to_build_dir() {
    if [ ! -z $INPUT_SUBDIR ]; then
        cd $INPUT_SUBDIR
    fi
}

check_if_meta_yaml_file_exists() {
    if [ ! -f meta.yaml ]; then
        echo "meta.yaml must exist in the directory that is being packaged and published."
        exit 1
    fi
}

build_package(){

    IFS=' ' read -ra PYTHON <<< "$INPUT_PYTHON"
    IFS=' ' read -ra PLATFORMS <<< "$INPUT_PLATFORMS"

    for python in "${PYTHON[@]}"; do 
        conda build -c conda-forge -c bioconda --output-folder . --python $python . 
    done
    for platform in "${PLATFORMS[@]}"; do 
        for filename in /$platform/*.tar.bz2; do
            conda convert /$platform/$filename -p $platform linux-64/*.tar.bz2 -o .
        done 
    done
}

upload_package(){

    IFS=' ' read -ra PLATFORMS <<< "$INPUT_PLATFORMS"

    export ANACONDA_API_TOKEN=$INPUT_ANACONDATOKEN

    for platform in "${PLATFORMS[@]}"; do 
        for filename in ./$platform/*.tar.bz2; do
            anaconda upload $filename
        done
    done 
}

go_to_build_dir
check_if_meta_yaml_file_exists
build_package
upload_package
