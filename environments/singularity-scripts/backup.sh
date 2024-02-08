#!/bin/bash
declare -a repos=(
    "verysure/rdkit-django"
    "verysure/postgres-alpine"
    "verysure/miniconda3-rdkit-dftb"
)

for url in "${repos[@]}"
do
    foldername=$(basename "$url")
    rm -rf $foldername
    git clone "https://github.com/$url.git"
    rm -rf $foldername/.git
done

