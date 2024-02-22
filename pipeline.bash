# \#!/usr/bin/env bash
date

stage=0
nocutoff=false

while getopts t:s:n flag
do
    case "${flag}" in
        t) pipeline_dir=${OPTARG};;
        s) stage=${OPTARG};;
        n) nocutoff=true;;
    esac
done

echo "Processing $pipeline_dir beginning at stage $stage"

if [ "$nocutoff" = true ] ; then
    echo "Processing with no cutoff..."
fi


if [ $stage -le 0 ]; then
    python -m trials.$pipeline_dir.config
fi

if [ $stage -le 1 ]; then
    if [ "$nocutoff" = true ] ; then
        python src/ribobase_counts_processing.py -i "trials/$pipeline_dir/ribo_raw.csv" -r "trials/$pipeline_dir/rnaseq_raw.csv" -m "paired" -o "trials/$pipeline_dir" --cpm_cut_off 0 --overall_cut_off 0
    else
        python src/ribobase_counts_processing.py -i "trials/$pipeline_dir/ribo_raw.csv" -r "trials/$pipeline_dir/rnaseq_raw.csv" -m "paired" -o "trials/$pipeline_dir"
    fi
fi


if [ $stage -le 2 ]; then    
    Rscript src/TE.R "trials/$pipeline_dir"
fi

if [ $stage -le 3 ]; then
    python src/transpose_TE.py -o "trials/$pipeline_dir"
fi

date