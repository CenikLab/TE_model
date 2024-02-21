# TE_model
This repository contains code supporting the manuscript titled 'xxxx' Included are scripts such as preprocessing ribosome profiling and RNA-seq data for the calculation of translation efficiency, benchmarking with GO terms or hu.MAP terms, and predicting novel gene functions.

## Overview 

There are 4 main folders in the repo:
- riboflow_scr/ - houses demo riboflow yaml files to generate ribo files
- data/ - houses all the backing data for the pipeline
- src/ - houses all the codes to run the pipeline for TE calculation
- trials/ - houses folders testing different preprocessing techniques
- other_scr/ - houses codes for benchmarking and gene function prediction.

The main code for winsorization is  `src/ribo_counts_to_csv.py`. This is the script that "flattens" the ribo files, preprocessing the counts. If you're looking for the standard capping procedure, capping all reads above the 99.5 percentile, that's in `trials/PAX_cap/config.py`. 

## Getting started
### Files you need to calculate TE

<!-- gene_correction_map_rna.csv - Same as gene_correction_map.csv. Trials that use this only use the 40 column -->

Ribo files required to reproduce these results are not included in the repo for size reasons and need to be requested by email ccenik@austin.utexas.edu, or generated yourself via scripts under riboflow_scr. What's required to run this out of the box is:
- Simlinking the directory of ribo files under `data/ribo/`, example HELA ribo files can be found via [Zenodo link](https://zenodo.org/records/10594392)

### Dependencies
- Python + libraries (pandas, ribopy, numpy, bioinfokit, Bio)
- R + libraries (R packages needed by `TE.R`)

## Workflow
Let's go through the process of trying out a new preprocessing variation. 
1. Create a new directory in the trials directory to house all the data files.
2. Update the sample information in the file `data/paxdb_filtered_sample.csv` to the specific samples you are currently working with.
3. Add a config.py file file to the directory.
4. That config.py file should call `main()` from `src.ribo_counts_to_csv.py`. 
5. Run `bash pipeline.bash -t <DIRECTORY_NAME>`

The bash script will first run your config.py file, and then the rest of the pipeline. You'll know you're done when you see a `model_results.txt` file in the directory. This can take several hours due to the `TE.R` file. See the section below on performance notes for how to speed this up.

### Naming Convention
Name trials that work on the same input experiments with a like prefix. For example, all the trials in this repo currently start with `"PAX_"`. 

### Flattening Arguments

```python
def main(workdir, sample_filter, ribo_dedup, rna_seq_dedup, process_coverage_fn=None, filter=None, custom_experiment_list=None, rnaseq_fn=None)
```
- `workdir` - The directory to spit out the flattened files. Just make it `os.path.dirname(os.path.realpath(__file__))` for the current directory.
- `sample_filter` - By default, this script pulls experiments from `data/paxdb_filtered_sample.csv`. This argument lets you specify a subset of these experiments by applying a filter to the Pandas Dataframe. Using `lambda df: df` includes all experiments from the `paxdb_filtered_sample` file.
- `ribo_dedup` - Boolean, whether to use dedup for the ribo data.
- `rna_seq_dedup` - Boolean, whether to use dedup for the rna_seq data.

### Pipeline Arguments
```
bash pipeline.bash -t <directory_name> [ -s <stage_number> -n ]
```
- `-t` - The name of the directory under `trials` to run the pipeline for. This directory must have a `config.py` file in it that runs the `main` function imported from `src.ribo_counts_to_csv.py` as detailed above. This is the only required argument.
- `-s <n>` - This argument lets you skip ahead to a specified stage number. Useful if the pipeline quits halfway unexpectedly.
- `-n` - This turns off the cutoff in `src/ribobase_counts_processing.py`, ensuring that all genes in the flattened file are processed. 

## Contact
If you have questions about any of this, or if there's any part of this repo you feel is lacking adequate documentation, please email jonathan.j.chacko@gmail.com or yliu5@utexas.edu.
