import numpy as np
from ribopy import Ribo
import ribopy
import pandas as pd
from src.utils import intevl
import glob
import multiprocessing
import concurrent.futures
from tqdm import tqdm

NUM_WORKERS = 48

def get_coverage_from_experiment(experiment, study, ribo_dedup, rna_seq_dedup, process_coverage_fn=None, filter=None, rnaseq_fn=None):
    
    # process ribosome profiling data
    ribo_ribo = Ribo(
        f"./data/ribo/{study}{'_dedup' if ribo_dedup else ''}/ribo/experiments/{experiment}.ribo",
        alias=ribopy.api.alias.apris_human_alias,
    )

    range_lower, range_upper, _ = intevl(ribo_ribo, experiment)

    if not process_coverage_fn:
        ribo_counts = ribo_ribo.get_region_counts(
            "CDS",
            sum_references=False,
            alias=True,
            range_lower=int(range_lower),
            range_upper=int(range_upper),
        )[experiment]
    else:
        coverage = ribo_ribo.get_coverage(
            experiment,
            alias=True,
            range_lower=int(range_lower),
            range_upper=int(range_upper),
        )
        if filter:
            modified_counts = {k: process_coverage_fn(v, k, ribo_ribo) for k, v in coverage.items() if k in filter}
        else:
            modified_counts = {k: process_coverage_fn(v, k, ribo_ribo) for k, v in coverage.items()}
        ribo_counts = pd.Series(modified_counts)
        ribo_counts.name = experiment

    # process rnaseq_counts
    ribo_rnaseq = Ribo(
        f"./data/ribo/{study}{'_dedup' if rna_seq_dedup else ''}/ribo/experiments/{experiment}.ribo",
        alias=ribopy.api.alias.apris_human_alias,
    )
    rnaseq_data = ribo_rnaseq.get_rnaseq().loc[experiment]
    rnaseq_data.index = map(
        ribopy.api.alias.apris_human_alias, rnaseq_data.index)
    rnaseq_counts = rnaseq_data["CDS"].rename(experiment)

    if rnaseq_fn:
        name = rnaseq_counts.name
        rnaseq_counts = {k: rnaseq_fn(v, k, ribo_rnaseq) for k, v in rnaseq_counts.items()}
        rnaseq_counts = pd.Series(rnaseq_counts)
        rnaseq_counts.name = name

    print("Processed " + experiment)
    return ribo_counts, rnaseq_counts


def worker(experiment, study, ribo_dedup, rna_seq_dedup, process_coverage_fn, filter, rnaseq_fn):
    try:
        return get_coverage_from_experiment(experiment, study, ribo_dedup, rna_seq_dedup, process_coverage_fn, filter, rnaseq_fn)
    except Exception as e:
        print(e)

def main(workdir, sample_filter, ribo_dedup, rna_seq_dedup, process_coverage_fn=None, filter=None, custom_experiment_list=None, rnaseq_fn=None):
    if custom_experiment_list:
        experiments = custom_experiment_list
    else:
        all_experiments = pd.read_csv("data/paxdb_filtered_sample.csv", index_col=0)
        filtered_experiments = sample_filter(all_experiments)
        experiments = filtered_experiments['experiment_alias'].tolist()

    # get all files listed as {study: experiment list}
    files = glob.glob("./data/ribo/*/ribo/experiments/*.ribo")
    files = [f for f in files if "dedup" not in f]
    study_experiments = {file.split(
        "/")[3]: [] for file in files if file.split("/")[6].split(".")[0] in experiments}

    for file in files:
        experiment = file.split("/")[6].split(".")[0]
        if experiment in experiments:
            study_experiments[file.split("/")[3]].append(experiment)

    # get data from all files
    executor = concurrent.futures.ProcessPoolExecutor(NUM_WORKERS)
    futures = [executor.submit(worker, experiment, study, ribo_dedup, rna_seq_dedup, process_coverage_fn, filter, rnaseq_fn)
            for study, experiments in study_experiments.items()
            for experiment in experiments]
    concurrent.futures.wait(futures)

    # get_coverage_from_experiment(experiment, study)
    ribo_series = [f.result()[0] for f in futures if f.result() is not None]
    rnaseq_series = [f.result()[1] for f in futures if f.result() is not None]

    combined_ribo = pd.concat(ribo_series, axis=1)
    combined_rnaseq = pd.concat(rnaseq_series, axis=1)

    if filter:
        combined_ribo = combined_ribo[combined_ribo.index.isin(filter)]
        combined_rnaseq = combined_rnaseq[combined_rnaseq.index.isin(filter)]
    
    combined_ribo = combined_ribo.loc[:,~combined_ribo.columns.duplicated()].copy()
    combined_rnaseq = combined_rnaseq.loc[:,~combined_rnaseq.columns.duplicated()].copy()

    combined_rnaseq.to_csv(f"{workdir}/rnaseq_raw.csv")
    combined_ribo.to_csv(f"{workdir}/ribo_raw.csv")
