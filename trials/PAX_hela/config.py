from src.ribo_counts_to_csv import main
from src.utils import get_cds_range_lookup, cap_outliers_cds_only
import os

workdir = os.path.dirname(os.path.realpath(__file__))
sample_filter = lambda df: df
ribo_dedup = False
rna_seq_dedup = True

def process_coverage_fn(coverage, gene, ribo):
    boundary_lookup = get_cds_range_lookup(ribo)
    return cap_outliers_cds_only(coverage, gene, boundary_lookup, 99.5).sum()

main(workdir, sample_filter, ribo_dedup, rna_seq_dedup, process_coverage_fn)