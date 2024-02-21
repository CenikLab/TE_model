import ribopy
from ribopy.core.get_gadgets import get_region_boundaries, get_reference_names
from functools import cache
import numpy as np
import itertools
import pandas as pd
from src.Fasta import FastaFile

@cache
def intevl(ribo_object, experiment_id):
    # ribo_object = ribo_data("%s.ribo"%experiment_id)
    data_tmp = ribo_object.get_length_dist("CDS")
    data_tmp.reset_index(inplace=True)
    data = data_tmp.iloc[6:26]
    pct_85 = sum(data["%s" % experiment_id]) * 0.85
    # pct_90=sum(data["%s"%experiment_id])*0.90
    value = data[data["%s" % experiment_id] == data["%s" % experiment_id].max()][
        "%s" % experiment_id
    ].values[0]
    mmin = mmax = data[data["%s" % experiment_id] == data["%s" % experiment_id].max()][
        "read_length"
    ].values[0]
    while value <= pct_85:
        if mmax < 40 and mmin > 21:
            if (
                data[data["read_length"] == mmax + 1]["%s" % experiment_id].values[0]
                >= data[data["read_length"] == mmin - 1]["%s" % experiment_id].values[0]
            ):
                mmax += 1
                value += data[data["read_length"] == mmax]["%s" % experiment_id].values[
                    0
                ]
            else:
                mmin -= 1
                value += data[data["read_length"] == mmin]["%s" % experiment_id].values[
                    0
                ]
        elif mmax == 40:
            mmin -= 1
            value += data[data["read_length"] == mmin]["%s" % experiment_id].values[0]
        elif mmin == 21:
            mmax += 1
            value += data[data["read_length"] == mmax]["%s" % experiment_id].values[0]
    # print(min,max)
    read_pct = value / sum(data["%s" % experiment_id])
    return mmin, mmax, read_pct

@cache
def get_cds_range_lookup(ribo):
    """
    Create a dict of gene to region ranges, so that the CDS range can be found for a given experiment.
    """
    names = get_reference_names(ribo._handle)
    if ribo.alias != None:
        names = map(ribo.alias.get_alias, names)
    boundaries = get_region_boundaries(ribo._handle)
    boundary_lookup = dict(zip(list(names), boundaries))
    return boundary_lookup


def cap_outliers(arr, thresh, filter_zeros=True, cap_min=0):
    arr = arr.copy()
    if filter_zeros:
        arr = arr[arr > 0]
    if len(arr) == 0: return arr
    cap = np.percentile(arr, thresh)
    if cap < cap_min:
        cap = cap_min
    arr[arr > cap] = cap
    return arr

def cap_outliers_cds_only(arr, gene, boundary_lookup, thresh=99, filter_zeros=False):
    start, end = boundary_lookup[gene][1]
    arr = arr.copy()
    arr = arr[start:end]
    if filter_zeros:
        arr = arr[arr > 0]
    if len(arr) == 0: return arr
    cap = np.percentile(arr, thresh)
    arr[arr > cap] = cap
    return arr

class BiasCorrection:
    def __init__(self, ribo, experiment, f3_length=3, f5_length=3):
        self.f3_length = f3_length
        self.f5_length = f5_length
        self.ribo = ribo
        self.experiment = experiment

        self.start_nmer_counts = self.init_nmer_counts(f5_length)
        self.start_baseline_counts = self.init_nmer_counts(f5_length)
        self.end_nmer_counts = self.init_nmer_counts(f3_length)
        self.end_baseline_counts = self.init_nmer_counts(f3_length)

        self.range_lookup = get_cds_range_lookup(ribo)
        self.intevl_start, self.intevl_end, _ = intevl(ribo, experiment)
        self.footprints = []

        fasta = FastaFile("data/appris_human_v2_selected.fa.gz")
        fasta_dict = {e.header: e.sequence for e in fasta}
        self.sequence_dict = {
            ribopy.api.alias.apris_human_alias(transcript): fasta_dict[transcript] for transcript in ribo.transcript_names
        }

    def process_gene(self, gene, arr=None):
        df = self.ribo.get_transcript_coverage(
            self.experiment, alias=(self.ribo.alias != None), transcript=gene
        )
        min_length = self.ribo.minimum_length
        arr = df.to_numpy()[self.intevl_start -
                            min_length: self.intevl_end + 1 - min_length]


        footprints = []
        max_count = -1
        total = 0
        start_cds, end_cds = self.range_lookup[gene][1]
        for i, row in enumerate(arr):
            for idx, count in enumerate(row):
                # make sure footprint is in cds range
                if idx >= start_cds and idx < end_cds:
                    read_length = i + self.intevl_start
                    start_nmer = self.get_start_nmer(gene, idx)
                    end_nmer = self.get_end_nmer(gene, idx, read_length)
                    # make sure footprint isn't on the edge of the range
                    if len(start_nmer) == self.f5_length and len(end_nmer) == self.f3_length:
                        total += count
                        if count > max_count:
                            max_count = count
                        footprints.append({
                            "idx": idx,
                            "COUNT": count,
                            "START_NMER": self.get_start_nmer(gene, idx),
                            "END_NMER": self.get_end_nmer(gene, idx, read_length),
                            "READ_LENGTH": read_length,
                            "gene": gene
                        })

        counts = [x['COUNT'] for x in footprints]
        mean = np.mean(counts)
        max = np.max(counts)
        p10 = np.percentile(counts, 10)
        p25 = np.percentile(counts, 25)
        p50 = np.percentile(counts, 50)
        p75 = np.percentile(counts, 75)
        p90 = np.percentile(counts, 90)
        p99 = np.percentile(counts, 99)
        std = np.std(counts)
                        
        # normalize by dividing my the max count
        for i in range(len(footprints)):
            # footprints[i]['COUNT'] *= 10000
            # footprints[i]['COUNT'] /= max_count
            footprints[i]['MEAN'] = mean
            footprints[i]['MAX'] = max
            footprints[i]['STD'] = std
            footprints[i]['P10'] = p10
            footprints[i]['P25'] = p25
            footprints[i]['P50'] = p50
            footprints[i]['P75'] = p75
            footprints[i]['P90'] = p90
            footprints[i]['P99'] = p99


        self.footprints.extend(footprints)
        return footprints

    def get_start_nmer(self, gene, i):
        n = self.f5_length
        transcript = self.sequence_dict[gene]
        return transcript[i:i+n]

    def get_end_nmer(self, gene, i, read_length):
        n = self.f3_length
        transcript = self.sequence_dict[gene]
        return transcript[i + read_length - n: i + read_length]

    def to_df(self):
        return pd.DataFrame(self.footprints)

    def add_dummy_footprints(self):
        footprints = []
        f5s = list(map(lambda x: "".join(x), itertools.product(['A', 'G', 'T', 'C'], repeat=self.f5_length)))
        f3s = list(map(lambda x: "".join(x), itertools.product(['A', 'G', 'T', 'C'], repeat=self.f3_length)))

        for start_nmer in f5s:
            footprints.append({
                "idx": -1,
                "COUNT": 0,
                "START_NMER": start_nmer,
                "END_NMER": 'A' * self.f5_length,
                "READ_LENGTH": 0,
                "gene": None
            })
            
        for end_nmer in f3s:
            footprints.append({
                "idx": -1,
                "COUNT": 0,
                "START_NMER": 'A' * self.f3_length,
                "END_NMER": end_nmer,
                "READ_LENGTH": 0,
                "gene": None
            })

        self.footprints.extend(footprints)

        
    @staticmethod
    def init_nmer_counts(n):
        return {k: 0 for k in map(lambda x: "".join(x), itertools.product(['A', 'G', 'T', 'C'], repeat=n))}
