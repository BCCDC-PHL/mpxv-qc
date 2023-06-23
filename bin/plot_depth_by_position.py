#!/usr/bin/env python3

import argparse
import csv
import glob
import os

from pathlib import Path
from itertools import islice

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from matplotlib.backends.backend_pdf import PdfPages


def parse_per_base_coverage(per_base_coverage_path: Path):
    """
    """
    per_base_coverage = pd.read_csv(per_base_coverage_path, sep='\t')

    return per_base_coverage


def batched(iterable, n):
    """
    Batch data into tuples of length n. The last batch may be shorter.
    """
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while (batch := tuple(islice(it, n))):
        yield batch
    

def main(args):
    if not os.path.exists(args.directory):
        print("Error, directory does not exist:", args.directory)
        exit(-1)

    all_libraries_per_base_coverage = None
    library_ids = []
    for per_base_coverage_path in glob.glob(os.path.join(args.directory, '*.per_base_coverage.bed')):
        filename = os.path.basename(per_base_coverage_path)
        library_id = filename.split('.')[0]
        library_ids.append(library_id)
        per_base_coverage = parse_per_base_coverage(per_base_coverage_path)
        per_base_coverage['library_id'] = library_id
        if all_libraries_per_base_coverage is None:
            all_libraries_per_base_coverage = per_base_coverage
        else:
            all_libraries_per_base_coverage = pd.concat([all_libraries_per_base_coverage, per_base_coverage], ignore_index=True)

    library_ids = sorted(library_ids)
    num_libraries = len(library_ids)
    max_plots_per_page = 8
    plots_per_page = min(max_plots_per_page, num_libraries)
    page_height = 16
    plot_width = 8.5
    plot_height = page_height / plots_per_page
    plot_aspect = plot_width / plot_height

    with PdfPages(args.output) as pdf:
        sns.set_style("whitegrid", {'grid.linestyle': '--'})
        sns.set_context("notebook", rc={"grid.linewidth": 0.25})
        for library_batch in batched(library_ids, plots_per_page):
            current_batch_libraries_per_base_coverage = all_libraries_per_base_coverage.loc[all_libraries_per_base_coverage['library_id'].isin(library_batch)]
            grid = sns.FacetGrid(current_batch_libraries_per_base_coverage, col="library_id", col_wrap=1, height=plot_height, aspect=plot_aspect)
            grid.map_dataframe(sns.lineplot, x='position', y='depth', color='black', linewidth=0.25)
            plt.yscale('log')
            plt.ylim(1, args.y_axis_limit)
            pdf.savefig()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', type=Path, help="Directory containing per_base_coverage.bed files")
    parser.add_argument('-o', '--output', type=Path, default="./per_base_coverage.pdf", help="")
    parser.add_argument('-y', '--y-axis-limit', type=int, default=10000, help="")
    args = parser.parse_args()
    main(args)
