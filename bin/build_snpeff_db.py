#!/usr/bin/env python

"""
Build a snpEff database within the current conda environment based on an NCBI accession.
"""

import argparse
import os

def main(args):
    snpeff_dir_prefix = 'snpeff'

    for dir_name in os.scandir('/'.join([os.environ['CONDA_PREFIX'], 'share'])):
        if dir_name.name.startswith(snpeff_dir_prefix):
            destination_data_dir = os.path.join(dir_name.path, 'data', args.accession)
            if not(os.path.exists(destination_data_dir)):
                os.chdir('/'.join([os.environ['CONDA_PREFIX'], 'share', dir_name.name]))
                os.system("/".join([os.environ["CONDA_PREFIX"], "share", dir_name.name, "scripts/buildDbNcbi.sh"]) + " " + args.accession)
            else:
                print("Database dir exists. Skipping database build for " + args.accession)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--accession')
    args = parser.parse_args()
    main(args)

