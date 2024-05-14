import argparse
import csv

def add_nextclade_version(input_file, output_file, dataset_version, nextclade_version):

    with open(input_file, 'r', newline='') as f_in:

        reader = csv.reader(f_in, delimiter='\t')
    
        header = next(reader)
        
        header.append('nextclade_version')
        header.append('nextclade_dataset_version')
        
        
        with open(output_file, 'w', newline='') as f_out:

            writer = csv.writer(f_out, delimiter='\t')

            writer.writerow(header)
            
            for row in reader:
 
                row.append(nextclade_version)
                row.append(dataset_version)
   
                writer.writerow(row)

def main(args):
    add_nextclade_version(args.nextclade_tsv, args.output_name, args.dataset_version, args.nextclade_version)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nextclade_tsv")
    parser.add_argument("--output_name")
    parser.add_argument("--dataset_version")
    parser.add_argument("--nextclade_version")
    args = parser.parse_args()
    main(args)