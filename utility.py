import argparse
from Bio import SeqIO
from pathlib import Path
'''
A file to house various functionalities.

'''

'''
Modify each record based on its content, output file name and index - determined by
its order in the file - the first sequence has index 0, the next 1 and so on.

- currently this just renames sequences
'''
def modify_record(record, index, output_name):
    record.id = output_name + '_' + str(index)
    record.name = ''
    record.description = ''
    
def modify(input_path, output_path):
    file_type = input_path[-5:]
    file_name = Path(output_path).stem
    records =  list(SeqIO.parse(input_path, file_type))
    index = 0
    for record in records:
        modify_record(record, index, file_name)
        index += 1
    SeqIO.write(records, output_path, file_type)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Some utility functionalities.')
    subparsers = parser.add_subparsers(dest='command')

    rename_parser = subparsers.add_parser("modify", help='modify sequences in a fasta/q file.')
    rename_parser.add_argument('-i', '--input', type=str, help='path to input fasta/q file.')
    rename_parser.add_argument('-o', '--output', type=str, help='path to output fasta/q file.')
    
    args = parser.parse_args()
    if args.command == 'modify':
        modify(args.input, args.output)
