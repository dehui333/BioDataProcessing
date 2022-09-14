


'''
To evaluate the quality of reads w.r.t. alignment to trimmed assembly.

1. put reads from different strains into different files using info in the description.
2. align each read set onto the trimmed assembly of that strain.
3. calculate accuracy.
'''









if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluate the quality of reads w.r.t. alignment to trimmed assembly.')
    parser.add_argument('-i', '--reads', type=str, help='Path to the reads.')
    parser.add_argument('-r', '--bam', type=str, help='Path to the sam/bam file with reads to trimmed assembly alignment.')
    args = parser.parse_args()
    evaluate_quality(args.reads, args.bam)