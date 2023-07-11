import Bio.SeqIO.FastaIO as FastaIO
from Bio.Seq import Seq
import gzip
import hashlib


def stream_fa(infile):
    if infile.endswith('fasta.gz') or infile.endswith('fa.gz'):
        with gzip.open(infile, 'rt') as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                yield (header, sequence)
    elif infile.endswith('fasta') or infile.endswith('fa'):
        with open(infile, 'rt') as handle:
            for (header, sequence) in FastaIO.SimpleFastaParser(handle):
                yield (header, sequence)
    else:
        raise Exception(f'{infile} not a sequence file.')


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Contig/Scaffold/Transcript filter')
    parser.add_argument('samplename', type=str, help='name of sample')
    parser.add_argument('seqtype', type=str, help='Type of input sequence', choices=['contigs', 'scaffolds', 'transcripts'])
    parser.add_argument('infile', type=str, help='Input Sequence file. Either fasta or fasta.gz.')
    parser.add_argument('outprefix', type=str, help='Prefix for output files.')
    args = parser.parse_args()


    samplename = args.samplename
    seqtype = args.seqtype[:-1]
    filters = [0, 500, 1000]
    infile = args.infile
    outprefix = args.outprefix + '/' + samplename
    sequences = []

    for cnt, (header, sequence) in enumerate(stream_fa(infile), 1):
        sequence = sequence.upper()
        seqlen = len(sequence)
        sequence_rev = str(Seq(sequence).reverse_complement())
        md5_fw = hashlib.md5(sequence.encode()).hexdigest()
        md5_rev = hashlib.md5(sequence_rev.encode()).hexdigest()
        seqname = f'{samplename}_{cnt} '
        #seqname = f'{samplename}-{seqtype}_{cnt} length={seqlen} orig={header}'
        sequences.append((seqname, sequence, md5_fw, md5_rev, seqlen))

    for filtersize in filters:
        with open(f'{outprefix}.{seqtype}s.min{filtersize}.fasta', 'w') as handle:
            for (seqname, sequence, md5_fw, md5_rev, seqlen) in sequences:
                if seqlen >= filtersize:
                    handle.write(f'>{seqname}\n{sequence}\n')
    with open(f'{outprefix}.{seqtype}s.hashes', 'w') as handle:
        for (seqname, sequence, md5_fw, md5_rev, seqlen) in sequences:
            if seqlen >= 500:
                handle.write(f'{seqname}\t{md5_fw}\t{md5_rev}\t{seqlen}\n')


if __name__ == '__main__':
    main()