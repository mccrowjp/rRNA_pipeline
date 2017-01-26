#!/usr/bin/env python
#
# fastq_filter - filter FASTQ for low quality, and chimeric sequences from uchime_ref
#
# Version: 0.4 (5/21/2016)
#
# Part of rRNA_pipeline - FASTQ filtering, and swarm OTU classification of 16/18S barcodes
#
# Original version: 4/11/2016 John P. McCrow (jmccrow [at] jcvi.org)
# J. Craig Venter Institute (JCVI)
# La Jolla, CA USA
#
import sys, re, os, getopt
import happyfile

dict_chimera_ids = {}
verbose = False
count_total = 0
count_short_seqs = 0
count_long_seqs = 0
count_non_acgt = 0
count_chimeras = 0
count_low_quality = 0
count_passed = 0

def read_chimeras(chimera_file):
    global dict_chimera_ids

    in_handle = happyfile.hopen_or_else(chimera_file)
        
    if verbose:
        print("Reading chimera file: " + chimera_file, file=sys.stderr)

    while 1:
        line = in_handle.readline()
        if not line:
            break
        line = line.rstrip()

        cols = line.split('\t')
        if cols[16] == 'Y':
            dict_chimera_ids[cols[1]] = 1

def filter_line(out_handle, id, seq, qual, min_quality, min_seq_len, max_seq_len):
    global count_total
    global count_short_seqs
    global count_long_seqs
    global count_non_acgt
    global count_chimeras
    global count_low_quality
    global count_passed
    skipline = False
    count_total += 1
    
    if len(seq) < min_seq_len:
        count_short_seqs += 1
        skipline = True
    elif len(seq) > max_seq_len:
        count_long_seqs += 1
        skipline = True
    elif re.search('[^acgtACGT]', seq):
        count_non_acgt += 1
        skipline = True
    elif id in dict_chimera_ids:
        count_chimeras += 1
        skipline = True
    else:
        lastbad = False
        thisbad = False
        for c in qual:
            thisbad = ord(c)-33 < min_quality
            if thisbad and lastbad:
                count_low_quality += 1
                skipline = True
                break
            lastbad = thisbad

    if not skipline:
        count_passed += 1
        print(">" + id + "\n" + seq, file=out_handle)

def filter_fastq(fastq_file, output_file, min_quality, min_seq_len, max_seq_len):
    in_handle = happyfile.hopen_or_else(fastq_file)
    
    if verbose:
        print("Reading FASTQ file: " + fastq_file, file=sys.stderr)

    out_handle = sys.stdout
    if output_file:
        out_handle = happyfile.hopen_write_or_else(output_file)

    if verbose:
        print("Writing FASTA file: " + output_file, file=sys.stderr)

    rnum = 1
    id = ""
    seq = ""
    qual = ""
    while 1:
        line = in_handle.readline()
        if not line:
            break
        line = line.rstrip()
        if rnum == 1:
            id = re.split('\s', line[1:])[0]
        elif rnum == 2:
            seq = line
        elif rnum == 4:
            qual = line
            filter_line(out_handle, id, seq, qual, min_quality, min_seq_len, max_seq_len)
        rnum += 1
        if rnum > 4:
            rnum = 1

def test_all():
    print("[fastq_filter] test_all: passed", file=sys.stderr)

###

def main(argv):
    help = "\n".join([
        "fastq_filter v0.4 (May 21, 2016)",
        "Filter FASTQ file for length, quality, and chimeras (usearch)",
        "",
        "Usage: " + os.path.basename(argv[0]) + " (options)",
        "   -f file        : FASTQ file (required)",
        "   -o file        : output FASTA file (default: stdout)",
        "   -c file        : usearch -uchime_ref output (optional)",
        "   -q int         : minimum quality score (default: 35)",
        "   -m int         : minimim sequence length (default: 50)",
        "   -x int         : maximum sequence length (default: Inf)",
        "   -h, --help     : help",
        "   -v, --verbose  : more information to stderr", ""])

    global verbose
    fastq_file = ""
    chimera_file = ""
    output_file = ""
    min_quality = 30
    min_seq_len = 50
    max_seq_len = float("Inf")
    unused_args = []
    
    try:
        opts, args = getopt.getopt(argv[1:], "f:o:c:q:m:x:hv", ["help", "verbose", "test"])
    except getopt.GetoptError:
        print(help, file=sys.stderr)
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(help, file=sys.stderr)
            sys.exit()
        elif opt == '--test':
            test_all()
            sys.exit()
        elif opt == '-f':
            fastq_file = arg
        elif opt == '-o':
            output_file = arg
        elif opt == '-c':
            chimera_file = arg
        elif opt == '-q':
            min_quality = int(re.sub('=','', arg))
        elif opt == '-m':
            min_seq_len = int(re.sub('=','', arg))
        elif opt == '-x':
            max_seq_len = int(re.sub('=','', arg))
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            unused_args.append(opt)

    if not fastq_file:
        print(help, file=sys.stderr)
        sys.exit(2)

    if len(unused_args) > 0:
        print(help, file=sys.stderr)
        for arg in unused_args:
            print("Unused argument: " + str(arg), file=sys.stderr)
        sys.exit(2)

    if verbose:
        print("\n".join([
            "fastq file:   " + fastq_file,
            "chimera file: " + chimera_file,
            "output file:  " + output_file,
            "min quality:  " + str(min_quality),
            "min seq len:  " + str(min_seq_len),
            "max seq len:  " + str(max_seq_len)]), file=sys.stderr)

    if chimera_file:
        read_chimeras(chimera_file)

    filter_fastq(fastq_file, output_file, min_quality, min_seq_len, max_seq_len)

    if verbose and count_total:
        print("\n".join([
            "seqs total:       " + str(count_total),
            "seqs too short:   " + str(count_short_seqs) + " (" + str(round(100.0*count_short_seqs/count_total, 1)) + "%)",
            "seqs too long:    " + str(count_long_seqs) + " (" + str(round(100.0*count_long_seqs/count_total, 1)) + "%)",
            "seqs bad chars:   " + str(count_non_acgt) + " (" + str(round(100.0*count_non_acgt/count_total, 1)) + "%)",
            "seqs chimera:     " + str(count_chimeras) + " (" + str(round(100.0*count_chimeras/count_total, 1)) + "%)",
            "seqs low quality: " + str(count_low_quality) + " (" + str(round(100.0*count_low_quality/count_total, 1)) + "%)",
            "seqs passed:      " + str(count_passed) + " (" + str(round(100.0*count_passed/count_total, 1)) + "%)"]), file=sys.stderr)

if __name__ == "__main__":
    main(sys.argv)
