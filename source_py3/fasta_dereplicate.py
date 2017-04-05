#!/usr/bin/env python
#
# fasta_dereplicate - Unique FASTA sequences (100% identity) in swarm or besthit formats
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
import gzip, bz2
import happyfile
import hashlib

verbose = False

dict_all_sample_names = {}
dict_sample_name = {}
dict_id_file_counts = {}
dict_id_counts = {}
dict_id_seq = {}
dict_id_map = {}
good_fasta_files = []

class Format:
    swarm = 1
    bestid = 2

def read_sample_names(sample_names_file):
    global dict_sample_name
    
    if sample_names_file:
        in_handle = happyfile.hopen_or_else(sample_names_file)
    
        if verbose:
            print("Reading sample names file: " + sample_names_file, file=sys.stderr)
        
        while 1:
            line = in_handle.readline()
            if not line:
                break
            line = line.rstrip()
            
            name, file = line.split("\t")
            if name in dict_all_sample_names:
                print("Duplicate sample name found: " + name, file=sys.stderr)
                sys.exit(2)
            
            dict_sample_name[file] = name
            dict_all_sample_names[name] = 1
    
            m = re.search('^(.+)\.filtered\.fa$', file)
            if m:
                dict_sample_name[m.group(1)] = name
            else:
                dict_sample_name[file + ".filtered.fa"] = name

        in_handle.close()

def derep_line(id, seq, filenum):
    global dict_id_file_counts
    global dict_id_counts
    global dict_id_seq
    
    if seq:
        seq = seq.lower()
        key = hashlib.sha1(seq).hexdigest()
        dict_id_counts[key] = dict_id_counts.get(key,0) + 1
        dict_id_file_counts[key, filenum] = dict_id_file_counts.get((key, filenum), 0) + 1
        dict_id_seq[key] = seq
        dict_id_map[id] = key

def derep_fasta(fasta_files, min_fasta):
    global good_fasta_files
    filenum = 0
    
    for fasta_file in fasta_files:
        total_seqs = 0
        in_handle = happyfile.hopen_or_else(fasta_file)
        
        if verbose:
            print("Reading FASTA file: " + fasta_file, file=sys.stderr)

        id = ""
        seq = ""
        while 1:
            line = in_handle.readline()
            if not line:
                break
            line = line.rstrip()

            if line.startswith(">"):
                total_seqs += 1
                derep_line(id, seq, filenum)
                id = line[1:]
                seq = ""
            else:
                seq += re.sub('\s', '', line)
        derep_line(id, seq, filenum)
        in_handle.close()
        
        # Remove counts for this file if below minimum
        if total_seqs < min_fasta:
            print("[fasta_dereplicate] Excluding: " + fasta_file, file=sys.stderr)
            for key in dict_id_counts:
                dict_id_counts[key] -= dict_id_file_counts.get((key, filenum), 0)
                dict_id_file_counts[key, filenum] = 0
        else:
            good_fasta_files.append(fasta_file)
            filenum += 1

def write_dereps(output_fasta_file, output_counts_file, output_map_file, id_format, min_samples, min_count):
    dict_bestid = {}
    dict_id_num_samples = {}
    
    for key in dict_id_counts:
        for filenum in range(len(good_fasta_files)):
            if dict_id_file_counts.get((key, filenum), 0) > 0:
                dict_id_num_samples[key] = dict_id_num_samples.get(key, 0) + 1

    out_handle1 = sys.stdout
    if output_fasta_file:
        out_handle1 = happyfile.hopen_write_or_else(output_fasta_file)

    if verbose and output_fasta_file:
        print("Writing FASTA file: " + output_fasta_file, file=sys.stderr)

    if id_format == Format.bestid:
        for id in dict_id_map:
            key = dict_id_map[id]
            if (not key in dict_bestid) and dict_id_counts.get(key, 0) > 0:
                dict_bestid[key] = id

    for key in dict_id_counts:
        if dict_id_num_samples.get(key, 0) >= min_samples and dict_id_counts[key] >= min_count and key in dict_id_seq:
            if id_format == Format.swarm:
                print(">" + key + "_" + str(dict_id_counts[key]) + "\n" + dict_id_seq[key], file=out_handle1)
            elif id_format == Format.bestid and key in dict_bestid:
                print(">" + dict_bestid[key] + "\n" + dict_id_seq[key], file=out_handle1)

    out_handle1.close()

    if output_counts_file:
        out_handle2 = happyfile.hopen_write_or_else(output_counts_file)

        if verbose:
            print("Writing counts file: " + output_counts_file, file=sys.stderr)

        column_names = ['id']
        for file in good_fasta_files:
            if file in dict_sample_name:
                column_names.append(dict_sample_name[file])
            else:
                column_names.append(re.sub('\.filtered\.fa$', '', file))

        print("\t".join(column_names), file=out_handle2)

        for key in dict_id_counts:
            if dict_id_num_samples.get(key, 0) >= min_samples and dict_id_counts[key] >= min_count:
                samplecounts = []
                id = key + "_" + str(dict_id_counts[key])
                if id_format == Format.bestid:
                    id = re.split('\s', dict_bestid[key])[0]
                for filenum in range(len(good_fasta_files)):
                    samplecounts.append(dict_id_file_counts.get((key, filenum), 0))
                print(id + "\t" + "\t".join(str(x) for x in samplecounts), file=out_handle2)

        out_handle2.close()

    if output_map_file:
        out_handle3 = happyfile.hopen_write_or_else(output_map_file)
            
        if verbose:
            print("Writing map file: " + output_map_file, file=sys.stderr)

        for id in sorted(dict_id_map, key=dict_id_map.get):
            key = dict_id_map[id]
            if dict_id_num_samples.get(key, 0) >= min_samples and dict_id_counts[key] >= min_count:
                if id_format == Format.swarm:
                    print(key + "_" + str(dict_id_counts[key]) + "\t" + id, file=out_handle3)
                elif id_format == Format.bestid:
                    print(re.split('\s', dict_bestid[key])[0] + "\t" + id, file=out_handle3)

        out_handle3.close()

def test_derep():
    retval = True
    seq = "acgtcatgcatctagctactacgagcacgatcatcgtagc"
    key = "6db096a7187e871007152ab79c856ec7d50236b6"
    derep_line("testid1", seq, 1)
    derep_line("testid2", seq, 2)
    if dict_id_counts.get(key,0) == 2 and dict_id_file_counts.get((key, 1),0) == 1 and dict_id_map.get("testid2","") == key and dict_id_seq.get(key, "") == seq:
        print("[fasta_dereplicate] test_derep: passed", file=sys.stderr)
    else:
        print("[fasta_dereplicate] test_derep: failed", file=sys.stderr)
        retval = False
    return retval

def test_all():
    if not test_derep():
        sys.exit(2)

###

def main(argv):
    help = "\n".join([
        "fasta_dereplicate v0.4 (May 21, 2016)",
        "Dereplicate FASTA",
        "",
        "Usage: " + os.path.basename(argv[0]) + " (options) [FASTA file(s)...]",
        "   -o file        : output FASTA file (default: stdout)",
        "   -c file        : output sample counts file",
        "   -m file        : output ID map table",
        "   -n file        : sample names file",
        "   -l int         : minimum samples (default: 1)",
        "   -t int         : minimum total count (default: 1)",
        "   -s, --swarm    : output format: swarm (default)",
        "   -b, --bestid   : output format: best ID",
        "   --fasta_min    : minimum sample sequences (default: 100)",
        "   -h, --help     : help",
        "   -v, --verbose  : more information to stderr", ""])

    global verbose
    fasta_files = []
    sample_names_file = ""
    output_fasta_file = ""
    output_counts_file = ""
    output_map_file = ""
    id_format = Format.swarm
    min_count = 1
    min_samples = 1
    min_fasta = 100
    
    try:
        opts, args = getopt.getopt(argv[1:], "o:c:m:n:t:l:sbhv", ["swarm", "bestid", "fasta_min", "help", "verbose", "test"])
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
        elif opt == '-o':
            output_fasta_file = arg
        elif opt == '-c':
            output_counts_file = arg
        elif opt == '-m':
            output_map_file = arg
        elif opt == '-n':
            sample_names_file = arg
        elif opt == '-l':
            min_samples = int(re.sub('=','', arg))
        elif opt == '-t':
            min_count = int(re.sub('=','', arg))
        elif opt in ("-s", "--swarm"):
            id_format = Format.swarm
        elif opt in ("-b", "--bestid"):
            id_format = Format.bestid
        elif opt == '--fasta_min':
            min_fasta = int(re.sub('=','', arg))
        elif opt in ("-v", "--verbose"):
            verbose = True

    if len(args) > 0:
        fasta_files = args
    else:
        print(help, file=sys.stderr)
        sys.exit(2)

    if verbose:
        if len(fasta_files) > 1:
            print("input fasta files:    " + fasta_files[0], file=sys.stderr)
            print("\n".join("                      " + x for x in fasta_files[1:]), file=sys.stderr)
        else:
            print("input fasta file:     " + fasta_files[0], file=sys.stderr)

        print("\n".join([
            "output fasta file:    " + output_fasta_file,
            "output counts file:   " + output_counts_file,
            "output map file:      " + output_map_file,
            "output id format:     " + ("swarm", "bestid")[id_format-1],
            "minimum total counts: " + str(min_count),
            "minimum samples:      " + str(min_samples),
            "minimum sequences:    " + str(min_fasta)]), file=sys.stderr)

    read_sample_names(sample_names_file)
    derep_fasta(fasta_files, min_fasta)
    write_dereps(output_fasta_file, output_counts_file, output_map_file, id_format, min_samples, min_count)

if __name__ == "__main__":
    main(sys.argv)
