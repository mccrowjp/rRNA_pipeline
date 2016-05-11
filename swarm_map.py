#!/usr/bin/env python
#
# swarm_map - swarm OTUs constructed and output swarm FASTA, map, and sample counts
#
# Version: 0.2 (5/3/2016)
#
# Part of rRNA_pipeline - FASTQ filtering, and swarm OTU classification of 16/18S barcodes
#
# Original version: 4/11/2016 John P. McCrow (jmccrow [at] jcvi.org)
# J. Craig Venter Institute (JCVI)
# La Jolla, CA USA
#
import sys, re, os, getopt
import happyfile

verbose = False

sample_list = []
dict_sample_name = {}
dict_id_sample_counts = {}
dict_id_counts = {}
dict_swarm_sample_counts = {}
dict_swarm_counts = {}
dict_swarm_seq = {}
dict_id_swarm = {}
dict_swarm_num_samples = {}

def read_sample_names(sample_names_file):
    global dict_sample_name
    
    if sample_names_file:
        in_handle = happyfile.hopen_or_else(sample_names_file)
    
        if verbose:
            print >>sys.stderr, "Reading sample names file: " + sample_names_file
        
        while 1:
            line = in_handle.readline()
            if not line:
                break
            line = line.rstrip()
            
            name, file = line.split("\t")
            dict_sample_name[file] = name
    
            m = re.search('^(.+)\.filtered\.fa$', file)
            if m:
                dict_sample_name[m.group(1)] = name
            else:
                dict_sample_name[file + ".filtered.fa"] = name

        in_handle.close()

def calc_swarm_counts():
    global dict_swarm_counts
    global dict_swarm_sample_counts
    global dict_swarm_num_samples
    
    for id in dict_id_swarm:
        swarm_id = dict_id_swarm[id]
        for i in range(len(sample_list)):
            dict_swarm_sample_counts[swarm_id, i] = dict_swarm_counts.get(swarm_id, 0) + dict_id_sample_counts.get((id, i),0)
            dict_swarm_counts[swarm_id] = dict_swarm_counts.get(swarm_id, 0) + dict_id_sample_counts.get((id, i),0)

    for swarm_id in dict_swarm_counts:
        for i in range(len(sample_list)):
            if (swarm_id, i) in dict_swarm_sample_counts:
                dict_swarm_num_samples[swarm_id] = dict_swarm_num_samples.get(swarm_id, 0) + 1


def read_counts(counts_file):
    global dict_id_counts
    global dict_id_sample_counts
    global sample_list
    
    if counts_file:
        in_handle = happyfile.hopen_or_else(counts_file)
        
        if verbose:
            print >>sys.stderr, "Reading counts file: " + counts_file
        
        firstline = 1
        while 1:
            line = in_handle.readline()
            if not line:
                break
            line = line.rstrip()
            
            if firstline:
                sample_list = line.split("\t")[1:]
            else:
                cols = line.split("\t")
                for i in range(1, len(cols)):
                    dict_id_sample_counts[cols[0], i-1] = int(cols[i])
                    dict_id_counts[cols[0]] = dict_id_counts.get(cols[0], 0) + int(cols[i])
            
            firstline = 0
                
        in_handle.close()

        calc_swarm_counts()

def get_swarms(fasta_file, swarm_file, cpus):
    global dict_id_swarm
    
    if fasta_file and not os.path.exists(swarm_file):
        print >>sys.stderr, "[swarm_map] running swarm"

        if cpus < 1:
            cpus = 1
        
        cmd = " ".join(["swarm -f -t", str(cpus), "-o", swarm_file, fasta_file])
        
        if verbose:
            print >>sys.stderr, cmd
        else:
            cmd += " &>/dev/null"
        
        rc = os.system(cmd)
        if rc != 0:
            print >>sys.stderr, "[swarm_map] ERROR: swarm"
            sys.exit(2)

    in_handle1 = happyfile.hopen_or_else(fasta_file)
    while 1:
        line = in_handle1.readline()
        if not line:
            break
        line = line.rstrip()

        if line.startswith(">"):
            id = line[1:]
            # set any IDs not returned by swarm, to their own cluster
            dict_id_swarm[id] = id
    in_handle1.close()

    in_handle2 = happyfile.hopen_or_else(swarm_file)
    if verbose:
        print >>sys.stderr, "Reading swarm file: " + swarm_file
    
    while 1:
        line = in_handle2.readline()
        if not line:
            break
        line = line.rstrip()

        id_list = re.split('\s', line)
        for id in id_list:
            dict_id_swarm[id] = id_list[0]
    in_handle2.close()


def read_swarm_fasta(fasta_file):
    in_handle = happyfile.hopen_or_else(fasta_file)
    
    if verbose:
        print >>sys.stderr, "Reading FASTA file: " + fasta_file
        
    id = ""
    seq = ""
    while 1:
        line = in_handle.readline()
        if not line:
            break
        line = line.rstrip()
        
        if line.startswith(">"):
            if seq:
                dict_swarm_seq[id] = seq
            id = line[1:]
            seq = ""
        else:
            seq += re.sub('\s', '', line)

    if seq:
        dict_swarm_seq[id] = seq
    in_handle.close()

def write_swarms(output_fasta_file, output_counts_file, output_map_file, min_samples, min_count):
    # set at least one sample where counts not given
    for swarm_id in dict_swarm_counts:
        if not swarm_id in dict_swarm_num_samples:
            dict_swarm_num_samples[swarm_id] = 1
    
    out_handle1 = sys.stdout
    if output_fasta_file:
        out_handle1 = happyfile.hopen_write_or_else(output_fasta_file)

    if verbose and output_fasta_file:
        print >>sys.stderr, "Writing FASTA file: " + output_fasta_file

    for swarm_id in dict_swarm_counts:
        if dict_swarm_num_samples[swarm_id] >= min_samples and dict_swarm_counts[swarm_id] >= min_count:
            print >>out_handle1, ">" + swarm_id + "\n" + dict_swarm_seq[swarm_id]

    out_handle1.close()

    if output_counts_file:
        out_handle2 = happyfile.hopen_write_or_else(output_counts_file)

        if verbose:
            print >>sys.stderr, "Writing counts file: " + output_counts_file

        column_names = ['id']
        for name in sample_list:
            if name in dict_sample_name:
                column_names.append(dict_sample_name[name])
            else:
                column_names.append(name)

        print >>out_handle2, "\t".join(column_names)

        for swarm_id in dict_swarm_counts:
            if dict_swarm_num_samples[swarm_id] >= min_samples and dict_swarm_counts[swarm_id] >= min_count:
                samplecounts = []
                for i in range(len(sample_list)):
                    samplecounts.append(dict_swarm_sample_counts.get((swarm_id, i), 0))
                print >>out_handle2, swarm_id + "\t" + "\t".join(str(x) for x in samplecounts)

        out_handle2.close()

    if output_map_file:
        out_handle3 = happyfile.hopen_write_or_else(output_map_file)
            
        if verbose:
            print >>sys.stderr, "Writing map file: " + output_map_file

        for id in sorted(dict_id_swarm, key=dict_id_swarm.get):
            swarm_id = dict_id_swarm[id]
            if dict_swarm_num_samples[swarm_id] >= min_samples and dict_swarm_counts[swarm_id] >= min_count:
                print >>out_handle3, swarm_id + "\t" + id

        out_handle3.close()

def test_all():
    print >>sys.stderr, "[swarm_map] test_all: passed"

###

def main(argv):
    help = "\n".join([
        "swarm_map v0.2 (May 3, 2016)",
        "swarm OTU FASTA, sample counts, and id mapping",
        "",
        "Usage: "+argv[0]+" (options)",
        "   -f file        : dereplicated FASTA (required)",
        "   -s file        : swarm file (required)",
        "   -d file        : dereplicated counts table (required if -c)",
        "   -o file        : output FASTA file (default: stdout)",
        "   -c file        : output swarm OTU counts file (requires -d)",
        "   -m file        : output ID map table",
        "   -n file        : sample names file",
        "   -l int         : minimum samples (default: 1, requires -d if > 1)",
        "   -t int         : minimum total count (default: 1)",
        "   -x, --cpus int : number of processes to run swarm (default: 1)",
        "   -h, --help     : help",
        "   -v, --verbose  : more information to stderr", ""])

    global verbose
    fasta_file = ""
    swarm_file = ""
    counts_file = ""
    sample_names_file = ""
    output_fasta_file = ""
    output_counts_file = ""
    output_map_file = ""
    min_count = 1
    min_samples = 1
    cpus = 1
    
    try:
        opts, args = getopt.getopt(argv[1:], "f:s:d:o:c:m:n:t:l:x:hv", ["cpus=", "help", "verbose", "test"])
    except getopt.GetoptError:
        print >>sys.stderr, help
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print >>sys.stderr, help
            sys.exit()
        elif opt == '--test':
            test_all()
            sys.exit()
        elif opt == '-f':
            fasta_file = arg
        elif opt == '-s':
            swarm_file = arg
        elif opt == '-d':
            counts_file = arg
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
        elif opt in ("-x", "--cpus"):
            cpus = int(re.sub('=','', arg))
        elif opt in ("-v", "--verbose"):
            verbose = True

    if not (fasta_file and swarm_file):
        print >>sys.stderr, help
        sys.exit(2)

    if (output_counts_file or min_samples > 1) and not counts_file:
        print >>sys.stderr, help + "\nDereplicated counts table required (-d)"
        sys.exit(2)

    if verbose:
        print >>sys.stderr, "input fasta file:     " + fasta_file

        print >>sys.stderr, "\n".join([
            "output fasta file:    " + output_fasta_file,
            "output counts file:   " + output_counts_file,
            "output map file:      " + output_map_file,
            "minimum total counts: " + str(min_count),
            "minimum samples:      " + str(min_samples)])

    read_sample_names(sample_names_file)
    get_swarms(fasta_file, swarm_file, cpus)
    read_swarm_fasta(fasta_file)
    read_counts(counts_file)
    write_swarms(output_fasta_file, output_counts_file, output_map_file, min_samples, min_count)

if __name__ == "__main__":
    main(sys.argv)
