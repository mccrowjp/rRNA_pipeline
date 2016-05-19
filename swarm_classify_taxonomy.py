#!/usr/bin/env python
#
# swarm_classify_taxonomy - merge best reference taxonomy with sample counts for OTUs
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
dict_swarm_best_hit = {}
dict_swarm_best_bs = {}
dict_id_taxonomy = {}
dict_swarm_sample_counts = {}
dict_swarm_counts = {}

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

def read_counts(counts_file):
    global dict_swarm_counts
    global dict_swarm_sample_counts
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
                    dict_swarm_sample_counts[cols[0], i-1] = int(cols[i])
                    dict_swarm_counts[cols[0]] = dict_swarm_counts.get(cols[0], 0) + int(cols[i])
            
            firstline = 0
                
        in_handle.close()

def get_taxonomy(fasta_file, ggsearch_file, database_file, cpus):
    global dict_swarm_best_hit
    global dict_swarm_best_bs
    global dict_id_taxonomy

    if fasta_file:
        if os.path.exists(ggsearch_file):
            if verbose:
                print >>sys.stderr, "Ignoring FASTA file: " + fasta_file
        else:
            print >>sys.stderr, "[swarm_classify_taxonomy] running ggsearch"
            
            if cpus < 1:
                cpus = 1
            
            cmd = " ".join(["glsearch36 -b 1 -m 8 -T", str(cpus), fasta_file, database_file, ">", ggsearch_file])
            
            if verbose:
                print >>sys.stderr, cmd
            
            rc = os.system(cmd)
            if rc != 0:
                print >>sys.stderr, "[swarm_classify_taxonomy] ERROR: ggsearch"
                sys.exit(2)

    in_handle1 = happyfile.hopen_or_else(ggsearch_file)
    if verbose:
        print >>sys.stderr, "Reading ggsearch file: " + ggsearch_file
    
    while 1:
        line = in_handle1.readline()
        if not line:
            break
        line = line.rstrip()

        qid, sid, pid, alen, mm, go, qs, qe, ss, se, e, bs = line.split("\t")[:12]
        if (not qid in dict_swarm_best_bs) or bs > dict_swarm_best_bs[qid]:
            dict_swarm_best_hit[qid] = sid
            dict_swarm_best_bs[qid] = bs
    in_handle1.close()

    in_handle2 = happyfile.hopen_or_else(database_file)
    if verbose:
        print >>sys.stderr, "Reading database file: " + database_file
        
    while 1:
        line = in_handle2.readline()
        if not line:
            break
        line = line.rstrip()
        
        if line.startswith(">"):
            m = re.match('>(\S+)\s+(.+)$', line)
            if m:
                id = m.group(1)
                taxstr = m.group(2)
                taxstr = re.sub('\|', ';', taxstr)
                taxstr = re.sub('\+', ' ', taxstr)
                dict_id_taxonomy[id] = taxstr
    in_handle2.close()

def write_swarms(output_counts_file):
    out_handle = sys.stdout
    if output_counts_file:
        out_handle = happyfile.hopen_write_or_else(output_counts_file)

    if verbose and output_counts_file:
        print >>sys.stderr, "Writing counts file: " + output_counts_file

    column_names = ['id', 'besthit', 'taxonomy']
    for name in sample_list:
        if name in dict_sample_name:
            column_names.append(dict_sample_name[name])
        else:
            column_names.append(name)

    print >>out_handle, "\t".join(column_names)

    for swarm_id in dict_swarm_counts:
        besthit = dict_swarm_best_hit.get(swarm_id, "")
        tax = ""
        if besthit:
            tax = dict_id_taxonomy.get(besthit, "")
        samplecounts = [swarm_id, besthit, tax]
        for i in range(len(sample_list)):
            samplecounts.append(dict_swarm_sample_counts.get((swarm_id, i), 0))
        print >>out_handle, "\t".join(str(x) for x in samplecounts)

    if output_counts_file:
        out_handle.close()

def test_all():
    print >>sys.stderr, "[swarm_classify_taxonomy] test_all: passed"

###

def main(argv):
    help = "\n".join([
        "swarm_classify_taxonomy v0.2 (May 3, 2016)",
        "swarm OTU add taxonomy to counts table",
        "",
        "Usage: " + os.path.basename(argv[0]) + " (options)",
        "   -f file        : swarm FASTA",
        "   -g file        : ggsearch -m8 file",
        "   -d file        : database FASTA file",
        "   -c file        : swarm counts file",
        "   -o file        : output counts file (default: stdout)",
        "   -n file        : sample names file (optional)",
        "   -t, --cpus int : number of processes to run ggsearch (default: 1)",
        "   -h, --help     : help",
        "   -v, --verbose  : more information to stderr", ""])

    global verbose
    fasta_file = ""
    ggsearch_file = ""
    database_file = ""
    counts_file = ""
    sample_names_file = ""
    output_counts_file = ""
    cpus = 1
    
    try:
        opts, args = getopt.getopt(argv[1:], "f:g:d:c:o:n:t:hv", ["cpus=", "help", "verbose", "test"])
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
        elif opt == '-g':
            ggsearch_file = arg
        elif opt == '-d':
            database_file = arg
        elif opt == '-c':
            counts_file = arg
        elif opt == '-o':
            output_counts_file = arg
        elif opt == '-n':
            sample_names_file = arg
        elif opt in ("-t", "--cpus"):
            cpus = int(re.sub('=','', arg))
        elif opt in ("-v", "--verbose"):
            verbose = True

    if not (ggsearch_file and database_file and counts_file):
        print >>sys.stderr, help
        sys.exit(2)

    if verbose:
        print >>sys.stderr, "\n".join([
            "input fasta file:     " + fasta_file,
            "input ggsearch file:  " + ggsearch_file,
            "input counts file:    " + counts_file,
            "input database file:  " + database_file,
            "input sample names:   " + sample_names_file,
            "output counts file:   " + output_counts_file,
            "cpus:                 " + str(cpus)])

    read_sample_names(sample_names_file)
    get_taxonomy(fasta_file, ggsearch_file, database_file, cpus)
    read_counts(counts_file)
    write_swarms(output_counts_file)

if __name__ == "__main__":
    main(sys.argv)
