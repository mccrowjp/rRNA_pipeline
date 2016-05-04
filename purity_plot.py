#!/usr/bin/env python
#
# purity_plot - classification of reads in swarm OTUs to measure and plot OTU purity
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
import gzip, bz2

prog_path = os.path.realpath(sys.argv[0])
prog_dir = os.path.dirname(prog_path)
R_script = os.path.join(prog_dir, "plot_OTU_purity.r")

verbose = False

dict_id_swarm = {}
dict_id_counts = {}
dict_derep_ids = {}
dict_id_best_hit = {}
dict_id_best_bs = {}
dict_id_taxonomy = {}

def hopen(infile):
    f = None
    try:
        if re.search('\.bz2$', infile):
            f = bz2.BZ2File(infile, 'r')
        elif re.search('\.gz$', infile):
            f = gzip.GzipFile(infile, 'r')
        else:
            f = open(infile)
    except IOError:
        return None
    return f

def hopen_or_else(infile):
    f = hopen(infile)
    if f:
        return f
    else:
        print >>sys.stderr, "Unable to open file: " + infile

def hopen_write(outfile):
    f = None
    try:
        f = open(outfile, 'w')
    except IOError:
        return None
    return f

def hopen_write_or_else(outfile):
    f = hopen_write(outfile)
    if f:
        return f
    else:
        print >>sys.stderr, "Unable to write to file: " + outfile

def read_swarms(swarm_file):
    in_handle = hopen_or_else(swarm_file)
    if verbose:
        print >>sys.stderr, "Reading swarm file: " + swarm_file
    
    while 1:
        line = in_handle.readline()
        if not line:
            break
        line = line.rstrip()
        
        id_list = re.split('\s', line)
        for id in id_list:
            dict_id_swarm[id] = id_list[0]
    in_handle.close()

def read_swarm_counts(swarm_counts_file, min_swarm_count):
    dict_swarm_counts = {}
    global dict_derep_ids

    in_handle = hopen_or_else(swarm_counts_file)
    
    if verbose:
        print >>sys.stderr, "Reading swarm counts file: " + swarm_counts_file
        
    firstline = 1
    while 1:
        line = in_handle.readline()
        if not line:
            break
        line = line.rstrip()
        
        if not firstline:
            cols = line.split("\t")
            for i in range(1, len(cols)):
                dict_swarm_counts[cols[0]] = dict_swarm_counts.get(cols[0], 0) + int(cols[i])
        
        firstline = 0

    in_handle.close()

    for id in dict_id_swarm:
        swarm_id = dict_id_swarm[id]
        if dict_swarm_counts.get(swarm_id, 0) >= min_swarm_count:
            dict_derep_ids[id] = 1

def write_swarm_content(fasta_file, swarm_content_fasta_file):
    in_handle = hopen_or_else(fasta_file)
        
    if verbose:
        print >>sys.stderr, "Reading FASTA file: " + fasta_file
        
    out_handle = hopen_write_or_else(swarm_content_fasta_file)

    if verbose:
        print >>sys.stderr, "Writing swarm content FASTA file: " + swarm_content_fasta_file
        
    write_out = False
    while 1:
        line = in_handle.readline()
        if not line:
            break
        line = line.rstrip()
        
        if line.startswith(">"):
            id = re.split('\s', line[1:])[0]
            if id in dict_derep_ids:
                write_out = True
            else:
                write_out = False

        if write_out:
            print >>out_handle, line
    
    in_handle.close()
    out_handle.close()

def get_taxonomy(swarm_content_fasta_file, swarm_content_ggsearch_file, database_file, cpus):
    global dict_id_best_hit
    global dict_id_best_bs
    global dict_id_taxonomy

    if swarm_content_fasta_file:
        if os.path.exists(swarm_content_ggsearch_file):
            if verbose:
                print >>sys.stderr, "Ignoring content FASTA file: " + swarm_content_fasta_file
        else:
            print >>sys.stderr, "[purity_plot] running ggsearch"
            
            if cpus < 1:
                cpus = 1
            
            cmd = " ".join(["glsearch36 -b 1 -m 8 -T", str(cpus), swarm_content_fasta_file, database_file, ">", swarm_content_ggsearch_file])
            
            if verbose:
                print >>sys.stderr, cmd
            
            rc = os.system(cmd)
            if rc != 0:
                print >>sys.stderr, "[purity_plot] ERROR: ggsearch"
                sys.exit(2)

    in_handle1 = hopen_or_else(swarm_content_ggsearch_file)
    if verbose:
        print >>sys.stderr, "Reading ggsearch file: " + swarm_content_ggsearch_file
    
    while 1:
        line = in_handle1.readline()
        if not line:
            break
        line = line.rstrip()

        qid, sid, pid, alen, mm, go, qs, qe, ss, se, e, bs = line.split("\t")[:12]
        if (not qid in dict_id_best_bs) or bs > dict_id_best_bs[qid]:
            dict_id_best_hit[qid] = sid
            dict_id_best_bs[qid] = bs
    in_handle1.close()

    in_handle2 = hopen_or_else(database_file)
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

def write_purity(output_swarm_content_tax_file, output_swarm_purity_file, output_purity_pdf):
    if output_swarm_content_tax_file:
        out_handle1 = hopen_write_or_else(output_swarm_content_tax_file)
        if verbose:
            print >>sys.stderr, "Writing swarm content taxonomy file: " + output_swarm_content_tax_file
        
        print >>out_handle1, "\t".join(['id', 'swarm', 'besthit', 'taxonomy'])

        for id in dict_derep_ids:
            swarm_id = dict_id_swarm.get(id, "")
            besthit = dict_id_best_hit.get(id, "")
            tax = ""
            if besthit:
                tax = dict_id_taxonomy.get(besthit, "")

            print >>out_handle1, "\t".join([id, swarm_id, besthit, tax])

        out_handle1.close()

    out_handle2 = hopen_write_or_else(output_swarm_purity_file)

    if verbose:
        print >>sys.stderr, "Writing swarm purity file: " + output_swarm_purity_file

    count_all = {}
    count_same_tax = {}
    for id in dict_derep_ids:
        id_key, id_size = id.split('_')[:2]
        swarm_id = dict_id_swarm.get(id, "")
        
        derep_size = int(id_size)
        if derep_size < 1:
            derep_size = 1
        
        if swarm_id:
            count_all[swarm_id] = count_all.get(swarm_id, 0) + derep_size
            if id == swarm_id:
                count_same_tax[swarm_id] = count_same_tax.get(swarm_id, 0) + derep_size
            else:
                besthit = dict_id_best_hit.get(id, "")
                besthit_swarm = dict_id_best_hit.get(swarm_id, "")
                if besthit and besthit_swarm:
                    id_tax = dict_id_taxonomy.get(besthit, "")
                    swarm_tax = dict_id_taxonomy.get(besthit_swarm, "")
                    if id_tax and id_tax == swarm_tax:
                        count_same_tax[swarm_id] = count_same_tax.get(swarm_id, 0) + derep_size

    print >>out_handle2, "\t".join(['swarm_id', 'taxonomy', 'size', 'same_taxon', 'purity'])

    print >>sys.stderr, str(len(count_all)) + " " + str(len(count_same_tax))

    count_pure_OTUs = 0
    for swarm_id in count_all:
        print >>out_handle2, "\t".join([swarm_id, dict_id_taxonomy.get(swarm_id, ""), str(count_all[swarm_id]), str(count_same_tax.get(swarm_id, 0)), str(1.0 * count_same_tax.get(swarm_id, 0) / count_all[swarm_id])])
        if count_same_tax.get(swarm_id, 0) == count_all[swarm_id]:
            count_pure_OTUs += 1

    out_handle2.close()

    print >>sys.stderr, "OTUs 100% purity: " + str(count_pure_OTUs) + " / " + str(len(count_all)) + " (" + str(round(100.0 * count_pure_OTUs / len(count_all), 1)) + "%)"

    cmd = " ".join([R_script, output_swarm_purity_file, output_purity_pdf])

    if verbose:
        print >>sys.stderr, cmd
        
    rc = os.system(cmd + " >/dev/null")
    if rc != 0:
        print >>sys.stderr, "[purity_plot] ERROR: " + R_script
        sys.exit(2)

###

def main(argv):
    help = "\n".join([
        "purity_plot v0.2 (May 3, 2016)",
        "swarm OTU purity calculate and plot PDF",
        "",
        "Usage: "+argv[0]+" (options)",
        "   -f file        : dereplicated FASTA",
        "   -s file        : swarm file",
        "   -c file        : swarm counts file",
        "   -d file        : database FASTA file",
        "   -m int         : minimum swarm OTU count (default: 100)",
        "   -o file        : output base filename",
        "   -t, --cpus int : number of processes to run ggsearch (default: 1)",
        "   -h, --help     : help",
        "   -v, --verbose  : more information to stderr", ""])

    global verbose
    fasta_file = ""
    swarm_file = ""
    swarm_counts_file = ""
    base_file = ""
    swarm_content_fasta_file = ""
    swarm_content_ggsearch_file = ""
    database_file = ""
    min_swarm_count = 100
    output_swarm_content_tax_file = ""
    output_swarm_purity_file = ""
    output_purity_pdf = ""
    cpus = 1
    
    try:
        opts, args = getopt.getopt(argv[1:], "f:s:c:d:m:o:t:hv", ["cpus=", "help", "verbose", "test"])
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
        elif opt == '-c':
            swarm_counts_file = arg
        elif opt == '-d':
            database_file = arg
        elif opt == '-m':
            min_swarm_count = int(re.sub('=','', arg))
        elif opt == '-o':
            base_file = arg
        elif opt in ("-t", "--cpus"):
            cpus = int(re.sub('=','', arg))
        elif opt in ("-v", "--verbose"):
            verbose = True

    if not (fasta_file and swarm_file and swarm_counts_file and base_file and database_file):
        print >>sys.stderr, help
        sys.exit()

    swarm_content_ggsearch_file = base_file + ".content.ggsearch"
    swarm_content_fasta_file = base_file + ".content.fa"
    output_swarm_content_tax_file = base_file + ".content.tax"
    output_swarm_purity_file = base_file + ".purity"
    output_purity_pdf = base_file + ".purity.pdf"

    if verbose:
        print >>sys.stderr, "\n".join([
            "input fasta file:           " + fasta_file,
            "input swarm file:           " + swarm_file,
            "input database file:        " + database_file,
            "minimum swarm count:        " + str(min_swarm_count),
            "swarm content fasta file:   " + swarm_content_fasta_file,
            "ggsearch file:              " + swarm_content_ggsearch_file,
            "content taxonomy file:      " + output_swarm_content_tax_file,
            "output swarm content table: " + output_swarm_purity_file,
            "output purity pdf:          " + output_purity_pdf,
            "cpus:                       " + str(cpus)])

    read_swarms(swarm_file)
    read_swarm_counts(swarm_counts_file, min_swarm_count)
    write_swarm_content(fasta_file, swarm_content_fasta_file)
    get_taxonomy(swarm_content_fasta_file, swarm_content_ggsearch_file, database_file, cpus)
    write_purity(output_swarm_content_tax_file, output_swarm_purity_file, output_purity_pdf)

if __name__ == "__main__":
    main(sys.argv)
