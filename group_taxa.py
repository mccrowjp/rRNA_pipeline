#!/usr/bin/env python
#
# group_taxa - merge sample counts by taxonomic group
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

verbose = False

sample_list = []
dict_taxonomy_group = {}
dict_taxa_counts = {}
dict_taxa_sample_counts = {}
dict_group_counts = {}
dict_group_sample_counts = {}

def read_groups(taxa_group_file):
    global dict_taxonomy_group
    
    in_handle = happyfile.hopen_or_else(taxa_group_file)

    if verbose:
        print >>sys.stderr, "Reading taxa groups file: " + taxa_group_file
    
    while 1:
        line = in_handle.readline()
        if not line:
            break
        line = line.rstrip()
        
        if line:
            group_name, taxstr = line.split("\t")[:2]
            if taxstr and group_name:
                dict_taxonomy_group[taxstr] = group_name

    in_handle.close()

def read_taxa_counts(swarm_tax_file):
    global dict_taxa_counts
    global dict_taxa_sample_counts
    global dict_group_counts
    global dict_group_sample_counts
    global sample_list
    
    in_handle = happyfile.hopen_or_else(swarm_tax_file)
    
    if verbose:
        print >>sys.stderr, "Reading taxa counts file: " + swarm_tax_file
    
    firstline = 1
    while 1:
        line = in_handle.readline()
        if not line:
            break
        line = line.rstrip()
        
        if firstline:
            sample_list = line.split("\t")[3:]
        else:
            cols = line.split("\t")
            taxstr = cols[2]
            for i in range(3, len(cols)):
                dict_taxa_sample_counts[taxstr, i-3] = dict_taxa_sample_counts.get((taxstr, i-3), 0) + int(cols[i])
                dict_taxa_counts[taxstr] = dict_taxa_counts.get(taxstr, 0) + int(cols[i])
        
        firstline = 0
            
    in_handle.close()

    for id_tax in dict_taxa_counts:
        best_grp_tax = ""
        for grp_tax in dict_taxonomy_group:
            sub_tax = id_tax[:len(grp_tax)]
            if sub_tax == grp_tax:
                if (not best_grp_tax) or len(grp_tax) > len(best_grp_tax):
                    best_grp_tax = grp_tax
    
        best_grp_name = dict_taxonomy_group.get(best_grp_tax, "Unclassified")

        for i in range(len(sample_list)):
            dict_group_sample_counts[best_grp_name, i] = dict_group_sample_counts.get((best_grp_name, i), 0) + dict_taxa_sample_counts.get((id_tax, i), 0)
            dict_group_counts[best_grp_name] = dict_group_counts.get(best_grp_name, 0) + dict_taxa_counts.get(id_tax, 0)

def write_group_counts(output_groups_file):
    out_handle = sys.stdout
    if output_groups_file:
        out_handle = happyfile.hopen_write_or_else(output_groups_file)

    if verbose and output_groups_file:
        print >>sys.stderr, "Writing group counts file: " + output_groups_file

    column_names = ['group']
    for name in sample_list:
        column_names.append(name)

    print >>out_handle, "\t".join(column_names)

    for group_name in sorted(dict_group_counts, key=lambda x: dict_group_counts.get(x), reverse=True):
        samplecounts = [group_name]
        for i in range(len(sample_list)):
            samplecounts.append(dict_group_sample_counts.get((group_name, i), 0))
        print >>out_handle, "\t".join(str(x) for x in samplecounts)

    if output_groups_file:
        out_handle.close()

def test_all():
    print >>sys.stderr, "[group_taxa] test_all: passed"

###

def main(argv):
    help = "\n".join([
        "group_taxa v0.4 (May 21, 2016)",
        "merge OTU counts by taxonomic groups",
        "",
        "Usage: " + os.path.basename(argv[0]) + " (options)",
        "   -f file        : swarm taxonomy file",
        "   -g file        : taxonomic groups file",
        "   -o file        : output group counts file (default: stdout)",
        "   -h, --help     : help",
        "   -v, --verbose  : more information to stderr", ""])

    global verbose
    taxa_group_file = ""
    swarm_tax_file = ""
    output_groups_file = ""
    
    try:
        opts, args = getopt.getopt(argv[1:], "f:g:o:hv", ["help", "verbose", "test"])
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
            swarm_tax_file = arg
        elif opt == '-g':
            taxa_group_file = arg
        elif opt == '-o':
            output_groups_file = arg
        elif opt in ("-v", "--verbose"):
            verbose = True

    if not (swarm_tax_file and taxa_group_file):
        print >>sys.stderr, help
        sys.exit(2)

    if verbose:
        print >>sys.stderr, "\n".join([
            "input taxa groups file: " + taxa_group_file,
            "input taxa counts file: " + swarm_tax_file,
            "output groups file:     " + output_groups_file])

    read_groups(taxa_group_file)
    read_taxa_counts(swarm_tax_file)
    write_group_counts(output_groups_file)

if __name__ == "__main__":
    main(sys.argv)
