#!/usr/bin/env python
#
# rRNA_pipeline - FASTQ filtering, and swarm OTU classification of 16/18S barcodes
#
# Version: 0.2 (5/3/2016)
#
# Original version: 4/11/2016 John P. McCrow (jmccrow [at] jcvi.org)
# J. Craig Venter Institute (JCVI)
# La Jolla, CA USA
#

import sys, re, os, getopt
import gzip, bz2

prog_path = os.path.realpath(sys.argv[0])
prog_dir = os.path.dirname(prog_path)
db_dir = os.path.join(prog_dir, 'db')

dict_database_path = {'16S' : os.path.join(db_dir,'db_16S.fa'), 'V4' : os.path.join(db_dir,'db_V4.fa'), 'V9' : os.path.join(db_dir,'db_V9.fa'), 'chloro' : os.path.join(db_dir,'db_plastid.fa')}

list_seq_file_pairs = []
cpus = 1
verbose = False
overwrite = False

class SequenceFilePair:
    fastq1 = ""
    fastq2 = ""
    basefile = ""
    pear = ""
    chimera = ""
    filtered = ""
    ispaired = True

    def __init__(self, file1, file2, ispaired):
        b = os.path.basename(file1)
        if ispaired:
            m = re.match('(.+)_R[12](_filtered)*\.f\w+$', b)
            if m:
                b = m.group(1)
        else:
            m = re.match('(.+)\.f\w+$', b)
            if m:
                b = m.group(1)
        
        self.fastq1 = file1
        self.fastq2 = file2
        self.ispaired = ispaired
        self.basefile = b
        if ispaired:
            self.pear = self.basefile + ".assembled.fastq"
        else:
            self.pear = file1
        self.chimera = self.basefile + ".uchime"
        self.filtered = self.basefile + ".filtered.fa"

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

def get_seq_file_pairs(fastq_dir):
    global list_seq_file_pairs
    
    for file in os.listdir(fastq_dir):
        f1 = os.path.join(fastq_dir, file)
        m1 = re.search('^(.+)_R([12])(_filtered)*\.(fastq|fq)$', file)
        if m1:
            if m1.group(2) == '1':
                f2 = os.path.join(fastq_dir, m1.group(1) + "_R2" + m1.group(3) + "." + m1.group(4))
                if os.path.exists(f2):
                    list_seq_file_pairs.append(SequenceFilePair(f1, f2, True))
                else:
                    list_seq_file_pairs.append(SequenceFilePair(f1, '', False))
                    print >>sys.stderr, "Pair not found: " + f2
        else:
            m2 = re.search('^(.+)\.(fastq|fq)$', file)
            if m2:
                list_seq_file_pairs.append(SequenceFilePair(f1, '', False))

def run_command(name, checkfile, cmd_exe, cmd_params, redirect_all):
    if overwrite or not os.path.exists(checkfile):
        print >>sys.stderr, "[rRNA_ORFs] running " + name + " " + checkfile

        cmd = cmd_exe + " "
        if verbose and not redirect_all:
            cmd += " -v "
        cmd += cmd_params
        if not verbose and redirect_all:
            cmd += " &>/dev/null"

        if verbose:
            print >>sys.stderr, cmd

        rc = os.system(cmd)
        if rc != 0:
            print >>sys.stderr, "[rRNA_ORFs] ERROR: " + name
            sys.exit(2)
    else:
        print >>sys.stderr, "[rRNA_ORFs] skipping " + name + " " + checkfile

def run_merge_fastq(fp):
    if fp.ispaired:
        cmd_params = " ".join(["-q 35 -t 50 --threads", str(cpus), "-f", fp.fastq1, "-r", fp.fastq2, "-o", fp.basefile])
        
        run_command('pear', fp.pear, "pear", cmd_params, True)

def run_usearch(fp, database_file):
    cmd_params = " ".join(["-threads", str(cpus), "-uchime_ref", fp.pear, "-db", database_file, "-uchimeout", fp.chimera, "-strand plus"])
    
    run_command('chimera', fp.chimera, "usearch", cmd_params, True)

def run_filter(fp):
    cmd_params = " ".join(["-f", fp.pear, "-o", fp.filtered, "-c", fp.chimera])
    
    run_command('filter', fp.filtered, os.path.join(prog_dir, "fastq_filter.py"), cmd_params, False)

def run_dereplicate(output_base_file, sample_names_file):
    derep_fa = output_base_file + ".derep.fa"
    derep_counts = output_base_file + ".derep.counts"
    cmd_params = "-l 2 -t 3 -o " + derep_fa + " -c " + derep_counts
    if sample_names_file:
        cmd_params += " -n " + sample_names_file
    for fp in list_seq_file_pairs:
        cmd_params += " " + fp.filtered

    run_command('dereplicate', derep_fa, os.path.join(prog_dir, "fasta_dereplicate.py"), cmd_params, False)

def run_swarm(output_base_file):
    derep_fa = output_base_file + ".derep.fa"
    derep_counts = output_base_file + ".derep.counts"
    swarm_file = output_base_file + ".swarm"
    swarm_fa = output_base_file + ".swarm.fa"
    swarm_counts = output_base_file + ".swarm.counts"
    cmd_params = " ".join(["-x", str(cpus), "-f", derep_fa, "-d", derep_counts, "-s", swarm_file, "-o", swarm_fa, "-c", swarm_counts])
    
    run_command('swarm', swarm_fa, os.path.join(prog_dir, "swarm_map.py"), cmd_params, False)

def run_classify(output_base_file, database_name):
    database_file = dict_database_path[database_name]
    swarm_fa = output_base_file + ".swarm.fa"
    swarm_counts = output_base_file + ".swarm.counts"
    ggsearch_file = output_base_file + ".swarm.ggsearch"

    outfile = output_base_file + ".swarm.tax"
    if database_name == '16S':
        outfile = output_base_file + ".swarm.tmp.tax"
    
    cmd_params = " ".join(["-t", str(cpus), "-f", swarm_fa, "-g", ggsearch_file, "-d", database_file, "-c", swarm_counts, "-o", outfile])
    
    run_command('classifiy', outfile, os.path.join(prog_dir, "swarm_classify_taxonomy.py"), cmd_params, False)

def remove_plastid_seqs(swarm_fa, swarm_tmp, swarm_16S_tax, swarm_plastid_fa):
    dict_plastid = {}
    
    if verbose:
        print >>sys.stderr, "Filtering chloroplast sequences"
    
    in_handle1 = hopen_or_else(swarm_tmp)
    out_handle1 = hopen_write_or_else(swarm_16S_tax)
    
    firstline = 1
    while 1:
        line = in_handle1.readline()
        if not line:
            break
        line = line.rstrip()
        cols = line.split('\t')

        m = re.match('Bacteria;Cyanobacteria;Chloroplast', cols[2])
        if firstline or not m:
            print >>out_handle1, line
        else:
            dict_plastid[cols[0]] = 1
        firstline = 0

    in_handle1.close()
    out_handle1.close()

    in_handle2 = hopen_or_else(swarm_fa)
    out_handle2 = hopen_write_or_else(swarm_plastid_fa)

    id = ""
    while 1:
        line = in_handle2.readline()
        if not line:
            break
        line = line.rstrip()
        
        if line.startswith(">"):
            id = re.split('\s', line[1:])

        if id and id in dict_plastid:
            print >>out_handle2, line

    in_handle2.close()
    out_handle2.close()

def run_classify_chloro(output_base_file, database_file):
    swarm_fa = output_base_file + ".swarm.fa"
    swarm_tmp = output_base_file + ".swarm.tmp.tax"
    swarm_16S_tax = output_base_file + ".swarm.tax"
    swarm_plastid_fa = output_base_file + ".plastid.swarm.fa"
    remove_plastid_seqs(swarm_fa, swarm_tmp, swarm_16S_tax, swarm_plastid_fa)

    swarm_counts = output_base_file + ".swarm.counts"
    ggsearch_file = output_base_file + ".plastid.swarm.ggsearch"
    outfile = output_base_file + ".plastid.swarm.tax"

    cmd_params = " ".join(["-t", str(cpus), "-f", swarm_plastid_fa, "-g", ggsearch_file, "-d", database_file, "-c", swarm_counts, "-o", outfile])
    
    run_command('classifiy_plastid', outfile, os.path.join(prog_dir, "swarm_classify_taxonomy.py"), cmd_params, False)

def run_purity(output_base_file, database_file):
    derep_fa = output_base_file + ".derep.fa"
    swarm_file = output_base_file + ".swarm"
    swarm_counts = output_base_file + ".swarm.counts"

    cmd_params = " ".join(["-t", str(cpus), "-f", derep_fa, "-s", swarm_file, "-c", swarm_counts, "-d", database_file, "-o", output_base_file + ".swarm"])

    run_command('purity', output_base_file + ".swarm.purity", os.path.join(prog_dir, "purity_plot.py"), cmd_params, False)

def run_plots(output_base_file, database_name):
    dict_title = {'16S' : '16S', 'V4' : '18S_V4', 'V9' : '18S_V9', 'chloro' : 'Plastid'}
    title = dict_title.get(database_name, "")
    
    infile1 = output_base_file + ".swarm.tax"
    outfile1 = output_base_file + ".swarm.sample_corr.pdf"

    cmd_params = " ".join([infile1, outfile1, title])

    run_command('plot_sample_correlations', outfile1, os.path.join(prog_dir, "plot_sample_correlations.r"), cmd_params, True)


def init():
    global dict_database_path
    init_file = os.path.join(prog_dir, 'init.txt')

    in_handle = hopen(init_file)
    if in_handle:
        while 1:
            line = in_handle.readline()
            if not line:
                break
            line = line.rstrip()
            
            m = re.match('^(\S+):\s+(.+)$', line)
            if m:
                key, value = m.group(1), m.group(2)
                db = re.sub('^db_', '', key)
                if db in ('16S', 'V4', 'V9', 'chloro'):
                    path = value
                    if not os.path.exists(path):
                        if os.path.exists(os.path.join(db_dir, path)):
                            path = os.path.join(db_dir, path)
                        elif os.path.exists(os.path.join(prog_dir, path)):
                            path = os.path.join(prog_dir, path)
                    if os.path.exists(path):
                        dict_database_path[db] = value
        
        in_handle.close()


###

def main(argv):
    help = "\n".join([
        "rRNA_pipeline v0.2 (May 3, 2016)",
        "Full ssu-rRNA, swarm OTU classification pipeline",
        "",
        "Usage: "+argv[0]+" (options)",
        "   -d name         : database name (16S, V4, V9)",
        "   -q dir          : FASTQ folder",
        "   -o file         : base filename for results (default: rrna)",
        "   -n file         : sample names file (optional)",
        "   -t, --cpus int  : number of processes (default: 1)",
        "   -W, --overwrite : overwrite files (default: No, run next step)",
        "   -h, --help      : help",
        "   -v, --verbose   : more information to stderr", ""])

    global verbose
    global overwrite
    global cpus
    database_name = ""
    database_file = ""
    fastq_dir = ""
    sample_names_file = ""
    output_base_file = "rrna"
    
    try:
        opts, args = getopt.getopt(argv[1:], "d:q:o:n:t:Whv", ["overwrite", "cpus=", "help", "verbose", "test"])
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
        elif opt == '-d':
            database_name = arg
        elif opt == '-q':
            fastq_dir = arg
        elif opt == '-o':
            output_base_file = arg
        elif opt == '-n':
            sample_names_file = arg
        elif opt in ("-t", "--cpus"):
            cpus = int(re.sub('=','', arg))
        elif opt in ("-W", "--overwrite"):
            overwrite = True
        elif opt in ("-v", "--verbose"):
            verbose = True

    if not (database_name):
        print >>sys.stderr, help
        sys.exit(2)

    output_base_file = output_base_file.rstrip('.')

    derep_counts = output_base_file + ".derep.counts"
    if not fastq_dir and not os.path.exists(derep_counts):
        print >>sys.stderr, help + "\nFASTQ folder -f required if not found: " + derep_counts
        sys.exit(2)

    init()

    if database_name in ('16S', 'V4', 'V9'):
        database_file = dict_database_path[database_name]
    else:
        print >>sys.stderr, help + "\nDatabase name must be one of: 16S, V4, V9"
        sys.exit(2)

    if verbose:
        print >>sys.stderr, "\n".join([
            "input fastq dir:    " + fastq_dir,
            "input sample names: " + sample_names_file,
            "database name:      " + database_name,
            "database file:      " + database_file,
            "output base file:   " + output_base_file,
            "overwrite files:    " + ("no", "yes")[overwrite],
            "cpus:               " + str(cpus)])

    if cpus < 1:
        cpus = 1

    if fastq_dir:
        get_seq_file_pairs(fastq_dir)

        print >>sys.stderr, "Found " + str(len(list_seq_file_pairs)) + " samples"
        for fp in sorted(list_seq_file_pairs, key=lambda fp: fp.basefile):
            print >>sys.stderr, fp.basefile + " " + ("", " [paired]")[fp.ispaired]

        for fp in sorted(list_seq_file_pairs, key=lambda fp: fp.basefile):
            run_merge_fastq(fp)
            run_usearch(fp, database_file)
            run_filter(fp)
    else:
        print >>sys.stderr, "[rRNA_ORFs] skipping FASTQ merge/chimera/filtering"
    
    run_dereplicate(output_base_file, sample_names_file)
    run_swarm(output_base_file)
    run_classify(output_base_file, database_name)

    # 16S contains plastid sequences that need to be classified seperately against PhytoRef
    if database_name == '16S':
        run_classify_chloro(output_base_file, dict_database_path['chloro'])
        run_purity(output_base_file + ".plastid", dict_database_path['chloro'])
        run_plots(output_base_file + ".plastid")
    
    run_plots(output_base_file, database_name)
    run_purity(output_base_file, database_file)


if __name__ == "__main__":
    main(sys.argv)
