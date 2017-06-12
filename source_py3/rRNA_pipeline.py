#!/usr/bin/env python
#
# rRNA_pipeline - FASTQ filtering, and swarm OTU classification of 16/18S barcodes
#
# Version: 0.4 (5/21/2016)
#
# Original version: 4/11/2016 John P. McCrow (jmccrow [at] jcvi.org)
# J. Craig Venter Institute (JCVI)
# La Jolla, CA USA
#

import sys, re, getopt
import os, shutil
import happyfile

prog_path = os.path.realpath(sys.argv[0])
prog_dir = os.path.dirname(prog_path)
db_dir = os.path.join(prog_dir, 'db')
taxa_groups_file = os.path.join(prog_dir, 'taxa_groups.txt')

dict_database_path = {'16S' : os.path.join(db_dir,'db_16S.fa'), '18S' : os.path.join(db_dir,'db_18S.fa'), 'V4' : os.path.join(db_dir,'db_V4.fa'), 'V9' : os.path.join(db_dir,'db_V9.fa'), 'chloro' : os.path.join(db_dir,'db_plastid.fa')}

list_seq_file_pairs = []
cpus = 1
verbose = False
overwrite = False
do_chimera_search = True

def xstr(s):
    if s is None:
        return ''
    return str(s)

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
            m = re.match('(.+)_R[12].*\.f\w+$', b)
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

# For compatability with different operating systems
def replace_file(src, dst):
    if os.path.exists(src):
        if os.path.exists(dst):
            os.remove(dst)
        shutil.move(src, dst)

def get_seq_file_pairs(fastq_dir):
    global list_seq_file_pairs
    
    for file in os.listdir(fastq_dir):
        f1 = os.path.join(fastq_dir, file)
        m1 = re.search('^(.+)_R([12])(.*)\.(fastq|fq)$', file)
        if m1:
            if m1.group(2) == '1':
                f2 = os.path.join(fastq_dir, m1.group(1) + "_R2" + xstr(m1.group(3)) + "." + m1.group(4))
                if os.path.exists(f2):
                    list_seq_file_pairs.append(SequenceFilePair(f1, f2, True))
                else:
                    list_seq_file_pairs.append(SequenceFilePair(f1, '', False))

        else:
            m2 = re.search('^(.+)\.(fastq|fq)$', file)
            if m2:
                list_seq_file_pairs.append(SequenceFilePair(f1, '', False))

def is_fasta(fasta_file):
    in_handle = happyfile.hopen(fasta_file)
    retval = False
    if in_handle:
        line = in_handle.readline()
        if line and re.match('^>', line):
            retval = True
        in_handle.close()
    return retval

def run_command(name, checkfile, cmd_exe, cmd_params, redirect_all):
    if overwrite or not os.path.exists(checkfile):
        print("[rRNA_pipeline] running " + name + " " + checkfile, file=sys.stderr)

        cmd = cmd_exe + " "
        if verbose and not redirect_all:
            cmd += " -v "
        cmd += cmd_params
        if not verbose and redirect_all:
            cmd += " &>/dev/null"

        if verbose:
            print(cmd, file=sys.stderr)

        rc = os.system(cmd)
        if rc != 0:
            print("[rRNA_pipeline] ERROR: " + name, file=sys.stderr)
            sys.exit(2)
    else:
        print("[rRNA_pipeline] skipping " + name + " " + checkfile, file=sys.stderr)

def run_merge_fastq(fp):
    if fp.ispaired:
        cmd_params = " ".join(["-q 25 -t 50 --threads", str(cpus), "-f", fp.fastq1, "-r", fp.fastq2, "-o", fp.basefile])
        
        run_command('pear', fp.pear, "pear", cmd_params, True)

def run_usearch(fp, database_file):
    global do_chimera_search
    cmd_params = " ".join(["-threads", str(cpus), "-uchime_ref", fp.pear, "-db", database_file, "-uchimeout", fp.chimera, "-strand plus"])

    if not (do_chimera_search and is_fasta(fp.pear)):
        # create empty file, so that step will be skipped, but reported
        open(fp.chimera, 'a').close()

    run_command('chimera', fp.chimera, "usearch", cmd_params, True)

def run_filter(fp, min_quality_score):
    cmd_params = " ".join(["-f", fp.pear, "-o", fp.filtered, "-c", fp.chimera, "-q", str(min_quality_score)])
    
    run_command('filter', fp.filtered, os.path.join(prog_dir, "fastq_filter.py"), cmd_params, False)

def run_dereplicate(output_base_file, sample_names_file):
    derep_fa = output_base_file + ".derep.fa"
    derep_counts = output_base_file + ".derep.counts"
    cmd_params = "-o " + derep_fa + " -c " + derep_counts + " -t 3"
    if len(list_seq_file_pairs) > 1:
        cmd_params += " -l 2"
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
    
    cmd_params = " ".join(["-t", str(cpus), "-f", swarm_fa, "-g", ggsearch_file, "-d", database_file, "-c", swarm_counts, "-o", outfile])
    
    run_command('classify', outfile, os.path.join(prog_dir, "swarm_classify_taxonomy.py"), cmd_params, False)

def remove_plastid_seqs(output_base_file):
    dict_plastid = {}
    swarm_tax = output_base_file + ".swarm.tax"

    derep_fa = output_base_file + ".derep.fa"
    derep_counts = output_base_file + ".derep.counts"
    swarm_table = output_base_file + ".swarm"
    swarm_fa = output_base_file + ".swarm.fa"
    swarm_counts = output_base_file + ".swarm.counts"
    
    derep_plastid_fa = output_base_file + ".plastid.derep.fa"
    derep_plastid_counts = output_base_file + ".plastid.derep.counts"
    swarm_plastid_table = output_base_file + ".plastid.swarm"
    swarm_plastid_fa = output_base_file + ".plastid.swarm.fa"
    swarm_plastid_counts = output_base_file + ".plastid.swarm.counts"

    tmp_derep_16S_fa = derep_fa + ".tmp"
    tmp_derep_16S_counts = derep_counts + ".tmp"
    tmp_swarm_16S_table = swarm_table + ".tmp"
    tmp_swarm_16S_fa = swarm_fa + ".tmp"
    tmp_swarm_16S_counts = swarm_counts + ".tmp"
    tmp_swarm_16S_tax = swarm_tax + ".tmp"

    if overwrite or not os.path.exists(swarm_plastid_fa):
        if verbose:
            print("Filtering chloroplast sequences", file=sys.stderr)

        # split 16S/Plastid swarm taxonomy table
        in_handle_tax = happyfile.hopen_or_else(swarm_tax)
        out_handle_16S_tax = happyfile.hopen_write_or_else(tmp_swarm_16S_tax)

        firstline = 1
        while 1:
            line = in_handle_tax.readline()
            if not line:
                break
            line = line.rstrip()
            cols = line.split('\t')

            if firstline:
                print(line, file=out_handle_16S_tax)
            else:
                m = re.match('Bacteria;Cyanobacteria;Chloroplast', cols[2])
                if m:
                    dict_plastid[cols[0]] = 1
                else:
                    print(line, file=out_handle_16S_tax)

            firstline = 0

        in_handle_tax.close()
        out_handle_16S_tax.close()

        # split 16S/Plastid swarm file
        in_handle_table = happyfile.hopen_or_else(swarm_table)
        out_handle_16S_table = happyfile.hopen_write_or_else(tmp_swarm_16S_table)
        out_handle_plastid_table = happyfile.hopen_write_or_else(swarm_plastid_table)
        
        while 1:
            line = in_handle_table.readline()
            if not line:
                break
            line = line.rstrip()
            id_list = re.split('\s', line)
            swarm_id = id_list[0]
            
            if swarm_id in dict_plastid:
                print(line, file=out_handle_plastid_table)
                for id in id_list:
                    dict_plastid[id] = 1
            else:
                print(line, file=out_handle_16S_table)

        in_handle_table.close()
        out_handle_16S_table.close()
        out_handle_plastid_table.close()

        # split 16S/Plastid derep FASTA
        in_handle_derep_fa = happyfile.hopen_or_else(derep_fa)
        out_handle_16S_derep_fa = happyfile.hopen_write_or_else(tmp_derep_16S_fa)
        out_handle_plastid_derep_fa = happyfile.hopen_write_or_else(derep_plastid_fa)
        
        id = ""
        while 1:
            line = in_handle_derep_fa.readline()
            if not line:
                break
            line = line.rstrip()
            
            if line.startswith(">"):
                id = re.split('\s', line[1:])[0]
            
            if id:
                if id in dict_plastid:
                    print(line, file=out_handle_plastid_derep_fa)
                else:
                    print(line, file=out_handle_16S_derep_fa)

        in_handle_derep_fa.close()
        out_handle_16S_derep_fa.close()
        out_handle_plastid_derep_fa.close()
        
        # split 16S/Plastid derep counts table
        in_handle_derep_counts = happyfile.hopen_or_else(derep_counts)
        out_handle_16S_derep_counts = happyfile.hopen_write_or_else(tmp_derep_16S_counts)
        out_handle_plastid_derep_counts = happyfile.hopen_write_or_else(derep_plastid_counts)
        
        firstline = 1
        while 1:
            line = in_handle_derep_counts.readline()
            if not line:
                break
            line = line.rstrip()
            cols = line.split('\t')
            
            if firstline:
                print(line, file=out_handle_16S_derep_counts)
                print(line, file=out_handle_plastid_derep_counts)
            else:
                if cols[0] in dict_plastid:
                    print(line, file=out_handle_plastid_derep_counts)
                else:
                    print(line, file=out_handle_16S_derep_counts)

            firstline = 0

        in_handle_derep_counts.close()
        out_handle_16S_derep_counts.close()
        out_handle_plastid_derep_counts.close()

        # split 16S/Plastid swarm FASTA
        in_handle_fa = happyfile.hopen_or_else(swarm_fa)
        out_handle_16S_fa = happyfile.hopen_write_or_else(tmp_swarm_16S_fa)
        out_handle_plastid_fa = happyfile.hopen_write_or_else(swarm_plastid_fa)

        id = ""
        while 1:
            line = in_handle_fa.readline()
            if not line:
                break
            line = line.rstrip()
            
            if line.startswith(">"):
                id = re.split('\s', line[1:])[0]

            if id:
                if id in dict_plastid:
                    print(line, file=out_handle_plastid_fa)
                else:
                    print(line, file=out_handle_16S_fa)

        in_handle_fa.close()
        out_handle_16S_fa.close()
        out_handle_plastid_fa.close()

        # split 16S/Plastid swarm counts table
        in_handle_counts = happyfile.hopen_or_else(swarm_counts)
        out_handle_16S_counts = happyfile.hopen_write_or_else(tmp_swarm_16S_counts)
        out_handle_plastid_counts = happyfile.hopen_write_or_else(swarm_plastid_counts)

        firstline = 1
        while 1:
            line = in_handle_counts.readline()
            if not line:
                break
            line = line.rstrip()
            cols = line.split('\t')

            if firstline:
                print(line, file=out_handle_16S_counts)
                print(line, file=out_handle_plastid_counts)
            else:
                if cols[0] in dict_plastid:
                    print(line, file=out_handle_plastid_counts)
                else:
                    print(line, file=out_handle_16S_counts)

            firstline = 0
        
        in_handle_counts.close()
        out_handle_16S_counts.close()
        out_handle_plastid_counts.close()

        # replace original swarm files with 16S only
        if os.path.exists(tmp_derep_16S_fa) and os.path.exists(tmp_derep_16S_counts) and os.path.exists(tmp_swarm_16S_table) and os.path.exists(tmp_swarm_16S_tax) and os.path.exists(tmp_swarm_16S_fa) and os.path.exists(tmp_swarm_16S_counts):
            replace_file(tmp_derep_16S_fa, derep_fa)
            replace_file(tmp_derep_16S_counts, derep_counts)
            replace_file(tmp_swarm_16S_table, swarm_table)
            replace_file(tmp_swarm_16S_tax, swarm_tax)
            replace_file(tmp_swarm_16S_fa, swarm_fa)
            replace_file(tmp_swarm_16S_counts, swarm_counts)
        else:
            print("Not all tmp_ files found", file=sys.stderr)
            sys.exit(2)

def run_classify_chloro(output_base_file, database_file):
    swarm_counts = output_base_file + ".plastid.swarm.counts"
    swarm_plastid_fa = output_base_file + ".plastid.swarm.fa"
    ggsearch_file = output_base_file + ".plastid.swarm.ggsearch"
    swarm_plastid_tax = output_base_file + ".plastid.swarm.tax"
    
    cmd_params = " ".join(["-t", str(cpus), "-f", swarm_plastid_fa, "-g", ggsearch_file, "-d", database_file, "-c", swarm_counts, "-o", swarm_plastid_tax])
    
    run_command('classify_plastid', swarm_plastid_tax, os.path.join(prog_dir, "swarm_classify_taxonomy.py"), cmd_params, False)

def run_purity(output_base_file, database_file):
    derep_fa = output_base_file + ".derep.fa"
    swarm_file = output_base_file + ".swarm"
    swarm_counts = output_base_file + ".swarm.counts"

    cmd_params = " ".join(["-t", str(cpus), "-f", derep_fa, "-s", swarm_file, "-c", swarm_counts, "-d", database_file, "-o", output_base_file + ".swarm"])

    run_command('purity', output_base_file + ".swarm.purity", os.path.join(prog_dir, "purity_plot.py"), cmd_params, False)

def run_plot_sample_correlations(output_base_file, title):
    infile = output_base_file + ".swarm.tax"
    outfile = output_base_file + ".swarm.sample_corr.pdf"

    cmd_params = " ".join([infile, outfile, title])

    run_command('plot_sample_correlations', outfile, os.path.join(prog_dir, "plot_sample_correlations.r"), cmd_params, True)

def run_plot_taxa_groups(output_base_file):
    swarm_tax_file = output_base_file + ".swarm.tax"
    group_counts_file = output_base_file + ".taxa_groups.txt"
    plot_file = output_base_file + ".taxa_groups.pdf"

    cmd_params1 = " ".join(["-f", swarm_tax_file, "-g", taxa_groups_file, "-o", group_counts_file])
    run_command('group_taxa', group_counts_file, os.path.join(prog_dir, "group_taxa.py"), cmd_params1, False)

    cmd_params2 = group_counts_file + " " + plot_file
    run_command('plot_taxa_groups', plot_file, os.path.join(prog_dir, "plot_taxa_groups.r"), cmd_params2, True)

def run_plot_diversity(output_base_file):
    infile = output_base_file + ".swarm.tax"
    outfile = output_base_file + ".swarm.diversity.pdf"
    
    cmd_params = infile + " " + outfile
    
    run_command('plot_diversity', outfile, os.path.join(prog_dir, "plot_diversity.r"), cmd_params, True)

def run_plot_heatmap(output_base_file):
    infile = output_base_file + ".swarm.tax"
    outfile = output_base_file + ".swarm.heatmap.pdf"
    
    cmd_params = infile + " " + outfile
    
    run_command('plot_heatmap', outfile, os.path.join(prog_dir, "plot_heatmap.r"), cmd_params, True)

def run_plots(output_base_file, database_name):
    dict_title = {'16S' : '16S', '18S' : '18S', 'V4' : '18S_V4', 'V9' : '18S_V9', 'chloro' : 'Plastid'}
    title = dict_title.get(database_name, "")

    run_plot_sample_correlations(output_base_file, title)
    run_plot_taxa_groups(output_base_file)
    run_plot_diversity(output_base_file)
    run_plot_heatmap(output_base_file)

def init():
    global dict_database_path
    global taxa_groups_file
    global do_chimera_search
    init_file = os.path.join(prog_dir, 'init.txt')

    in_handle = happyfile.hopen(init_file)
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
                if db in ('16S', '18S', 'V4', 'V9', 'chloro'):
                    path = value
                    if not os.path.exists(path):
                        if os.path.exists(os.path.join(db_dir, path)):
                            path = os.path.join(db_dir, path)
                        elif os.path.exists(os.path.join(prog_dir, path)):
                            path = os.path.join(prog_dir, path)
                    if os.path.exists(path):
                        dict_database_path[db] = value
                if key == 'taxa_groups':
                    if os.path.exists(value):
                        taxa_groups_file = value
                    else:
                        taxa_groups_file = os.path.join(prog_dir, value)
                if key == 'chimera':
                    if re.match('^(off|no)', value.lower()):
                        do_chimera_search = False
    
        in_handle.close()

def test_each_dependency(cmd, name):
    failed = 0
    try:
        rc = os.system("which " + cmd + " &>/dev/null")
    except Exception:
        rc = 1
    if rc:
        failed = 1
        print("[rRNA_pipeline] test_dependencies: " + name + " failed", file=sys.stderr)
    else:
        print("[rRNA_pipeline] test_dependencies: " + name + " passed", file=sys.stderr)
    return failed

def test_dependencies():
    failed = 0
    failed += test_each_dependency("pear", "PEAR")
    failed += test_each_dependency("usearch", "USEARCH")
    failed += test_each_dependency("swarm", "SWARM")
    failed += test_each_dependency("glsearch36", "FASTA36")
    if failed:
        print("[rRNA_pipeline] test_dependencies: " + str(failed) + " test(s) failed", file=sys.stderr)
    else:
        print("[rRNA_pipeline] test_dependencies: All tests passed", file=sys.stderr)
    return failed

def test_databases():
    failed = 0
    try:
        init()
    except Exception:
        failed = 1
    for db in ('16S', '18S', 'V4', 'V9', 'chloro'):
        if os.path.exists(dict_database_path.get(db, "")):
            print("[rRNA_pipeline] test_databases: " + db + " found", file=sys.stderr)
        else:
            print("[rRNA_pipeline] test_databases: " + db + " database not found", file=sys.stderr)
            failed += 1
    if failed:
        print("[rRNA_pipeline] test_databases: " + str(failed) + " test(s) failed", file=sys.stderr)
    else:
        print("[rRNA_pipeline] test_databases: All tests passed", file=sys.stderr)
    return failed

def test_each_script(script):
    failed = 0
    full_script = os.path.join(prog_dir, script)
    try:
        rc = os.system(full_script + " --test")
    except Exception:
        rc = 1
    if rc:
        failed = 1
    return failed

def test_scripts():
    failed = 0
    failed += test_each_script("fastq_filter.py")
    failed += test_each_script("fasta_dereplicate.py")
    failed += test_each_script("swarm_map.py")
    failed += test_each_script("swarm_classify_taxonomy.py")
    failed += test_each_script("purity_plot.py")
    failed += test_each_script("group_taxa.py")
    if failed:
        print("[rRNA_pipeline] test_scripts: " + str(failed) + " test(s) failed", file=sys.stderr)
    else:
        print("[rRNA_pipeline] test_scripts: All tests passed", file=sys.stderr)
    return failed

def test_all():
    failed = 0
    failed += test_dependencies()
    failed += test_databases()
    failed += test_scripts()
    if failed:
        print("[rRNA_pipeline] test_all: " + str(failed) + " test(s) failed", file=sys.stderr)
        sys.exit(2)
    else:
        print("[rRNA_pipeline] test_all: All tests passed", file=sys.stderr)

###

def main(argv):
    help = "\n".join([
        "rRNA_pipeline v0.4 (May 21, 2016)",
        "Full ssu-rRNA, swarm OTU classification pipeline",
        "",
        "Usage: " + os.path.basename(prog_path) + " (options)",
        "   -d name          : database name (16S, 18S, V4, V9)",
        "   -q dir           : FASTQ folder",
        "   -o file          : base filename for results (default: rrna)",
        "   -n file          : sample names file (optional)",
        "   -p               : calculate/plot OTU purity",
        "   -m int           : minimum quality score for FASTQ (default: 30)",
        "   -s, --steps list : run only the steps in list (default: All)",
        "   -t, --cpus int   : number of processes (default: 1)",
        "   -w               : no overwrite of files, skip completed steps (default)" ,
        "   -W, --overwrite  : overwrite files (default if -s)",
        "   -h, --help       : help",
        "   -v, --verbose    : more information to stderr",
        "",
        "   If -s is used, then -W overwrite is the default",
        "   Steps list can include any of the following: ",
        "      All, merge_fastq, chimera, filter_fasta, derep",
        "      swarm, classify, plots, purity, split_plastid",
        "      plastid_classify, plastid_plots, plastid_purity",
        "",
        "Example: "+ os.path.basename(prog_path) + " -d 16S -q ./fastq -o rrna", ""])

    global verbose
    global overwrite
    global cpus
    global do_chimera_search
    database_name = ""
    database_file = ""
    fastq_dir = ""
    sample_names_file = ""
    calc_purity = False
    output_base_file = "rrna"
    min_quality_score = 30
    run_all_steps = True
    request_no_overwrite = False
    run_steps_str = ""
    dict_steps = {}
    
    try:
        opts, args = getopt.getopt(argv[1:], "d:q:o:n:pm:s:t:Wwhv", ["steps", "overwrite", "cpus=", "help", "verbose", "test"])
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
        elif opt == '-d':
            database_name = arg
        elif opt == '-q':
            fastq_dir = arg
        elif opt == '-o':
            output_base_file = arg
        elif opt == '-n':
            sample_names_file = arg
        elif opt == '-p':
            calc_purity = True
        elif opt == '-m':
            min_quality_score = int(re.sub('=','', arg))
        elif opt in ("-s", "--steps"):
            run_steps_str = arg
            run_all_steps = False
            if not request_no_overwrite:
                overwrite = True
        elif opt in ("-t", "--cpus"):
            cpus = int(re.sub('=','', arg))
        elif opt in ("-W", "--overwrite"):
            overwrite = True
        elif opt == '-w':
            overwrite = False
            request_no_overwrite = True
        elif opt in ("-v", "--verbose"):
            verbose = True

    if not (database_name):
        print(help, file=sys.stderr)
        sys.exit(2)

    output_base_file = output_base_file.rstrip('.')

    for s in re.split('[^\w_]+', run_steps_str.lower()):
        dict_steps[s] = 1
        if s in ('purity', 'plastid_purity'):
            calc_purity = True
        if s == 'all':
            run_all_steps = True

    derep_counts = output_base_file + ".derep.counts"
    if not fastq_dir and not os.path.exists(derep_counts):
        print(help + "\nFASTQ folder -f required if not found: " + derep_counts, file=sys.stderr)
        sys.exit(2)

    init()

    if database_name in ('16S', '18S', 'V4', 'V9'):
        database_file = dict_database_path[database_name]
    else:
        print(help + "\nDatabase name must be one of: 16S, 18S, V4, V9", file=sys.stderr)
        sys.exit(2)

    if verbose:
        print("\n".join([
            "input fastq dir:    " + fastq_dir,
            "input sample names: " + sample_names_file,
            "database name:      " + database_name,
            "database file:      " + database_file,
            "output base file:   " + output_base_file,
            "overwrite files:    " + ("no", "yes")[overwrite],
            "chimera search:     " + ("no", "yes")[do_chimera_search],
            "min fastq quality:  " + str(min_quality_score),
            "cpus:               " + str(cpus)]), file=sys.stderr)

    if cpus < 1:
        cpus = 1

    if fastq_dir:
        get_seq_file_pairs(fastq_dir)

        print("Found " + str(len(list_seq_file_pairs)) + " samples", file=sys.stderr)
        for fp in sorted(list_seq_file_pairs, key=lambda fp: fp.basefile):
            print(fp.basefile + " " + ("", " [paired]")[fp.ispaired], file=sys.stderr)

        for fp in sorted(list_seq_file_pairs, key=lambda fp: fp.basefile):
            if run_all_steps or 'merge_fastq' in dict_steps:
                run_merge_fastq(fp)
            if run_all_steps or 'chimera' in dict_steps:
                run_usearch(fp, database_file)
            if run_all_steps or 'filter_fasta' in dict_steps:
                run_filter(fp, min_quality_score)
    else:
        print("[rRNA_pipeline] skipping FASTQ merge/chimera/filtering", file=sys.stderr)
    
    if run_all_steps or 'derep' in dict_steps:
        run_dereplicate(output_base_file, sample_names_file)

    if run_all_steps or 'swarm' in dict_steps:
        run_swarm(output_base_file)

    if run_all_steps or 'classify' in dict_steps:
        run_classify(output_base_file, database_name)

    # 16S contains plastid sequences that need to be classified seperately against PhytoRef
    if database_name == '16S':
        if run_all_steps or 'split_plastid' in dict_steps:
            remove_plastid_seqs(output_base_file)
        
        if run_all_steps or 'plastid_classify' in dict_steps:
            run_classify_chloro(output_base_file, dict_database_path['chloro'])

        if run_all_steps or 'plastid_plots' in dict_steps:
            run_plots(output_base_file + ".plastid", 'chloro')

        if calc_purity and (run_all_steps or 'plastid_purity' in dict_steps):
            run_purity(output_base_file + ".plastid", dict_database_path['chloro'])

    if run_all_steps or 'plots' in dict_steps:
        run_plots(output_base_file, database_name)
    
    if calc_purity and (run_all_steps or 'purity' in dict_steps):
        run_purity(output_base_file, database_file)


if __name__ == "__main__":
    main(sys.argv)
