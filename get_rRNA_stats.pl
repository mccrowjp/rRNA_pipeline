#!/usr/bin/env perl
use strict;

my $rrnadir = shift;
my $fastqdir = shift;
my $otufile = shift;
my $sampfile = shift;

unless($rrnadir && $fastqdir && $otufile) {
    die "Usage: $0 [rRNA_pipeline output directory] [FASTQ directory] [OTU count table .swarm.counts] ([sample names file])\n";
}

my %sampname;
my %samplib;
my %alllibs;
my %libseqs;
my %libbp;
my %libfiltseqs;
my %libfiltbp;
my %libassigned;
my %libotus;

if($sampfile) {
    open(IN, $sampfile) or die "Unable to open file $sampfile\n";
    while(<IN>) {
	chomp;
	my ($new, $old) = split(/\t/);
	$sampname{$old} = $new;
	$samplib{$new} = $old;
    }
    close(IN);
}

print STDERR "reading rRNA_pipeline output folder: $rrnadir\n";
foreach my $file (glob($rrnadir."/*")) {
    if($file =~ /\/([^\/_]+)_[^\/]+\.filtered\.fa/) {
	my $lib = $1;
	if(exists($sampname{$lib})) {
	    $lib = $sampname{$lib};
	}
	$alllibs{$lib} = 1;
	
	if($file =~ /\.gz$/) {
	    open(IN, "zcat $file 2>/dev/null |");
	} else {
	    open(IN, $file) or die "Unable to open file $file\n";
	}
	while(<IN>) {
	    chomp;
	    if(/^>/) {
		$libfiltseqs{$lib}++;
	    } else {
		s/[^a-zA-Z]//g;
		$libfiltbp{$lib} += length($_);
	    }
	}
	close(IN);
    }
}

print STDERR "reading FASTQ folder: $fastqdir\n";
foreach my $file (glob($fastqdir."/*")) {
    if($file =~ /\/([^\/_]+)_[^\/]+\.fastq/) {
	my $lib = $1;
	if(exists($sampname{$lib})) {
	    $lib = $sampname{$lib};
	}
	$alllibs{$lib} = 1;
	
	if($file =~ /\.gz$/) {
	    open(IN, "zcat $file 2>/dev/null |");
	} else {
	    open(IN, $file) or die "Unable to open file $file\n";
	}
	my $sn = 0;
	while(<IN>) {
	    chomp;
	    $sn++;
	    if($sn == 2) {
		$libseqs{$lib}++;
		$libbp{$lib} += length($_);
	    }
	    if($sn >= 4) {
		$sn = 0;
	    }
	}
	close(IN);
    }
}

open(IN, $otufile) or die "Unable to open file $otufile\n";
print STDERR "reading $otufile\n";

my $fl = 1;
my @head;
while(<IN>) {
    chomp;
    my ($id, @counts) = split(/\t/);
    if($fl) {
	for(my $i=0; $i<scalar(@counts); $i++) {
	    my $lib = $counts[$i];
	    if(exists($sampname{$lib})) {
		$lib = $sampname{$lib};
	    }
	    push(@head, $lib);
	    $alllibs{$lib} = 1;
	}
	
    } else {
	for(my $i=0; $i<scalar(@counts); $i++) {
	    $libassigned{$head[$i]} += $counts[$i];
	    if($counts[$i] > 0) {
		$libotus{$head[$i]}++;
	    }
	}
    }
    $fl = 0;
}
close(IN);

print join("\t", ('sample', 'library', 'total_sequences', 'total_bp', 'filtered_sequences', 'filtered_bp', 'otu_assigned_reads', 'otu_count'))."\n";
foreach my $lib (sort keys %alllibs) {
    my $samp = $lib;
    if(exists($samplib{$samp})) {
	$lib = $samplib{$samp};
    }
    printf join("\t", ("%s", "%s", "%d", "%d", "%d", "%d", "%d", "%d"))."\n", $samp, $lib, $libseqs{$samp}, $libbp{$samp}, $libfiltseqs{$samp}, $libfiltbp{$samp}, $libassigned{$samp}, $libotus{$samp};
}
