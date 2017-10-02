#!/usr/bin/perl
use strict ;
#use MJL_dataUtils ;

(@ARGV == 1) || die "usage: $0 fastaFile\n" ;
my $inFile = shift ;

my $outFile = $inFile ;
$outFile =~ s/\.\w+$// ;
$outFile .= 'NR0.fa' ;

my @seqs = readFasta ($inFile) ;
my $n = scalar @seqs ;

# restrict to non-redundant sequences:
my @seqsNR = nonRedundantSeqS (@seqs) ;
my $nNR = scalar @seqsNR ;

print "compacting $n sequences to $nNR non-redundant sequences...\n" ;

# write to file:
open (OUT, ">$outFile") || die ;
foreach my $seq (@seqsNR) {
	print OUT ">$seq->[0]   $seq->[2]\n$seq->[1]\n" ;
}
close OUT ;


sub nonRedundantSeqS (@) {
	my @seqs = @_ ;

	my $counter = my $percentage = 0 ;
	for (my $i1=0; $i1<@seqs-1; $i1++) {
		my $seq1 = substr $seqs[$i1][1], 1;
		#my $seq1 = $seqs[$i1][1];
		for (my $i2=$i1+1; $i2<@seqs; $i2++) {
			#my $seq2 =$seqs[$i2][1];
			my $seq2 = substr $seqs[$i2][1], 1;
			if (index ($seq1, $seq2) >= 0) {
				$seqs[$i2][2] = 'redundant' ;
			} elsif (index ($seq2, $seq1) >= 0) {
				$seqs[$i1][2] = 'redundant' ;
			}
		}
		if (int($counter++ / scalar (@seqs) * 100) > $percentage) {
			printf STDERR "\r%3d%%", $percentage ;
			$percentage = int($counter++ / scalar (@seqs) * 100) ;
		}
	}
	print STDERR "\n" ;
	# remove redundant proteins and write non-redundant set to file:
	return grep ($_->[2] ne 'redundant', @seqs) ;
}



sub readFasta {
	my $file = shift ;
	my $verbose = shift ;

	unless (open (IN, $file)) {
		warn "couldn't open Fasta file $file for reading" ;
		return () ;
	}
	
	my $header = '' ;
	while ($header = <IN>) {
		last if $header =~ /^>/ ;
	}

	my $seq = '' ;
	my @sequences = () ;
	my $counter = 0 ;
	while (1) {
		my $line = <IN> ;
		unless (($line =~ /^>/) || !$line) {
			$seq .= $line ;
			next ;
		}
		if (($line =~ /^>/) || eof IN) {
			my ($name, $desc) = ($header =~ /^>(\S*)\s+(.*)/) ;
			$seq =~ s/\s//g ;
			#$seq = substr $seq, 1; ;
			my @sequence = ($name, $seq, $desc) ;
			push @sequences, \@sequence ;
			$verbose && !($counter++ % 100) && printf STDERR "\r%6d seqs. read", $counter ;

			last if eof IN ;
			$header = $line ;
			$seq = '' ;
		}
	}
	close IN ;
	$verbose && printf (STDERR "\r%6d seqs. read", $counter) ;

	return @sequences ;
}
