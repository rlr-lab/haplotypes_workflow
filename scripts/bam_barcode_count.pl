#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "
bam_barcode_count.pl

Requires samtools to be in your PATH

Required:
  -b    path to BAM file
  -s    barcode start position in alignment (1-based)
  -e    barcode end position in alignment (1-based)

Options:
  -n    Reference sequence name
        (default: will output separate results for each reference in the alignment)
  -l    Slop amount, in bp. Reads starting and or ending less than this distance upstream 
        or downstream of the barcode start or stop position will be ignored.
        (default: 5)

";

use Getopt::Std;
use vars qw( $opt_b $opt_s $opt_e $opt_n $opt_l $opt_d);
getopts('b:s:e:n:l:d:');

die $usage unless ($opt_b and $opt_s and $opt_e);

my $input   = $opt_b;
my $istart  = $opt_s;
my $istop   = $opt_e;
my $ileng = $istop - $istart + 1;

my $name    = $opt_n if $opt_n;
my $slop    = $opt_l ? $opt_l : 5;
my $debug   = $opt_d ? $opt_d : " ";

## Get reference names
my %rseen;
my @refs;
open (my $bin, "samtools view $input | ") or die "ERROR: $!\n";
while (my $line = <$bin>){
    chomp $line;
    my $ref = (split("\t", $line))[2];
    unless ($rseen{$ref}){
        push @refs, $ref;
    }
    $rseen{$ref} = 1;
}

if ($opt_n){
    my $nseen;
    foreach (@refs){
        $nseen = 1 if $_ eq $opt_n;
    }
    die "ERROR: reference '$opt_n' not found in $input\n" unless $nseen;
}

my ($reads_in, $barcodes_out) = (0) x 2;
foreach my $ref (@refs){
    open (my $in, "samtools view $input $ref:$istart-$istop | ") or die "ERROR: $!\n";
    my %hash;
    my %seen;
    while (my $line = <$in>){
        chomp $line;
        $reads_in++;
        my ($id, $flag, $start, $cigar, $seq) = (split("\t", $line))[0,1,3,5,9];
        next if $seen{$id};
        my $iseq = "";
        ## add some slop to both of the stops below? I.e. make sure the read starts or ends no less than 3 or 5 bp up or down from the istart / istop?
        next if $start > ($istart - $slop);
        next if ($start + length($seq)) - 1 < ($istop + $slop);
        #next if $flag & 2048; #skip supplementary alignments. This is a blunt criteron, may miss some good alignments, but primary alignment is usually longest and best
        my @cigars = split(/(?<=\D)/, $cigar);
        my $refpos = $start; ## 1-based 
        my $readpos = 0; ## zero-based to make substring operation easier
        my $is_first; 
        foreach my $part (@cigars){
            if ($part =~ m/(\d+)([MIDSH])/){
                my ($num, $op) = ($1, $2);
                next if $op eq "H"; ## skip hard-clipped reads as the read sequence won't be in the BAM file.
                if ($op eq "S"){
                    if (!$is_first){
                        $is_first = 1;
                        $refpos -= $num;
                        last if $refpos + $num >= ($istart - $slop); ## stop if soft clipping goes into the barcode + slop. 
                        next unless $refpos <= ($istart - $slop);
                        for my $i (1 .. $num){
                            $refpos++;
                            $readpos++;
                            last if $refpos == $istart;  ## unnecessary?
                        } 
                    } else {
                        last;
                    }
                }
                if ($op =~ m/[MID]/){
                    for my $i (1 .. $num){
                        $refpos++ unless $op eq "I";
                        $readpos++ unless $op eq "D";
                        if ((sort{$a <=> $b}($refpos, $istart, $istop))[1] == $refpos){
                            my $base = substr($seq, $readpos, 1);
                            $iseq .= $base unless $op eq "D" and $i != $num; ## second part of the and statement is to include first base overlapping barcode after a deletion (from trial and error)
                        }
                    }
                }
                if  ($op eq "N"){
                    print STDERR "uh oh. 'N' cigar.\n";
                }
            }
        }
        next if ($refpos <= $istop + $slop);
        if ($iseq){
            if (length($iseq) == $ileng){
                $barcodes_out++;
                $hash{$iseq}++;
                $seen{$id} = 1;
                print STDERR "$iseq\n" if $id eq $debug;
            }
        }
    }
    close ($in);

    print "ref\tbarcode\tcount\n";
    foreach my $bar (sort keys %hash){
        my $val = $hash{$bar};
        print "$ref\t$bar\t$val\n";
    }

    print STDERR "$ref\treads:$reads_in\treads_w_bc:$barcodes_out\n";
}

