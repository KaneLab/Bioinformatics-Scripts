# filters blast table for single-copy genes. This is for Kyle Keeper's pipeline using a SPADES assembly of lichen WGS blasted against fungal BUSCO genes
# column1 is the SPADES assembly contig, column 2 is the fungal BUSCO hit. Column 12 is the bitscore
# usage:
# perl filter_blast_table.pl BLAST_OUT OUTFILE
$BLAST = $ARGV[0];

$OUT = $ARGV[1];


$length_tot=0;
$count=0;
$lengthxcov=0;
@lengths = ();
@covs = ();

open BLAST or die "No file $BLAST";
while (<BLAST>){
        @tabs = split "\t", $_;
	@tabs2 = split "_", $tabs[0];
	$count=$count+1;
	$length_tot=$length_tot + $tabs2[3];
	$lengthxcov=$lengthxcov + ($tabs2[3] * $tabs2[5]);
	push @lengths, $tabs2[3];
	if ($tabs2[3]>1000){
		push @covs, $tabs2[5];	
	}
}
close BLAST;

$avg_length = $length_tot / $count;
$coverage = $lengthxcov / $length_tot;

print "$count contigs in blast table. Average length is $avg_length. Average coverage is $coverage.\n";
print "Median length is ";
print median(@lengths);
print ". Median coverage is ";
print median(@covs);
print ".\n";

open OUTFILE, ">$OUT";

$mediancovs=median(@covs);
$min = 0.8 * $mediancovs;
$max = 1.2 * $mediancovs;
print "Minimum coverage to keep is $min and maximum is $max.\n";

$last = "";
$count = 0;
open BLAST or die "No file $BLAST";
while (<BLAST>){
	@tabs = split "\t", $_;
        @tabs2 = split "_", $tabs[0];
        if ($tabs2[3]>1000 && $tabs2[5]>$min && $tabs2[5]<$max){
                print OUTFILE $_;
        }

}
close BLAST;
close OUTFILE;








#######################################

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

