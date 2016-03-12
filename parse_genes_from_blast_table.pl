#Takes a blast table in outfmt 6 and searches against a fasta file, and produces a filtered fasta file containing only the subsequences that are blast hits.
# usage :
# perl parse_genes_from_blast_table.pl BLASTOUTPUT_OUTFMT6 FASTA OUTFILE
$BLAST = shift @ARGV;
$FASTA = shift @ARGV;
$OUT = shift @ARGV;
@ids = ();
%fasta = ();
@fasta2 = ();
@seqs = ();
%subseqs = ();

#open BLAST or die "No file $BLAST";		#READ ID into an array
#while (<BLAST>){
#    @tabs = split "\t", $_;
#	push (@ids, $tabs[0]);
#}
#close BLAST;


#then do the fasta loop with undefine
#then a foreach loop with id array grep against @fasta

open FASTA or die "No file $FASTA";

for ($junk = 0; $junk <1;  $junk++){
	$/ = '>';
	while (<FASTA>) {
	        $_ =~ s/>//g;
		push @fasta2, $_;
	}	
}
close FASTA;

shift @fasta2;
$/ = "\n";
foreach $seq (@fasta2) {
	@tabs = split "\n", $seq;
	$name = shift @tabs;
	@tabs2 = split " ", $name;
	$name2 = shift @tabs2;
	$seq2 = join '', @tabs;
	#print "$name2\n";
	$fasta{$name2} = $seq2;
}



$last = "";
$count = 0;
$seq = "";
%allseqs = ();
$busconame="";
open BLAST or die "No file $BLAST";
while (<BLAST>){
	@tabs = split "\t", $_;
	$busconame=$tabs[1];	

	if (exists $fasta{$tabs[0]}){
		if ($busconame eq $last){
                        $count=$count+1;
                }
		elsif ($seq eq ""){}
                else {
 #                       print OUTFILE '>', "$tabs[1]_$count\n$seq\n";
  			if(length($seq) > length($allseqs{$busconame})){
				$allseqs{$last}=$seq;
			}			
                        $count=0;
                        $seq = "";
                }

		if ($tabs[6] < $tabs[7]){
			$subseq = substr $fasta{$tabs[0]}, $tabs[6] -1, ($tabs[7] - $tabs[6] +1);
			$seq = "$seq$subseq";
		}
		else{
			$subseq = reverse (substr $fasta{$tabs[0]}, $tabs[7] -1, (1+ $tabs[6] - $tabs[7] ));
			$subseq =~ tr /atcgATCG/tagcTAGC/;
                        $seq = "$subseq$seq";
		}
		$last = $busconame;	
	}	
	else {        print "$tabs[0] has no corresponding fasta in your fasta file.\n";
	}
}
 #                      print OUTFILE '>', "$tabs[1]_$count\n$seq\n";

if(length($seq) > length($allseqs{$busconame})){
      $allseqs{$busconame}=$seq;
}

open OUTFILE, ">$OUT";
foreach (sort keys %allseqs) {
    print OUTFILE '>', "$_\n$allseqs{$_}\n";
}	
close OUTFILE;



