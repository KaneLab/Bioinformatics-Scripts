# reverse complement a fastq file

$INFILE = $ARGV[0]; # fastq file


open INFILE;
while(<INFILE>){
        chomp $_;
        $header = $_;
        $seq = (<INFILE>);
        chomp $seq;
        $temp = (<INFILE>);
        $qual = (<INFILE>);
        chomp $qual;

        $rc_seq = reverse $seq;
        $rc_seq =~ tr/ATGCatgc/TACGtacg/;
        $r_qual = reverse $qual;

        print "$header\n$rc_seq\n+\n$r_qual\n";

}

close INFILE;
exit;
