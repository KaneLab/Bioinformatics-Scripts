# reverse complement a fastq file

$seq = $ARGV[0]; # fastq file

        $rc_seq = reverse $seq;
        $rc_seq =~ tr/ATGCatgc/TACGtacg/;
        $r_qual = reverse $qual;

        print "$rc_seq";

exit;
