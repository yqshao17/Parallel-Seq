infile=$1
outfile=$2
cut -f 1-5 $infile | uniq | cut -f 1-4 |sort| uniq -c |awk '{printf ("%s\t%d\t%d\t%d\t%s\n",$2,$3,$4,$1,$5)}'>$outfile