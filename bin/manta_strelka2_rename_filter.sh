#! /usr/bin/env bash
SAMPLEID=$3
export SAMPLEID

if [[ $1 =~ ".gz$" ]];then
  gunzip -c $1 | \
  perl -ane 'if($F[0]=~m/^#/){if($_=~m/^#CHROM/){
      $_=~s/TUMOR/$ENV{'SAMPLEID'}/;
      print $_;next;
    }
      else{print $_;next;}
    }
  else{print $_;}' > $2
else
  cat $1 | \
  perl -ane 'if($F[0]=~m/^#/){if($_=~m/^#CHROM/){
      $_=~s/TUMOR/$ENV{'SAMPLEID'}/;
      print $_;next;
    }
      else{print $_;next;}
    }
  else{print $_;}' > $2
fi
