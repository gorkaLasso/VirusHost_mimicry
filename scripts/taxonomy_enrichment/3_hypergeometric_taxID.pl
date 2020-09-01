#!/usr/bin/perl

=description
1. get the user-defined input directory
2. get the background size for each taxonomy category
3. read the protein neighbor file
   3.1 count the total number of protein neighbors
   3.2 count the total number of proteins within each taxonomy category
4. compute the hypergeometric test and store all pvalues in an array
5. correct the p-values for multiple testing
=cut

use strict;
use warnings;
use Statistics::R;
my $R = Statistics::R->new();

################################################################################
#    GLOBAL
################################################################################

my $backgroundFile = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/Databases/Skan_summary/6_taxid_prevalence.txt";
$backgroundFile =~ s/\/home\/gl2411// if (!-e $backgroundFile);
my $uniprotNeigFilePreffix = "uniprot_neighbor_space";

################################################################################
#    MAIN
################################################################################

&Main;

################################################################################
#    FUNCTIONS
################################################################################

sub Main{
  my ($inputDir,$sasThresh) = &GetInputOptions(@ARGV);
  print "inputDir\t$inputDir\n";
  my ($backSize,%backCat,%queryCat,%taxid2os);
  &ReadBackground($backgroundFile,\$backSize,\%backCat,\%queryCat,\%taxid2os);
  
  my $querySize;
  &ReadProtNeig("$inputDir\/$uniprotNeigFilePreffix\_$sasThresh\.txt",\$querySize,\%queryCat);
  
  my %commandList;
  my @pvalList;
  my %cat2pval;
  foreach my $cat (sort keys %backCat){
    my $catQuerySize = $queryCat{$cat};
    my $catQueryPercen = $catQuerySize / $querySize;
    my $catBackSize = $backCat{$cat};
    #next if ( (not defined $backSize) || ($catBackSize == 0) || ($backSize == 0) || ($catQuerySize == 0));
    my $catBackPercen = $catBackSize / $backSize;
    #print "\t$cat\t" . sprintf("%.3f",$catQueryPercen) . "\t" . sprintf("%.3f",$catBackPercen) . "\n";
    my ($pval,$command) = &ComputeHypergeometricPval($catQuerySize,$catBackSize,$backSize-$catBackSize,$querySize);
    $commandList{$cat} = $command;
    push(@pvalList,$pval);
    #print "pval is $pval\n";
  }
  #print "size of pval arr is " . scalar(@pvalList) . "\n";
  my @pvalCorList = &CorrectPvalues(\@pvalList);
  
  my $outFile = $inputDir . "/enrichment_" . $sasThresh . ".txt";
  open (OUTFILE,">",$outFile) or die $!;
  print OUTFILE "background size: $backSize\n";
  foreach my $cat (sort keys %backCat){
    my $command = $commandList{$cat};
    print OUTFILE "$cat\t$command\n";
  }
  print OUTFILE "\n";
  
  print OUTFILE "#organism\ttaxid\tNeigProtsInCat\tAllProtsInCat\tAllProts\tAllNeig\tObservedPercen\tExpectedPercen\tpval\tcorrectedPval\t-log10CorrPval\n";
  my $ct = 0;
  foreach my $cat (sort keys %backCat){
    my $catQuerySize = $queryCat{$cat};
    my $catQueryPercen = $catQuerySize / $querySize;
    my $catBackSize = $backCat{$cat};
    #next if ( (not defined $backSize) || ($catBackSize == 0) || ($backSize == 0) || ($catQuerySize == 0));
    my $catBackPercen = $catBackSize / $backSize;
    
    my $pval = $pvalList[$ct];
    my $pvalCor = $pvalCorList[$ct];
    my $lg = 10000;
    $lg = -log($pvalCor)/log(10) if ($pvalCor > 0);
    
    my $osName = $taxid2os{$cat};
    
    print OUTFILE "$osName\t$cat\t$catQuerySize\t$catBackSize\t$backSize\t$querySize\t" . sprintf("%.3f",$catQueryPercen) . "\t" . sprintf("%.3f",$catBackPercen) ."\t" . $pval . "\t" . $pvalCor . "\t" . $lg . "\n";
    $ct++;
  }
  
  close OUTFILE;
}

################################################################################

sub CorrectPvalues{
  my @pvalList = @{$_[0]};
  #print "pvalues: " . scalar(@pvalList) . "\n";
  $R->set('arr2',"nan");
  $R->set('arr2',\@pvalList);
  $R->run(
  q`correct<-p.adjust(arr2,method="bonferroni")`
  );
  my $corrected = $R->get('correct');
  #foreach my $item (@{$corrected}){
  #  print "$item\n";
  #}
  return(@{$corrected});
}

################################################################################

sub ComputeHypergeometricPval{
  my $q = $_[0];
  my $m = $_[1];
  my $n = $_[2];
  my $k = $_[3];
  
  #print "q is $q\n";
  
  $R->set('q',$q-1);
  $R->set('m',$m);
  $R->set('n',$n);
  $R->set('k',$k);
  my $command = "phyper(" . ($q-1) . ",$m\,$n\,$k\,lower.tail=FALSE)";
  my $out1 = $R->run(
  q`phyper(q,m,n,k,lower.tail=FALSE)`
  );
  my ($temp,$pval) = split(/\s/,$out1);
  #print "joder..... #$pval\#\n";
  return($pval,$command);
}

################################################################################

sub ReadProtNeig{
  my $file = $_[0];
  my $querySizeRef = $_[1];
  my $queryCatRef = $_[2];
  
  open (INFILE,"<",$file) or die "cannot open $file\n";
  
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    next if (length($s) < 5);
    my ($uni,$taxid,$os,$taxGroup) = split(/\t/,$s);
    $$querySizeRef++ if ($taxid ne "nan");
    if (defined ${$queryCatRef}{$taxid}){
      ${$queryCatRef}{$taxid}++;
    }
    
  }
  
  close INFILE;
  
}

################################################################################

sub ReadBackground{
  my $file = $_[0];
  my $backSizeRef = $_[1];
  my $backCatRef = $_[2];
  my $queryCatRef = $_[3];
  my $taxid2osRef = $_[4];
  
  open (INFILE,"<",$file) or die "cant open $file\n";
  
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    next if (substr($s,0,1) eq "#");
    last if (substr($s,0,2) eq "//");
    next if (length($s)<5);
    my @w = split(/\t/,$s);
    if ($w[0] eq "dbSize"){
      $$backSizeRef = $w[1];
    }elsif ($w[0] eq "taxGroup"){
      ${$backCatRef}{$w[1]} = $w[3];
      ${$queryCatRef}{$w[1]} = 0;
      ${$taxid2osRef}{$w[1]} = $w[2];
    }
  }
  
  close INFILE;
}

################################################################################

sub GetInputOptions{
    my @arguments = @_;
    my $PGM;
    my $usage;
    my $s;
    my $inputDir;
    my $sasThresh = 3.5;
    
    $PGM = $0;
    $PGM =~ s#.*/##;                #remove part up to last slash
    
$usage = <<USAGE;
Usage:
$PGM -d input directory -s 
  -d 
  -s SAS threshold (optional) by default -s 3.5
USAGE

    while(@arguments){
        $s = shift(@arguments);
        if ($s){
            if ($s eq "-d") {$inputDir = shift(@arguments); next;}
            if ($s eq "-s") {$sasThresh = shift(@arguments); next;}
            }
    }
    if (!$inputDir) {&PrintCommandError("error: input params not defined",$usage)}
    return ($inputDir,$sasThresh);
}

################################################################################

sub PrintCommandError{
    my $errorMessage = $_[0];
    my $usage = $_[1];
    print "\n". $errorMessage . "\n";
    print "\n". $usage . "\n";
    exit;
}
