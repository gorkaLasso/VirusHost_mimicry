#!/usr/bin/perl

=description

=cut

use strict;
use warnings;

##################################################################################
#    Globals
##################################################################################
my $filesuffix = "_david_2.5.txt";
my %fam2file = (
#  "AallVir" => "all_humanVir" . $filesuffix,
  "dsDNA Adenoviridae" => "adenoviridae" . $filesuffix,
  "n_ssRNA Arenaviridae" => "arenaviridae" . $filesuffix,
  "p_ssRNA Astroviridae" => "astroviridae" . $filesuffix,
  "p_ssRNA Caliciviridae" => "caliciviridae" . $filesuffix,
  "ssDNA Circoviridae" => "circoviridae" . $filesuffix,
  "p_ssRNA Coronaviridae" => "coronaviridae" . $filesuffix,
  "n_ssRNA Filoviridae" => "filoviridae" . $filesuffix,
  "p_ssRNA Flaviviridae" => "flaviviridae" . $filesuffix,
  "dsDNA Herpesviridae" => "herpesviridae" . $filesuffix,
  "p_ssRNA Hepeviridae" => "hepeviridae" . $filesuffix,
  "n_ssRNA Orthomyxoviridae" => "orthomyxoviridae" . $filesuffix,
  "n_ssRNA Paramyxoviridae" => "paramyxoviridae" . $filesuffix,
  "ssDNA Parvoviridae" => "parvoviridae" . $filesuffix,
  "p_ssRNA Picornaviridae" => "picornaviridae" . $filesuffix,
  "n_ssRNA Pneumoviridae" => "pneumoviridae" . $filesuffix,
  "dsDNA Polyomaviridae" => "polyomaviridae" . $filesuffix,
  "dsDNA Poxviridae" => "poxviridae" . $filesuffix,
  "dsRNA Reoviridae" => "reoviridae" . $filesuffix,
  "n_ssRNA Rhabdoviridae" => "rhabdoviridae" . $filesuffix,
  "ssDNA Smacoviridae" => "smacoviridae" . $filesuffix,
  "p_ssRNA Togaviridae" => "togaviridae" . $filesuffix,
);

my %cat2class =(
  "GAD_DISEASE" => "disease",
  "GAD_DISEASE_CLASS" => "disease",
  "OMIM_DISEASE" => "disease",
  "GOTERM_MF_DIRECT" => "molfun",
  "GOTERM_BP_DIRECT" => "pathway",
  "BIOCARTA" => "pathway",
  "KEGG_PATHWAY" => "pathway",
  "REACTOME_PATHWAY" => "pathway",
  "INTERPRO" => "domain",
  "PIR_SUPERFAMILY" => "domain",
);

my %cat2abv =(
  "GAD_DISEASE" => "GD_",
  "GAD_DISEASE_CLASS" => "GC_",
  "OMIM_DISEASE" => "OM_",
  "GOTERM_MF_DIRECT" => "GO_",
  "GOTERM_BP_DIRECT" => "GO_",
  "BIOCARTA" => "BC_",
  "KEGG_PATHWAY" => "KG_",
  "REACTOME_PATHWAY" => "",
  "INTERPRO" => "IN_",
  "PIR_SUPERFAMILY" => "PR_",
);

my $pvalThresh = 0.01;
my $targetClass = "molfun";

##################################################################################
#    Main
##################################################################################
my $outputFile = "matrix_$targetClass\_$pvalThresh\.txt"; 
&Main;

##################################################################################
#    Functions
##################################################################################

sub Main{
  my %targetTerm; #terms enriched with a corrected pval < pvalThresh
  my %term2family; # terms to family to -log(10) Pval
  my %reactomeMapping = &ReadReactomeFile("reactome.txt");
  my %fam2sigPathCt;
  foreach my $family (sort keys %fam2file){
    $fam2sigPathCt{$family} = 0;
    #next if ($family ne "Adenoviridae");
    my $file = $fam2file{$family};
    die "cant find $file" if (!-e $file);
    open (INFILE,"<",$file) or die $!;
    my $sigPathCt = 0;
    while(<INFILE>){
      my $s = $_;
      chomp($s);
      my @w = split(/\t/,$s);
      next if ($w[0] eq "Category");
      die "wrong number of columns in file $file\n" if (scalar(@w) != 13);
      my $cat = $w[0];
      my $term = $w[1];
      my $corrPval = $w[10];
      next if (not defined $cat2class{$cat});
      my $class = $cat2class{$cat};
      next if ($class ne $targetClass);
      $sigPathCt++ if ($corrPval < $pvalThresh);
      if ($term =~ m/R-HSA-/){
        $term =~ s/\:.*//;
      }
      $term =~ s/\"//g;
      $term = $cat2abv{$cat} . $term;
      my $logPval = -log($corrPval)/log(10);
      #print "$family\t$class\t$term\t$corrPval\t$logPval\n";
      $term2family{$term}{$family} = $logPval;
      $targetTerm{$term}++ if ($logPval > 2);
      
    }
    close INFILE;
    print "$family\t$sigPathCt significant terms\n";
    $fam2sigPathCt{$family} = $sigPathCt;
  }
  
  print "\nSignificant terms in >=1 viral family: " . keys(%targetTerm) . "\n";
  
  open (OUTFILE,">",$outputFile) or die $!;
  
  my $s = "term";
  foreach my $family (sort keys %fam2file){
    next if ($fam2sigPathCt{$family} == 0); #not including the viral families for which there is not a single enriched term
    $s .= "\t$family";
  }
  print OUTFILE "$s\n";
  foreach my $term (sort keys %targetTerm){
    my $editTerm = $term;
    $editTerm =~ s/GO\:.*\~//;
    $editTerm =~ s/hsa.*\://;
    $editTerm =~ s/h_.*\://;
    $editTerm =~ s/IPR.*\://;
    $editTerm =~ s/PIR.*\://;
    if (defined $reactomeMapping{$term}){
      $editTerm = "RE_" . $reactomeMapping{$term};
    }
    my $s = $editTerm;
    foreach my $family (sort keys %fam2file){
      next if ($fam2sigPathCt{$family} == 0); #not including the viral families for which there is not a single enriched term
      if (not defined $term2family{$term}{$family}){
        $s .= "\t0.0";
      }else{
        my $logPval = $term2family{$term}{$family};
        $s .= "\t" . sprintf("%.1f",$logPval);
      }
    }
    print OUTFILE "$s\n";
  }
  
  close(OUTFILE);
}

##################################################################################

sub ReadReactomeFile{
  my $file = $_[0];
  my %term2name;
  
  open(INFILE,"<",$file) or die $!;
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    my ($term,$name,$os) = split(/\t/,$s);
    $term2name{$term} = $name;
  }
  close INFILE;
  
  return(%term2name);
}