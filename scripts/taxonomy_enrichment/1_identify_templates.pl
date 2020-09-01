#!/usr/bin/perl

=description

1. get the user-defined host taxonomic group
2. create dir if does not exist
3. identify all viral proteins belonging to viruses infecting the host taxonomic group
4. identify all modeling templates
5. identify the subset of modeling templates used to run skan
   5.1 write output

host taxonomic groups 
L1: Bacteria
L1: Invert&Bacteria
L1: Invertebrate
L1: Invert&PlantFungi
L1: Invert&Vert
L1: PlantFungi
L1: Vertebrate
L2: Invert&Mammal
L2: Mammal
L2: NonMammal
L3: Human
L3: Invert&Human
L3: Invert&NonHuman
L3: NonHuman



=cut

use strict;
use warnings;

################################################################################
#    GLOBALS
################################################################################

my $virProtMapFile = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/PHI/virushostdb/viral_mimicry/0_virushostDB_by_hostDiv/2_merged_all/glcA_mapping_2.txt";
my $virProtMapWithRedundantFile = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/PHI/virushostdb/viral_mimicry/0_virushostDB_by_hostDiv/2_merged_all/glcA_mapping_3_withRedundant.txt";
#my $allVirProtMapFile = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/PHI/virushostdb/viral_mimicry/0_virushostDB_by_hostDiv/2_merged_all/merged_fasta_mapping.txt";
#my $clusterFile = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/PHI/virushostdb/viral_mimicry/0_virushostDB_by_hostDiv/2_merged_all/merged_fasta_nr100.txt.clstr";

my $genomeDir = "/home/gl2411/ifs/data/c2b2/bh_lab/shares/databases/hfpd/genomes/";
my $templateFile = "templates/templates_query_1e+00_1e+00.txt";
my $templateMapFile = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/PHI/virushostdb/viral_mimicry/scr/merge_templates/templates_query_1e+00_1e+00_merged_mapping.txt";
my $pdb2uniSkan60File = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/Databases/Skan_summary/2_pdb_in_skan60_withUni.txt";
my $uni2hostGroupFile = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/Databases/Skan_summary/5_unique_uniprots_inSKAN60_withTaxid_Group.txt";

################################################################################
#    MAIN
################################################################################

&Main;

################################################################################
#    FUNCTION
################################################################################

sub Main{
  my ($hostGroup,$onlyVirtemp,$hostTaxid,$virTaxid) = &GetInputOptions(@ARGV);
  if (defined $hostGroup){
    print "taxonomy group of host: $hostGroup | only viral templates: $onlyVirtemp\n";
  }else{
    print "taxonomy id of host: $hostTaxid | only viral templates: $onlyVirtemp\n" if (defined $hostTaxid);
    print "taxonomy id of virus: $virTaxid | only viral templates: $onlyVirtemp\n" if (defined $virTaxid);
  }
  #die;
  my $newDir;
  if (defined $hostGroup){
    $newDir = "$hostGroup\_$onlyVirtemp";
    $newDir =~ s/\&/\_/g;
  }elsif (defined $hostTaxid){
    $newDir = "host_$hostTaxid\_$onlyVirtemp";
  }elsif (defined $virTaxid){
    $newDir = "vir_$virTaxid\_$onlyVirtemp";
  }
  
  system("mkdir $newDir") if (!-d $newDir);
  my $summaryOutFile = $newDir . "/template_summary.txt";
  open (SUMOUTFILE,">",$summaryOutFile) or die $!;
  
  #get the subset of viral proteins
  my @glcbList;
  if (defined $hostGroup){
    @glcbList = &GetProteins($virProtMapFile,$hostGroup);
  }elsif (defined $hostTaxid){
    @glcbList = &GetProteinsFromHostTaxid($virProtMapWithRedundantFile,$hostTaxid); #ATTENTION! IT ALLOWS REDUNDANT PROTEINS TO BE CONSIDERED
  }elsif (defined $virTaxid){
    @glcbList = &GetProteinsFromViralTaxid($virProtMapWithRedundantFile,$virTaxid); #ATTENTION! IT ALLOWS REDUNDANT PROTEINS TO BE CONSIDERED
  }
  print SUMOUTFILE "viral proteins: " . scalar(@glcbList) . "\n";

  #get pdb2uniprot and uniprot2hostgroup mappings
  my (%pdb2uniprot,%viralUniprot,%uniprot2hostGroup);
  &ReadPdb2Uniprot($pdb2uniSkan60File,\%pdb2uniprot);
  &ReadViralUniprot($uni2hostGroupFile,\%viralUniprot,\%uniprot2hostGroup);
  
  #read templates
  my (%tempall,%tempFile,%glcbWithTemp);
  my $ct = 0;
  foreach my $glcb (@glcbList){
    $ct++;
    print "processing $ct / " . scalar(@glcbList) . "\n" if ($ct%1000==0);
    my ($genomeName,$hfpd) = &ParseGlcb($glcb);
    #print "$glcb\t$genomeName\t$hfpd\n";
    &GetTemplates($genomeDir,$genomeName,$hfpd,$templateFile,\%tempall,\%tempFile,\%glcbWithTemp,$onlyVirtemp,\%pdb2uniprot,\%viralUniprot);
  }
  print SUMOUTFILE "viral proteins with templates: " . keys(%glcbWithTemp) . "\n";
  print SUMOUTFILE "templates (pdb_chain_start_end): " . keys(%tempall) . "\n";
  open (HFPDWITHTEMP,">","$newDir\/virprotsWithTemplate.txt") or die "I could not open file ";
  foreach my $glcb (keys %glcbWithTemp){
    print HFPDWITHTEMP "$glcb\n";
  }
  close(HFPDWITHTEMP);
  
  
  #get mapping template to representative template
  my (%temp2rep,%reptemp);
  &ReadTempMapFile($templateMapFile,\%temp2rep);
  
  #get the list of representative templates used for skan search
  foreach my $temp (sort keys %tempall){
    #print "$temp\n";
    die "this templates is not in the mapping file $temp\n" if (not defined $temp2rep{$temp});
    my $rep = $temp2rep{$temp};
    $reptemp{$rep} = 1;
  }
  print SUMOUTFILE "representative templates (pdb_chain_start_end): " . keys(%reptemp) . "\n";
  
  my (%tempHostGroup2uni);
  
  #print representative templates with their uniprot and hostGroup
  my $outFile = $newDir . "/representative_template_space.txt";
  open (OUTFILE,">",$outFile) or die $!;
  foreach my $temp (sort keys %reptemp){
    my @w = split(/\t/,$temp);
    my $pdb = $w[0];
    #print "$pdb\n";
    my $uni = "nan";
    $uni = $pdb2uniprot{$pdb} if (defined $pdb2uniprot{$pdb});
    my $hostGroup = "nan";
    $hostGroup = $uniprot2hostGroup{$uni} if (defined $uniprot2hostGroup{$uni});
    #print "$uni\n";
    print OUTFILE "$temp\t$uni\t$hostGroup\n";
    $tempHostGroup2uni{$hostGroup}{$uni} = 1 if ($uni ne "nan");
  }
  close OUTFILE;
  
  #print a summary of uniprots founds as templates
  my $uniCt = 0;
  foreach my $hostGroup (sort keys %tempHostGroup2uni){
    $uniCt = $uniCt + keys(%{ $tempHostGroup2uni{$hostGroup} });
  }
  print SUMOUTFILE "Unique uniprots used as templates: $uniCt\n";
  foreach my $hostGroup (sort keys %tempHostGroup2uni){
    my $ct = keys(%{ $tempHostGroup2uni{$hostGroup} });
    print SUMOUTFILE $hostGroup . "\t" . $ct  . "\t" . sprintf("%.2f",$ct/$uniCt) . "\n";
  }

  close SUMOUTFILE;
  
}

################################################################################

sub ReadViralUniprot{
  my $file = $_[0];
  my $viralUniRef = $_[1];
  my $uni2hostGroupRef = $_[2];
  
  open (INFILE,"<",$file) or die "cant open $file\n";
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    next if (length($s) < 5);
    #print "$s\n";
    my ($uni,$taxid,$os,$hostGroup) = split(/\t/,$s);
    #print "uni $uni\n";
    #print "taxid $taxid\n";
    ${$uni2hostGroupRef}{$uni} = $hostGroup;
    next if ($hostGroup ne "virus");
    ${$viralUniRef}{$uni} = 1;
  }
  close(INFILE);
}

################################################################################

sub ReadPdb2Uniprot{
  my $file = $_[0];
  my $pdb2uniRef = $_[1];
  
  open (INFILE,"<",$file) or die $!;
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    my ($pdb,$uni) = split(/\t/,$s);
    ${$pdb2uniRef}{$pdb} = $uni;
  }
  close INFILE;
}

################################################################################

sub ReadTempMapFile{
  my $file = $_[0];
  my $temp2repRef = $_[1];
  
  open (INFILE,"<",$file) or die "cant open $file\n";
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    next if (substr($s,0,1) eq "#");
    next if (length($s)<5);
    my ($temp,$tempStart,$tempEnd,$rep,$repStart,$repEnd) = split(/\t/,$s);
    ${$temp2repRef}{"$temp\t$tempStart\t$tempEnd"} = "$rep\t$repStart\t$repEnd";
  }
  close (INFILE);
}

################################################################################

sub GetTemplates{
  my $genomeDir = $_[0];
  my $genomeName = $_[1];
  my $hfpd = $_[2];
  my $templateFile = $_[3];
  my $tempallRef = $_[4];
  my $tempFileRef = $_[5];
  my $hfpdWithTempRef = $_[6];
  my $onlyVirtemp = $_[7];
  my $pdb2uniRef = $_[8];
  my $viralUniRef = $_[9];
  
  my $fullFile = $genomeDir . $genomeName . "/" . $templateFile;
  
  if (not defined ${$tempFileRef}{$fullFile}){
    if (-e $fullFile){
      ${$tempFileRef}{$fullFile} = 1;
    }else{
      die "cant find template file $fullFile\n";
    }
  }
  
  open (INFILE,"<",$fullFile) or die $!;
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    my ($query,$method,$queryRegion,$template,$xrayStart,$xrayEnd,$eval,$evalCutOff,$queryStart,$queryEnd) = split(/\t/,$s);
    #print "comparing $query\tVs\t$hfpd\n";
    next if (substr($query,0,6) ne $hfpd); # I also need to consider the domains not just the full sequences!
    #print "yesss:\t$query\t$hfpd\n";
    my $templ = $template . "\t" . $xrayStart . "\t" . $xrayEnd;
    if (lc($onlyVirtemp) eq "false"){
      ${$tempallRef}{$templ} = 1;
      ${$hfpdWithTempRef}{"$genomeName\_$hfpd"} = 1;
    }else{
      next if (not defined ${$pdb2uniRef}{$template});
      my $uni =  ${$pdb2uniRef}{$template};
      next if (not defined ${$viralUniRef}{$uni});
      ${$tempallRef}{$templ} = 1;
      ${$hfpdWithTempRef}{"$genomeName\_$hfpd"} = 1;
    }
  }
  close INFILE;
  
}

################################################################################

sub ParseGlcb{
  my $glcb = $_[0];
  my $hfpd = substr($glcb,-6);
  my $genomeName = $glcb;
  $genomeName =~ s/\_$hfpd//;
  return($genomeName,$hfpd);
}

################################################################################

sub GetProteinsFromHostTaxid{
  my $file = $_[0];
  my $hostTaxid = $_[1];
  
  my @glcbList;
  
  open (INFILE,"<",$file) or die $!;
  
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    #print "$s\n";
    next if (substr($s,0,1) eq "#");
    my ($glca,$glcb,$taxid,$l1,$l2,$l3,$hostTaxids,$hostNames) = split(/\t/,$s);
    #print "$l3\n";
    $glcb =~ s/repr_//;
    my $temp = "#" . $hostTaxid . "#";
    next if ($hostTaxids !~ m/$temp/);
    push(@glcbList,$glcb);
  }
  
  close(INFILE);
  
  return(@glcbList);
}

################################################################################

sub GetProteinsFromViralTaxid{
  my $file = $_[0];
  my $virTaxid = $_[1];
  
  my @glcbList;
  
  open (INFILE,"<",$file) or die $!;
  
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    #print "$s\n";
    next if (substr($s,0,1) eq "#");
    my ($glca,$glcb,$taxid,$l1,$l2,$l3,$hostTaxids,$hostNames) = split(/\t/,$s);
    #print "$l3\n";
    $glcb =~ s/repr_//;
    #my $temp = "#" . $hostTaxid . "#";
    #next if ($hostTaxids !~ m/$temp/);
    next if ($taxid ne $virTaxid);
    push(@glcbList,$glcb);
  }
  
  close(INFILE);
  
  return(@glcbList);
}

################################################################################

sub GetProteins{
  my $file = $_[0];
  my $hostGroup = $_[1];
  my @glcbList;
  
  open (INFILE,"<",$file) or die $!;
  
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    next if (substr($s,0,1) eq "#");
    my ($glca,$glcb,$taxid,$l1,$l2,$l3) = split(/\t/,$s);
    if ( (uc($l1) eq uc($hostGroup)) || (uc($l2) eq uc($hostGroup)) || (uc($l3) eq uc($hostGroup)) ){
      push(@glcbList,$glcb);
    }else{
      next;
    }
  }
  
  close(INFILE);
  
  return(@glcbList);
}

################################################################################

sub GetInputOptions{
    my @arguments = @_;
    my $PGM;
    my $usage;
    my $s;
    my $hostGroup;
    my $onlyVirtemp = "false";
    my $hostTaxid;
    my $virTaxid;
    
    $PGM = $0;
    $PGM =~ s#.*/##;                #remove part up to last slash
    
$usage = <<USAGE;
Usage:
$PGM -g host group -t hostTaxid -p viralTaxid -v false/true
  host group options : Bacteria | Invert\\&Bacteria | Invertebrate
  (case insentitive)   Invert\\&PlantFungi | Invert\&Vert | PlantFungi
                       Vertebrate | Invert\\&Mammal | Mammal | NonMammal
                       Human | Invert\\&Human | Invert\\&NonHuman | NonHuman
  -v                 : true (only viral templates) / false
USAGE

    while(@arguments){
        $s = shift(@arguments);
        if ($s){
            if ($s eq "-g") {$hostGroup = shift(@arguments); next;}
            if ($s eq "-v") {$onlyVirtemp = shift(@arguments); next;}
            if ($s eq "-t") {$hostTaxid = shift(@arguments); next;}
            if ($s eq "-p") {$virTaxid = shift(@arguments); next;}
            }
    }
    $onlyVirtemp = lc($onlyVirtemp);
    if ( (not defined $hostGroup) && (not defined $hostTaxid) && (not defined $virTaxid) ) {&PrintCommandError("error: input params not defined",$usage)}

    if ( ($hostGroup) && ($hostTaxid) && ($virTaxid) ) {&PrintCommandError("error: you cannot only defined one out of hostGroup (-g), hostTaxid (-t), viralTaxid (-p)",$usage)}
    if ( ($hostGroup) && ($hostTaxid) ) {&PrintCommandError("error: you cannot only defined one out of hostGroup (-g), hostTaxid (-t), viralTaxid (-p)",$usage)}
    if ( ($hostGroup) && ($virTaxid) ) {&PrintCommandError("error: you cannot only defined one out of hostGroup (-g), hostTaxid (-t), viralTaxid (-p)",$usage)}
    if ( ($hostTaxid) && ($virTaxid) ) {&PrintCommandError("error: you cannot only defined one out of hostGroup (-g), hostTaxid (-t), viralTaxid (-p)",$usage)}

    if ( ($onlyVirtemp ne "false") && ($onlyVirtemp ne "true" ) ) {&PrintCommandError("-v can only be true or false",$usage)}
    return ($hostGroup,$onlyVirtemp,$hostTaxid,$virTaxid);
}

################################################################################

sub PrintCommandError{
    my $errorMessage = $_[0];
    my $usage = $_[1];
    print "\n". $errorMessage . "\n";
    print "\n". $usage . "\n";
    exit;
}
