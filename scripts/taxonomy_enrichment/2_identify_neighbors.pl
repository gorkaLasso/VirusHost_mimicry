#!/usr/bin/perl

=description
1. get the user-defined input directory
2. read list of representative templates
3. map the templates to their hfpds
4. read neighbor file
     ignore those that cant be mapped to uniprot
     save in a hash those with sas <= user-defined criteria
5. map the pdb neighbors to uniprots
6. map the uniprots to host groups
=cut

use strict;
use warnings;

################################################################################
#    GLOBALS
################################################################################

my $repTempFile = "representative_template_space.txt";
my $pdbGenome = "/home/gl2411/ifs/data/c2b2/bh_lab/shares/databases/hfpd/genomes/glc_pdb_1/";

#no Sequence Filtering In Skan
#my $seqFilterTag = "";
#my $pdb2uniSkan60File = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/Databases/Skan_summary/2_pdb_in_skan60_withUni.txt";
#my $uni2hostGroupFile = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/Databases/Skan_summary/5_unique_uniprots_inSKAN60_withTaxid_Group.txt";

#sequence filtering in Skan at 90%
#my $seqFilterTag = "_filtAT90";
#my $pdb2uniSkan60File = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/Databases/Skan_summary/filtered_at_90/2_pdb_toUni_filteredAT90.txt";
#my $uni2hostGroupFile = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/Databases/Skan_summary/filtered_at_90/5_unique_uniprots_inSKAN60_withTaxid_Group_filteredAT90.txt";

#sequence filtering in Skan at 70%
my $seqFilterTag = "_filtAT70";
my $pdb2uniSkan60File = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/Databases/Skan_summary/filtered_at_70/2_pdb_toUni_filteredAT70.txt";
my $uni2hostGroupFile = "/home/gl2411/ifs/data/c2b2/bh_lab/gl2411/Databases/Skan_summary/filtered_at_70/5_unique_uniprots_inSKAN60_withTaxid_Group_filteredAT70.txt";

################################################################################
#    MAIN
################################################################################

&Main;

################################################################################
#    FUNCTION
################################################################################

sub Main{
  my ($inputDir,$sasThresh) = &GetInputOptions(@ARGV);
  print "input dir: $inputDir | sasThresh: $sasThresh\n";
  my %repTemp = &ReadTemplates("$inputDir\/$repTempFile");
  print "representative templates: " . keys(%repTemp) . "\n";
  my %pdb2hfpd = &ReadMapList($pdbGenome . "/fasta/map_list");
  print "map list loaded, entries: " . keys(%pdb2hfpd) . "\n";
  
  my (%pdb2uniprot);
  &ReadPdb2Uniprot($pdb2uniSkan60File,\%pdb2uniprot);
  
  my (%nbrUniprot);
  my $ct = 0;
  my $tempCt = keys(%repTemp);
  foreach my $repTemp (sort keys %repTemp){
    $ct++; print "procesing representative template $ct / $tempCt\n" if ($ct%1000 == 0);
    my ($pdb,$pdbStart,$pdbEnd) = split(/\t/,$repTemp);
    die "cant map pdb to hfpd: $pdb\n" if (not defined $pdb2hfpd{$pdb});
    my $hfpd = $pdb2hfpd{$pdb};
    my $fullHfpd = $hfpd . ".d" . $pdbStart . "_" . $pdbEnd;
    die "$fullHfpd does not exist" if (!-e "$pdbGenome\/fasta/$fullHfpd");
    my $nbrFile = $pdbGenome . "Seqs/" . $fullHfpd . "/Nbr/" . $fullHfpd . ".neigh.compact";
    if (!-e $nbrFile){
      print STDERR "nbr file does not exists:\n$nbrFile\n";
      next;
    }
    &ReadNeigh($nbrFile,$sasThresh,\%pdb2uniprot,\%nbrUniprot);
  }
  print "structural neighbors (uniprots): " . keys(%nbrUniprot) . "\n";
  
  my %uni2tax = &ReadUniTaxonomy($uni2hostGroupFile);
  my $outFile = $inputDir . "/uniprot_neighbor_space_" . $sasThresh . $seqFilterTag . ".txt";
  open (OUTFILE,">",$outFile) or die "cant open $outFile\n";
  foreach my $uni (sort keys %nbrUniprot){
    die "cant find taxonomic info for uniprot $uni" if (not defined $uni2tax{$uni});
    my $sas = $nbrUniprot{$uni};
    my $taxid = $uni2tax{$uni}{taxid};
    my $os = $uni2tax{$uni}{os};
    my $taxGroup = $uni2tax{$uni}{taxGroup};
    print OUTFILE "$uni\t$taxid\t$os\t$taxGroup\t$sas\n";
  }
  close OUTFILE; 
}

################################################################################

sub ReadUniTaxonomy{
  my $file = $_[0];
  my %uni2tax;
  
  open (INFILE,"<",$file) or die "cant open $file\n";
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    my ($uni,$taxid,$os,$taxGroup) = split(/\t/,$s);
    $uni2tax{$uni}{taxid} = $taxid;
    $uni2tax{$uni}{os} = $os;
    $uni2tax{$uni}{taxGroup} = $taxGroup;
  }
  close INFILE;
  
  return %uni2tax;
}

################################################################################

sub ReadNeigh{
  my $file = $_[0];
  my $sasThresh = $_[1];
  my $pdb2uniRef = $_[2];
  my $nbrUniRef = $_[3];
  
  open (INFILE,"<",$file) or die "cant open $file\n";
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    my ($pdb,$psd,$rms,$sas) = split(/\t/,$s);
    next if (not defined ${$pdb2uniRef}{$pdb});
    next if ($sas >= $sasThresh);
    #next if ($psd > 0.6);
    my $uni = ${$pdb2uniRef}{$pdb};
    if (not defined ${$nbrUniRef}{$uni}){
      ${$nbrUniRef}{$uni} = $sas;
    }elsif (${$nbrUniRef}{$uni} > $sas){
      ${$nbrUniRef}{$uni} = $sas;
    }
  }
  close INFILE;
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

sub ReadMapList{
  my $mapList = $_[0];
  my %pdb2hfpd;
  
  open(INFILE,"<",$mapList) or die $!;
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    my ($hfpd,$pdb) = split(/\t/,$s);
    $pdb =~ s/\>//;
    $pdb2hfpd{$pdb} = $hfpd;
  }
  close INFILE;
  return %pdb2hfpd;
}

################################################################################

sub ReadTemplates{
  my $file = $_[0];
  my %repTemp;
  
  open(INFILE,"<",$file) or die "cant open $file\n";
  while(<INFILE>){
    my $s = $_;
    chomp($s);
    my @w = split(/\t/,$s);
    $repTemp{"$w[0]\t$w[1]\t$w[2]"} = 1;
  }
  close INFILE;
  return(%repTemp);
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
