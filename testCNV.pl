
use strict;
use warnings;

use Algorithm::ViterbiC;
use Data::Dumper;
use Math::Trig;

## Probability density function for a normal distribution with a given mean and standard deviation
sub pdf_normal {
  my ($x, $mean, $sd) = @_;
  return exp( (-($x-$mean)**2) / (2*$sd**2) ) / ($sd * sqrt(2*pi));
}


## TODO
## 1. read PBF file, make position => PBF hash
## 2. assign PBF to every observation

my $pfb;
print "Reading PFB file ... ";
open FILE, "scotland-PFB.txt" or die $!;
my @lines = <FILE>;
close FILE;

foreach my $line (@lines) {
  chop $line;
  my @bits = split(/\t/, $line);
  $pfb->{$bits[0]} = $bits[3];
}
print "done\n";

undef @lines;


print "Reading observations file ... ";
open FILE, "FR-1.txt" or die $!;
my @lines = <FILE>;
close FILE;

my @headers = split(/\t/, substr(shift @lines, 0, -1));

my @observations;
for my $line (@lines) {
  chop $line;
  my $i = 0;
  my %hash = map { $headers[$i++] => $_ } split(/\t/, $line);

  push @observations, \%hash
    if scalar(keys %hash) == scalar(@headers) && $hash{Position} && $hash{LRR} =~ /\d+/ && $hash{BAF} =~ /\d+/;

  if ($hash{SNPName}) {
    if (defined $pfb->{$hash{SNPName}}) {
      $hash{pfb} = $pfb->{$hash{SNPName}};
    } else {
      $hash{pfb} = 0.5;
      warn "No PFB for $hash{SNPName}, setting as 0.5";
    }
  }
}

@observations = sort { $a->{Position} <=> $b->{Position} } @observations;

print "done\n";

my $E = [
  ## Two copies deletion
  {
    LRR => sub{ pdf_normal($_[0], -3.527211, 1.329152) },
    BAF => sub{ pdf_normal($_[0], 0.5, 0.304243) },
  },
  
  ## One copy deletion
  {
    LRR => sub{ pdf_normal($_[0], -0.664184, 0.284338) },
    BAF => sub{ (1 - $_[1]->{pfb}) * pdf_normal($_[0], 0, 0.016372) 
                    + $_[1]->{pfb} * pdf_normal($_[0], 1, 0.016372) },
  },
  
  ## Normal (2 copies)
  {
    LRR => sub{ pdf_normal($_[0], 0, 0.159645) },
    BAF => sub{    ((1-$_[1]->{pfb})**2) * pdf_normal($_[0], 0, 0.016372) 
                +  2 * (1-$_[1]->{pfb}) * $_[1]->{pfb} * pdf_normal($_[0], 0.5, 0.034982)
                +  ($_[1]->{pfb}**2) * pdf_normal($_[0], 1, 0.016372) },
  },
  
  ## Copy-neutral LOH
  {
    LRR => sub{ pdf_normal($_[0], 0, 0.211396) },
    BAF => sub{     (1 - $_[1]->{pfb}) * pdf_normal($_[0], 0, 0.016372)
                 +  $_[1]->{pfb} * pdf_normal($_[0], 1, 0.016372) },
  },
  
  ## Single copy duplication
  {
    LRR => sub{ pdf_normal($_[0], 0.395621, 0.209089) },
    BAF => sub{     ((1 - $_[1]->{pfb})**3) * pdf_normal($_[0], 0, 0.016372)
                 +  3 * ((1 - $_[1]->{pfb})**2) * $_[1]->{pfb} * pdf_normal($_[0], 0.333333, 0.045126)
                 +  3 * (1-$_[1]->{pfb}) * ($_[1]->{pfb}**2) * pdf_normal($_[0], 1-0.333333, 0.045126)
                 +  ($_[1]->{pfb}**3) * pdf_normal($_[0], 1, 0.016372) },
  },
  
  ## Double copy duplication
  {
    LRR => sub{ pdf_normal($_[0], 0.678345, 0.191579) },
    BAF => sub{    ((1 - $_[1]->{pfb})**4) * pdf_normal($_[0], 0, 0.016372)
                +  4 * ((1 - $_[1]->{pfb})**3) * $_[1]->{pfb} * pdf_normal($_[0], 0.250000, 0.042099)
                +  6 * ((1-$_[1]->{pfb})**2) * ($_[1]->{pfb}**2) * pdf_normal($_[0], 0.5, 0.034982)
                +  4 * (1-$_[1]->{pfb}) * ($_[1]->{pfb}**3) * pdf_normal($_[0], 1-0.250000, 0.042099)
                +  ($_[1]->{pfb}**4) * pdf_normal($_[0], 1, 0.016372) },
  },
];


my $TM = [
  [ 0.905850086, 0.000000001, 0.048770575, 0.045379330, 0.000000010, 0.000000003 ], 
  [ 0.000000001, 0.950479016, 0.048770575, 0.000750402, 0.000000007, 0.000000003 ],
  [ 0.000001064, 0.000024530, 0.998795591, 0.001165429, 0.000012479, 0.000000912 ],
  [ 0.000049998, 0.000049998, 0.000049998, 0.999793826, 0.000049998, 0.000006187 ],
  [ 0.000000001, 0.000000001, 0.048770575, 0.001248044, 0.949981383, 0.000000001 ],
  [ 0.000000001, 0.000000001, 0.017682158, 0.000000001, 0.000297693, 0.982020152 ],
];


my $T = sub {
  my ($state1, $state2, $observation1, $observation2) = @_;

  my $dist = $observation2->{Position} - $observation1->{Position};

  my $TMadj;

  for (my $i = 0; $i < scalar(@$TM); $i++) {
    my $offdiagonal_sum = 0;
    for (my $j = 0; $j < scalar(@$TM); $j++) {
      if ($i != $j) {
        if ($i == 3) {
          $TMadj->[$i][$j] = $TM->[$i][$j] * (1-exp(-$dist/100))  / (1-exp(-5000/100));
        } else {
          $TMadj->[$i][$j] = $TM->[$i][$j] * (1-exp(-$dist/100000)) / (1-exp(-5000/100000));
        }
        if ($TMadj->[$i][$j] > 1) {
          warn "WARNING: Off-diagonal cell TM[$i][$j] (". $TM->[$i][$j] ." to ". $TMadj->[$i][$j] ." by $dist) in transition matrix is over boundary of 1 (HMM model is not optimized). Assign 0.999 as the value instead.\n";
          $TMadj->[$i][$j] = 0.999;     ## maximum possible off-diagonal value
        }
        
        $offdiagonal_sum += $TMadj->[$i][$j];
      }
    }
    
    if ($offdiagonal_sum >= 1) {
      for (my $j = 0; $j < scalar(@$TM); $j++) {
        next if $i == $j;
        ##die "TMadj[$i][$j] isn't defined" unless defined $TMadj->[$i][$j];
        $TMadj->[$i][$j] /= ($offdiagonal_sum/0.999);
      }
      $offdiagonal_sum = 0.999;
    }
    $TMadj->[$i][$i] = 1 - $offdiagonal_sum;    
  }
  
  return $TMadj->[$state1][$state2];
};


my $START = [ 1e-9, 0.000500, 0.999000, 1e-9, 0.000500, 1e-9 ];


my $v = Algorithm::ViterbiC->new(
  start      => $START,
  emission   => $E,
  transition => $T,
);

my ($v_path, $v_prob) = $v->forward_viterbi(\@observations);

for (my $i = 0; $i < scalar(@observations); $i++) {
  $observations[$i]->{state} = $v_path->[$i];
}

open FILE, ">out.txt" or die $!;
print FILE $_->{state} ."\t". $_->{BAF} ."\t". $_->{LRR} ."\t". $_->{Position}. "\n" for @observations;
close FILE;
