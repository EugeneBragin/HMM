package Algorithm::ViterbiC;

use vars qw/$VERSION/;
$VERSION = '0.01';

use Data::Dumper;

use strict;
use warnings;
use constant TINY_FLOAT => 1e-308;

sub new {
  my $class = shift; 
  my $self = {@_};
  bless $self, $class;

  $self->{states} = [ 0 .. scalar(@{$self->{start}})-1 ]
    if defined $self->{start};

  return $self;
}


sub train {
  die "train isn't implemented";
}


sub start {
  my $self = shift;
  $self->{start} = shift if @_;

  ## TODO: check input consistency
  $self->{states} = [ 0 .. scalar(@{$self->{start}})-1 ];
  
  return $self->{start};
}


sub emission {
  my $self = shift;
  $self->{emission} = shift if @_;

  ## TODO: check input consistency
  $self->{states} = [ 0 .. scalar(@{$self->{emission}})-1 ];
  
  return $self->{emission};
}


sub transition {
  my $self = shift;
  $self->{transition} = shift if @_;
  
  return $self->{transition};
}


sub forward_viterbi {
  my ($self, $o) = @_;

  my @V;
  my $path;

  # Initialize start cases (t == 0)
  for my $state (@{$self->{states}}) {
    $V[0]->{$state} = log($self->{start}[$state]) + log($self->get_emission($o->[0], $state));
    $path->{$state} = [ $state ];
  }
  
  my $t;
  for ($t = 1; $t < scalar(@$o); $t++) {
    print(($t+1) ."\t". $o->[$t]->{BAF} ."\t". $o->[$t]->{LRR} ."\t");
    my $newpath;
    
    my $maxE = -10000000000000000;
    my $maxEstate = undef;
    for my $y (@{$self->{states}}) {
      my $E = $self->get_emission($o->[$t], $y);
      if ($E > $maxE) {
        $maxE = $E;
        $maxEstate = $y;
      }
  
      $E = TINY_FLOAT if $E < TINY_FLOAT;
      
      my $max = undef;
      my $state;
      for my $y0 (@{$self->{states}}) {
        my $T = $self->get_transition($y0, $y, $o->[$t-1], $o->[$t]);
        $T = TINY_FLOAT if $T < TINY_FLOAT;

        die "T isn't defined" unless defined $T;
        die "E isn't defined" unless defined $E;
        die "V[".($t-1)."]prev p for state $y0 isn't defined" unless defined $V[$t-1]->{$y0};

        my $p = $V[$t-1]->{$y0} + log($T) + log($E);
        if (!defined($max) || $p > $max) {
          $max = $p;
          $state = $y0;
        }
      }
      
      $V[$t]->{$y} = $max;
      $newpath->{$y}  = [ @{ $path->{$state} }, $y ];
    }

    $path = $newpath;

    for my $y (@{$self->{states}}) {
      print "\t" . $V[$t]->{$y};
    }
    print "\n";
  }

  my $max = 0;
  my $state;
  $t--;
  
  for my $y (@{$self->{states}}) {
    if ($V[$t]->{$y} > $max) {
      $max   = $V[$t]->{$y};
      $state = $y;
    }
  } 

  return ($path->{$state}, $max);
}


sub get_emission {
  my ($self, $output, $state) = @_;

  my $p = 1;
  my $e = $self->{emission}[$state];  
  
  foreach my $key (keys %$e) {
    $p *= $e->{$key}($output->{$key}, $output)
      if defined($e->{$key}) && ref($e->{$key}) eq 'CODE';
  }

  $p = TINY_FLOAT if $p < TINY_FLOAT;

  ## UF adjustment (like in PennCNV)  
  return 0.01+0.99*$p;
}


sub get_transition {
  my ($self, $state1, $state2, $o1, $o2) = @_;
  
  my $p = $self->{transition}($state1, $state2, $o1, $o2);

  $p = TINY_FLOAT if $p < TINY_FLOAT;
  
  ## UF adjustment (like in PennCNV)  
  return 0.01+0.99*$p;
}


1;
