#!/usr/bin/perl -w

use strict;

# opens a generic .mdp file and replaces the values of
# 'init_lambda_state' with values at increments of 1

# If command-line arguments are not given, exit program
unless (@ARGV) {
    die "Usage: $0 input.mdp\n";
}

# define mdp file as first command-line argument
# $ is a scalar variable i.e. one element
my $mdp = $ARGV[0];

my $n_lam = $ARGV[1];

# split mdp file name by fullstop
# @ means array
my @temp = split('\.', $mdp);

# get the filename without file extension
my $base = $temp[0];

open(IN, "<$mdp>");
my $in = <IN>;
close(IN);

for (my $i=0; $i<$n_lam; $i++) {
	
	my $filename = "${base}_${i}.mdp";
	
}

exit;
