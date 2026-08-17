#!/usr/bin/perl
#by Dina Zhabinskaya
#This code is called by compete_many.thread.pl to run IR_finder and SIDD on single-transition code
#at input temperature for an input sequence file

use strict; use warnings;
use Getopt::Long;
use FindBin qw($RealBin);
my ($file,$trans, $out_file);
my $temp = 310;
my $sig = 0.06;
my $theta = 12;
my $salt = 0.01;
my $shape = "linear";
my $c = "";
my $b = "";
my $p = "";
my $r = "";
my $n = "";

my $usage = "\nusage: $0 -f <sequence file> -a <algorithm_type> (choose algorithm_type: M, Z, C, or A) [options]\n\n".
"Input requires a sequence file and algorithm type.\n".
"Script analyzes superhelically induced structural transition probabilities for each base pair.\n".
"Algorithm types: M (melting), Z (Z-DNA), C (cruciforms), and A (competition between all three).\n".
"Sequence file will be converted to the format required by the algorithm.\n".
"Code in directory trans_three/ will handle M, Z, and C algorithm types.\n".
"Code in directory trans_compete/ will handle A algorithm type.\n".
"For algorithm type options -a C and -a A user will need an Inverted Repeat Finder (IRF) executable compatible with user's operating system.\n".
"IRF download page: http://tandem.bu.edu/irf/irf.download.html.\n".
"Selected output will be printed to the screen.\n\n".
"Options:\n".
"-f       Required: specify sequence file \n".
"-a       Required: specify algorithm type: M, Z, C, or A\n".
"-t       Optional: set temperature (default 310)\n".
"-s       Optional: set superhelical density (default -0.06)\n".
"-i       Optional: set ionic strength (default 0.01)\n".
"-th      Optional: set energy threshold (default 12), -t 10 is recommended for -a A (competition algorithm)\n".
"-c       Optional: flag to set molecular type to circular (default linear)\n".
"-n       Optional: flag to set melting energetics to nearest neighbor (default copolymeric)\n".
"-b       Optional: flag to print base pair for each position (default null)\n".
"-p       Optional: flag to print algorithm parameters (default null)\n".
"-r       Optional: flag to print ensemble average results (default null)\n".
"-o       Optional: specify output file name\n";  


if($#ARGV < 2) {
	die $usage;
}

GetOptions (
"f=s"   => \$file,   # sequence file
"a=s"  => \$trans,   # algorithm
"T=s" => \$temp,     # temperature
"s=s" => \$sig,      # superhelical density
"i=s" => \$salt,     # ionic concentration
"th=s" => \$theta,   # energy threshold
"c" => \$c,          # molecular shape
"b" => \$b,          # print sequence
"n" => \$n,          # nearest neighbor energetics for melting
"p" => \$p,          # print parameters
"r" => \$r,          # print average results
"o=s" => \$out_file) # output file
or die("Error in command line arguments\n$usage\n");

$b="-b" if($b);
$p="-p" if($p);
$r="-r" if($r);
$n="-n" if($n);
if($c) {
    $c="-c";
    $shape="circular";
}


my @name = split("/",$file);
my $single_exe = "$RealBin/trans_three/qsidd";
my $compete_exe = "$RealBin/trans_compete/qsidd";

my $output_IR;
if ($trans eq "M") {
    #run single transition algorithm for melting
    $output_IR = `$single_exe $b $p $r $n $c -T $temp -s $sig -i $salt -t $theta -f $file`;
}
if ($trans eq "Z") {
    #run single transition algorithm for Z-DNA
    $output_IR = `$single_exe $b $p $r $c -T $temp -s $sig -i $salt -t $theta -Z -f $file`;
}
if ($trans eq "C" or $trans eq "A") {
    my $code_IR = "$RealBin/IR_finder.pl";
    my $IR_results = `perl $code_IR $temp $shape $file`; #run IR_finder
    if($trans eq "C") {
        #run single transition algorithm for cruciforms
        $output_IR = `$single_exe $b $p $r $c -T $temp -s $sig -i $salt -t $theta -C -X "$IR_results" -f $file`;
    }
    if ($trans eq "A") {
        #run competition algorith: melting, Z-DNA, and cruciforms competing
        $output_IR = `$compete_exe $b $p $r $n $c -T $temp -s $sig -i $salt -t $theta -X "$IR_results" -f $file`;
    }
}
if ($out_file) {
    open(my $out,">$out_file") || die "error creating $out_file\n";
    print $out $output_IR;
}
else {
     print  "$output_IR";
}
