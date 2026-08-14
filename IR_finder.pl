#!/usr/bin/perl
#by Dina Zhabinskaya
#analyze a sequence with IR finder and output length of sequence, and position, length, and energy string used as input for SIDD-code for each IR
#this code can analyze circular plasmids in order to find possible IR's at the junction between the start and end position of the original sequence

use strict; use warnings;

my $usage = "\nusage: $0 <temperature> <shape> <sequence_file>\n\n".
"The script runs Inverted Repeat Finder (IRF), which can be downloaded at: http://tandem.bu.edu/irf/irf.download.html.\n".
"See below and description in reference (Nucleic Acids Res, 41(21), 9610) for IRF parameters used.\n".
"The output is a string of all possible IRs, with their start positions, lengths, and energies of extrusion.\n".
"Example of output: 1,0,10000,|2992,5,17.2738039658653,|2991,7,17.2738039658653,|2990,9,17.2738039658653,|...\n".
"Input requires temperature, shape of DNA, and sequence file:\n".
"temperature: numerical value in units of Kelvin\n".
"shape: circular or linear\n".
"sequence_file: provide a sequenece file\n".
"Sequence file (file_name) will be converted to a format appropriate for IRF: file one_line.file_name will be created.\n".
"An additional file circ.one_line.file_name will be created for circular sequences.\n";

if($#ARGV < 2) {
	die $usage;
}

my $code = "irf";

my $temp = $ARGV[0]; #temperature
my $shape = $ARGV[1]; #linear or circular
my $file = $ARGV[2]; #sequence fasta file
my @seq_name = split("/",$file);
my $in_file = $seq_name[-1];

use constant PI => 4*atan2(1,1);
my %bp = ('A','T','T','A','C','G','G','C');

#convert fasta file to fit IR_finder format
open(my $in,$in_file) || die "error opening $in_file\n";
my $seq_file = "one_line.$in_file";
open(my $out,">$seq_file") || die "error creating $seq_file\n";

my $first = 1;
while (my $line = <$in>) {
    if ($first) {
        if($line =~ m/>/) {
            chomp($line);
            print $out "$line\n";
        }
        else {
            print $out ">$seq_file\n";
        }
        $first = 0;
        next;
    }
    chomp($line);
    print $out "$line";
}
close($in);
close($out);

#IR finder parameters
my $match = 2;
my $mismatch = 10;
my $delta = 10;
my $pm = 80;
my $pi = 10;
my $minscore = 20;
my $maxlength = 10000;
my $maxloop = 100;
my $minloop = 3;

#energy parameters
my $salt = 0.01; #salt concentration
my $RT = 1.9872*$temp/1000.0;
my $Eb_c = 2.42; #imperfection energy
my $Ecr = 192.5-$temp*0.565-4*energy_melt("A")-2*2.44*$RT*log(4); #cruciform initiation energy

#set and initialize parameters
my @imper;
my @R = (); my @L = ();
my ($s_loop,$l_loop,$length_IR,$start_IR,$end_IR,$El);
$s_loop=$l_loop=$length_IR=$start_IR=$end_IR=$El=0;
my ($Et,$Eimp);
$Eimp=$Et=0;
my $shorten_arm;
my @IR_array;

my ($length_seq, $shift_seq);
my $circ_seq=0;
if ($shape eq "circular") {
    ($length_seq,$shift_seq)=get_info($seq_file);
    #create sequence file shifting the start position to the middle
    $circ_seq = "circ.$seq_file";
    open(my $out,">$circ_seq") || die "error creating $circ_seq\n";
    print $out ">$circ_seq\n";
    print $out "$shift_seq";
    close($out);
}

my $string_IR = "1,0,10000,|"; #format of string required for SIDD-code
$string_IR = get_result($seq_file,0);
if ($shape eq "circular") {
    $string_IR = get_result($circ_seq,1);
}
print "$string_IR\n";

#Functions:

#obtain a string containing IR start position, length, and energy as required for SIDD input
sub get_result {
    my ($in_file,$shape) = @_;
    my $skip = 0;
    #PARAMETERS: Match: $match Mismatch: $mismatch Delta: $delta pm: $pm pi: $pi minscore:$minscore maxlength: $maxlength maxloop: $maxloop \n";
    system("$code $in_file $match $mismatch $delta $pm $pi $minscore $maxlength $maxloop > /dev/null");
    my $num_files = `ls -l $in_file.$match.$mismatch.$delta.$pm.$pi.$minscore.$maxlength.$maxloop.*.txt.html | wc -l`;
    for(my $i = 1; $i <= $num_files; $i++) {
        my $out_file = "$in_file.$match.$mismatch.$delta.$pm.$pi.$minscore.$maxlength.$maxloop.$i.txt.html";
        open(my $output,$out_file) || die "error opening $out_file\n";
        my $read = 0;
        my @energies;
        while (my $line = <$output>) {
            $read = 0 if($line =~ m/frequent/);
            last if ($line =~m/Done/); #no more IRs
            if ($line =~ m/Loop:/) {
                $skip = 0;
                my @par_IR = get_IR($line); 	
                $l_loop = $par_IR[0]; #Loop length	
                next if ($l_loop > $maxloop);  #ignore loops of in this range
                $s_loop = $par_IR[1]; #Start position of loop 
                $start_IR = $par_IR[2]; #start of IR
                $end_IR = $par_IR[3]; #end of IR
                $length_IR = $par_IR[4]; #length of IR
                $IR_array[$s_loop][$l_loop][$start_IR]=$length_IR; #make array of positions and lengths
                $read = 1;
            }
            if ($read) {
                if ($line =~m/>>/ and $line !~m/LF/) { #Left arm
                    @L = get_arm($line,\@L);
                }
                if ($line =~m/<</ and $line !~m/RF/) { #Right arm
                    @R = get_arm($line,\@R);
                }
                if($line =~m/Statistics/) {
                    open(my $in,$in_file) || die "error opening $in_file\n";
                    my @l_seq = get_sequence($in,$s_loop,$l_loop);	#loop bp
                    $shorten_arm = 0;
                    my $val = -1;
                    #increase loop size to satify minloop requirement
                    while ($l_loop < $minloop) {
                        if ($L[$val] =~ m/A|T|C|G/) {
                            if ($R[$val] eq "-") {
                                push(@l_seq, $L[$val]);
                            }
                            elsif ($R[$val] eq "*") {
                                push(@l_seq, ($L[$val],$bp{$L[$val]}));
                            }
                            else {
                                push(@l_seq, ($L[$val],$R[$val]));
                            }
                        }
                        else {
                            push(@l_seq, $R[$val]);
                            $shorten_arm++;
                        }
                        $l_loop = @l_seq;
                        $s_loop-=1;
                        $shorten_arm++;
                        $val--;
                    }
                    #when circular redefine coordinates
                    if ($shape==1) {
                        if ($s_loop>int($length_seq/2)) {
                            $s_loop-=int($length_seq/2)+$length_seq%2;;
                            $start_IR-=int($length_seq/2)+$length_seq%2;;
                            
                        }
                        else {
                            $s_loop+=int($length_seq/2);
                            $start_IR+=int($length_seq/2);
                        }
                        if ($IR_array[$s_loop][$l_loop][$start_IR]) {
                            if ($IR_array[$s_loop][$l_loop][$start_IR]==$length_IR) {
                                $skip = 1;
                                @L = (); #initialize arm arrays
                                @R = ();
                            }
                        }
                    }
                    next if ($skip); #for circular sequences
                    $El = energy_loop(@l_seq);	#energy of loop
                    my $length_arm = @L - $shorten_arm;
                    @L = splice(@L,0,$length_arm);
                    @R = splice(@R,0,$length_arm);
                    @energies = energy_imperfections(\@R,\@L);
                    @energies = reverse(@energies);
                    $Et = $Ecr + $El;  #energy of imperfections
                    my $s_IR = $s_loop;
                    my $l_IR = $l_loop;
                    my $j=0;
                    for(my $i = 0; $i < $length_arm; $i++) {
                        $Et += 2*$energies[$i];
                        if ($L[$length_arm-1-$i] ne "-") {
                            $s_IR = $s_loop - ($j+1);
                            $l_IR = $l_loop + 2*($j+1);
                            $j++;
                        }
                        #string used as input in SIDD
                        $string_IR = join("",$string_IR,$s_IR,",",$l_IR,",",$Et,",","|");
                    }
                    @L = (); #initialize arm arrays
                    @R = ();
                }
            }
        }
        close($output);
    }
    return $string_IR;
}

#get length of sequence and shifted sequence when circular
sub get_info {
	my $seq_file = $_[0];
    open(my $in,$seq_file) || die "error opening $seq_file\n";
	my $l_seq;
    my $sf=0;
	while (my $line = <$in>) {
		next if($line =~ m/>/); #first line of fasta file
		$l_seq = length($line);
        if ($shape eq "circular") {
            my $s1 = substr($line,0,int($l_seq/2));
            my $s2 = substr($line,int($l_seq/2));
            #this generates a string where the start position of the sequence is shifted to the middle
            $sf = join("",$s2,$s1);
        }
	}
	return ($l_seq,$sf);
    close($in);
}

#get loop bp
sub get_sequence {
	my $in = $_[0];
	my $start = $_[1];
	my $l_loop = $_[2];
	my @seq;
	my @loop;
	while (my $line = <$in>) {
		next if($line =~ m/>/);
		@seq = split("",$line);
	}
	for(my $i = 0; $i < $l_loop; $i++) {
		$loop[$i]=$seq[$start+$i];
	}
	return @loop;
}

# convert an IR arm into an array
sub get_arm {
	my ($line,$A) = @_;
	my @results = split(" ",$line);
	my $arm = $results[2];
	my @left = split("",$arm);
    my $l_arm = @left;
    @left = splice(@left,0,$l_arm);
	push(@$A,@left);
	return @$A;
}

#obtain the length of arms and loop
sub get_IR {
	my $line = $_[0];
	my @results = split(" ",$line);
	my $l_loop = $results[3];	#lenght of loop
	my $arms = $results[1];
	my @rrl = split(",",$arms);
	my @la = split("--",$rrl[0]);
	my @ra = split("--",$rrl[1]);
	my $start = $la[0]; #start position of IR
	my $end = $ra[1]; #end of IR
	my $s_loop = $la[1]+1; #start position of loop
    my $l_IR = $end-$start+1; #length of IR
	return ($l_loop,$s_loop,$start,$end,$l_IR);
}

#calculate loop free energy
sub energy_loop {
	my @seq = @_;
    my $length = @seq;
	my $El = 0;
	for(my $i = 0; $i < $length; $i++) {
		$El += energy_melt($seq[$i]);
    }
    $El+=2*2.44*$RT*log($length); #loop entropy term
    return $El;
}

#calculate imperfections energy
sub energy_imperfections {
	my ($RA,$LA) =@_;
	my @R = @$RA;
	my @L = @$LA;
    my @energies;
    for(my $i = 0; $i < @L; $i++) {
        $energies[$i] = 0;
    }
	my $Eb = 0;
	my $E_miss = 0;
	for(my $i = 0; $i < @L; $i++) { #bulges on right arm
		if ($L[$i] eq "-") {
            my $melt = energy_melt($R[$i]);
            $energies[$i] = $Eb_c +$melt;
		}
	}
	for(my $i = 0; $i < @R; $i++) {
		if ($i == 0 and $R[$i] ne '*') { #test
			die "Doesn't start with a star\n";
		}
        if ($R[$i] eq "-") { #bulges on left arm
            my $melt = energy_melt($L[$i]);
			$energies[$i] = $Eb_c + $melt;
		}
        if ($R[$i] ne '*' and $R[$i] ne '-') {
            next if($L[$i] eq '-');  #case accounted for above
            if ($R[$i-1] ne '*' or $R[$i+1] ne '*') { #internal loop with 2bp open
                $energies[$i] = 2*$Eb_c + energy_melt($R[$i]) + energy_melt($L[$i]);
            }
            # missmatch
            if ($R[$i-1] eq '*' and $R[$i+1] eq '*' and $L[$i-1] ne '-' and $L[$i+1] ne '-') {
				my $tri1 = join("",$L[$i-1],$L[$i],$L[$i+1]);
				my $tri2 = join("",$bp{$L[$i-1]},$R[$i],$bp{$L[$i+1]});
				$energies[$i]= energy_missmatch($tri1,$tri2);
            }
        }
    }
	return @energies;
}

#energy of melting each base pair
sub energy_melt {
	my $base = $_[0];
	my $Ta = 354.65 + 16.6*log($salt)/log(10);
	my $Tg = $Ta + 41;
	my $Ea = 7.2464*(1-$temp/$Ta);
	my $Eg = 9.0172*(1-$temp/$Tg);
	if ($base eq 'a' or $base eq 't'or $base eq 'A' or $base eq 'T') {
		return $Ea;
	}
	else {
		return $Eg;
	}
}

#calculate the energy penalty of a missmatch
sub energy_missmatch {
	my ($tri1,$tri2) = @_;
	my @s1 = split("",$tri1);
	my @s2 = split("",$tri2);
	my $d1a = join('',$s1[0],$s1[1]);
	my $d1b = join('',$s1[1],$s1[2]);
	my %WC = energy_WC();
	my $e_WC=$WC{$d1a}+$WC{$d1b};
    my $d2a = join('',$s2[0],$s2[1]);
	my $d1r = join('',$s1[2],$s1[1]);
	my $d2r = join('',$s2[2],$s2[1]);
	my $m1 = join('.',$d1a,$d2a);
	my $m2 = join('.',$d2r,$d1r);
	my %miss = energy_miss();
	my $e_miss=$miss{$m1}+$miss{$m2};
	my $E_miss = $e_miss-$e_WC;
	return $E_miss;
}

#Watson-Crick base-specific energies (Table 1 in paper)
sub energy_WC {
	my %WC = ("AA",-1.0,"TT",-1.0,"AT",-0.88,"TA",-0.58,
    "CA",-1.45,"AC",-1.45,"GT",-1.44,"TG",-1.44,
    "CT",-1.28,"TC",-1.28,"GA",-1.3,"AG",-1.3,
    "CG",-2.17,"GC",-2.24,"GG",-1.84,"CC",-1.84
    );
	return %WC;
}

#Missmatch base-specific energies (Table 2 in paper)
sub energy_miss {
	my %miss = ("GA.CA",0.17,"GA.CC",0.81,"GA.CG",-0.25,
    "GC.CA",0.47,"GC.CC",0.79,"GC.CT",0.62,
    "GG.CA",-0.52,"GG.CG",-1.11,"GG.CT",0.08,
    "GT.CC",0.98,"GT.CG",-0.59,"GT.CT",0.45,
    "CA.GA",0.43,"CA.GC",0.75,"CA.GG",0.03,
    "CC.GA",0.79,"CC.GC",0.70,"CC.GT",0.62,
    "CG.GA",0.11,"CG.GG",-0.11,"CG.GT",-0.47,
    "CT.GC",0.40,"CT.GG",-0.32,"CT.GT",-0.12,
    "AA.TA",0.61,"AA.TC",0.88,"AA.TG",0.14,
    "AC.TA",0.77,"AC.TC",1.33,"AC.TT",0.64,
    "AG.TA",0.02,"AG.TG",-0.13,"AG.TT",0.71,
    "AT.TC",0.73,"AT.TG",0.07,"AT.TT",0.69,
    "TA.AA",0.69,"TA.AC",0.92,"TA.AG",0.42,
    "TC.AA",1.33,"TC.AC",1.05,"TC.AT",0.97,
    "TG.AA",0.74,"TG.AG",0.44,"TG.AT",0.43,
    "TT.AC",0.75,"TT.AG",0.34,"TT.AT",0.68,
	);
	return %miss;
}
