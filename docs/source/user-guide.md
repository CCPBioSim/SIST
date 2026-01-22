# SIST: stress-induced structural transitions in superhelical DNA 

## Purpose of SIST

The codes in this repository are for analyzing three types of structural transitions in superhelical DNA molecules of specified base sequences and kilobase lengths.  These are strand separation, BZ transitions and cruciform extrusion.  More types of transitions may be added as their energetics become known.  The statistical mechanical methods and algorithms used in these analyses are described in the papers cited below. 

## Codes

1. *master.pl*: This is a pipeline to run all algorithms.

2. *IR_finder.pl*: analyzes the output of Inverted Repeat Finder (see below), and outputs a list of inverted repeats and their start positions, all their possible extrusion lengths, and cruciform formation energies for each such length.

3. *trans_three*: C++ code containing the algorithms for analyzing strand separation,  Z-DNA, and cruciform extrusion independently.

4. *trans_compete*: C++ code containing the algorithm for analyzing the competition between strand separation, Z-DNA, and cruciform extrusion.

## Instructions for running *master.pl*

The code *master.pl* contains a pipeline that allows the user to run either trans_three or trans_compete with various parameters provided in the menu of the code. Run “perl master.pl” for more detailed instructions and to see all available parameters.  Sequence file fasta format is recommended. If another format is provided, the sequence is converted into an appropriate format.  Multiple sequences in one file are not permitted.

**Algorithm types:**

1. -a M: melting transition only (SIDD).

2. -a Z: Z-DNA transition only.

3. -a C: cruciform transition only.

4. -a A: competition between melting, Z-DNA, and cruciform transitions.


**To compile and run follow these steps:**

1. Download *master.pl*, *IR_finder.pl*, trans_three/ and trans_compete/ folders into a working directory.

2. Compile the C++ codes: go to both trans_three/ and trans_compete/ directories, and type **make** on the command line.

3. Download Inverted Repeats Finder (IRF) that suits your operating system and move it to your working directory. This code can be downloaded here: https://github.com/Benson-Genomics-Lab/IRF/releases/tag/IRFv3.09.

4. Ignore this step if you are a Linux user.  Otherwise, replace the name of IRF executable in the line “my $code = "irf308.linux.exe;" in IR_finder.pl with the executable appropriate for your operating system.

5. Move a sequence file you want to analyze into the working directory.

6. To execute the melting transition (SIDD) with default parameters, type:  perl master.pl -a M -f sequence_file.

**More examples of step 6 above:**

Competition with default parameters and a specified output file: 


```ruby
perl master.pl -a A -f sequence_file -o output_file
```

Competition code at superhelix density of σ = -0.07 for a circular plasmid, displaying all available output:

```ruby
perl master.pl -a A -s 0.07 -c -b -p -r -f sequence_file
```

## Usage of C++ codes

Type make on command line to compile.  Once this is done, qsidd is the executable. 
Run “./qsidd” for more detailed instructions and to see all available parameters. To run with default parameters: ./qsidd -f sequence_file. 

The following steps are to run cruciform analysis with trans_three and the competition analysis with trans_compete:

1. First run *IR_finder.pl*:
```ruby
perl IR_finder.pl temperature shape sequence_file
```

2. For cruciforms using trans_three run: 
```ruby
./qsidd -C -X “string” -f sequence_file
```

3. For competition using trans_compete run: 
```ruby
./qsidd -X “string” -f sequence_file
```

Here “string” is the output from *IR_finder.pl*.  Run “ perl IR_finder.pl” for a more detailed explanation of step 1.  Note: *master.pl* can do the above workflow automatically.

## Algorithm limitations and recommended parameter ranges

If the codes are executed outside of the ranges listed below, a warning will be printed in the output.  We advise that you do not perform the analysis outside the recommended guidelines. Master.pl options corresponding to each parameter are listed below.


> [!WARNING]
> Algorithm warnings

1. Sequence length in input sequence file should be greater than 1,500 and less than 10,000.

2. Energy Threshold (-th) below 9 many yield inaccurate results, and -th above 15 many results in very long execution times.

3. Absolute value of superhelical density (-s) greater than 0.15 may be outside of physiological range.

4. Temperatures (-t) less than 220 K and greater than 320 K may be outside physiological range.

5. Salt concentration less than 0.0001 M may be outside physiological range.

## Citations

When using these algorithms you must cite the first paper below, and some or all of the others, depending on which types of analyses you perform:

Zhabinskaya, D., Madden, S., & Benham, C. J. (2015). SIST: stress-induced structural transitions in superhelical DNA. Bioinformatics, 31(3), 421-422.

Fye, R. M. and Benham, C. J. (1999), “Exact method for numerically analyzing a model of local denaturation in superhelically stressed DNA”, Phys Rev E, 59, 3408-3426.

Zhabinskaya, D. and Benham, C. J. (2011), “Theoretical Analysis of the Stress Induced BZ Transition in Superhelical DNA”, PLoS Comput Biol, 7, 1-14.

Zhabinskaya, D. and Benham, C. J. (2013), “Competitive superhelical transitions involving cruciform extrusion”, Nucleic Acids Res, 41(21), 9610-9621.

## Contact Information

For questions or problems, please contact either Dina Zhabinskaya (dzhabinskaya@ucdavis.edu), Craig Benham (cjbenham@ucdavis.edu), or Sally Madden (sallymadden@gmail.com).
