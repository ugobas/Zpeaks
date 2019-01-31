###################### ZPeaks ################################

The program ZPeaks calls peaks from a genomic (NGS) experiment, comparing experiment and control over optimized windows.
Author: Ugo Bastolla <ubastolla@cbm.csic.es>
Centro de Biologia Molecular Severo Ochoa (CSIC-UAM), Madrid Spain

###############################################################


ZPeaks calls peaks from genomic experiments. It computes the Z score of experimental reads with respect to control reads (if the control is absent, averaged experiment reads are used) over windows.
The window size is optimized in such a way to maximize the discrimination power of the method.

Installation (linux):
====================

Download ZPeaks.zip to your computer and execute the following commands:

>mkdir <dir-name> (create a directory where ZPeaks will be stored)
>mv ZPeaks.zip <dir-name>
>cd <dir-name>
>unzip ZPeaks.zip
>make
>cp ZPeaks <your-path-directory>

Configuration file:
==================

You have to build a configuration file structured as the sample file "Input_ZPeaks_example.in" provided in the package.
The configuration file must specify the files that contain the experimental data in WIG format (EXPER parameter in configuration file), the data used for control (CONTROL parameter in configuration file; if the control is absent the program will use averaged experimental reads as a control), an optional BED format file containing the peaks that you want to compare with (PREDICTION parameter in configuration file) and profiles of genomic or epigenomic data measured across the chromosomes in GR or WIG formats (PROF parameter in configuration file).
All types of data can be contained in a single file or in a list of files, one for each chromosome. Both individual and groups of files must be terminated with the record END at the beginning of a line.
The name of the output file may be specified (NAME parameter in configuration file).
The parameters THR, SIZE_MIN, DCLUST and DTOL are the same for all experimental profiles.

GR format:      https://stackoverflow.com/questions/28880086/what-is-gr-file-format
BED format:     https://genome.ucsc.edu/FAQ/FAQformat.html#format1
WIG format:     https://genome.ucsc.edu/FAQ/FAQformat.html#format6

FORMAT of the configuration file:

THR=<Threshold for positives> [default=2.50]
DCLUST=<Distance threshold for joining fragments> [default=200] (Fragments above threshold that are at a bp distance < DCLUST are joined)
SIZE_MIN=<Minimum size for calling a peak> (Fragments smaller than SIZE_MIN are ignored, default=100)
DTOL=<Tolerance for comparison> [default=50] (Tolerance for identifying called peaks and reference peaks)
PEAKSIZE=<Peak size in bp> [default=150] (The property of called peaks are computed at +/-PEAKSIZE around the midpoint)
EXPER:
DIR=<path/to/exper/files>
<exper_file_1> (either one file or one file for each chromosome in WIG format, MANDATORY)

<exper_file_2> (OPTIONAL)

... 

<exper_file_n> (OPTIONAL)

END

CONTROL:

DIR=<path/to/control/files>

<control_file> (either one file or one file for each chromosome, OPTIONAL)

END

PREDICTION:

DIR=<path/to/reference/peaks>

<reference/peaks> (just one file, bed format OPTIONAL)

END

PROF1 <profile_name_1>

DIR=<path/to/prof_1>

<prof_1_file_1> (either single file or one for each chromosome, OPTIONAL)

<prof_1_file_2> (OPTIONAL)

... 

<prof_1_file_n> (OPTIONAL)

END

PROF2 <profile_name_2>

DIR=<path/to/prof_2>

<prof_2_file_1> (either single file or one for each chromosome, OPTIONAL)

<prof_2_file_2> (OPTIONAL)

... 

<prof_2_file_n> (OPTIONAL)

END

....

PROFn <profile_name_n>

DIR=<path/to/prof_n>

<prof_n_file_1> (either single file or one for each chromosome, OPTIONAL)

<prof_n_file_2> (OPTIONAL)

... 
<prof_n_file_n> (OPTIONAL)
END

Run:
===

>ZPeaks <config_file>

where config_file is the configuration file

Output:
======

The Zscore output comprises several files:

<name>_score.wig                            Zscore in wig format
<name>_ALL_NotCentered_<Parameters>.bed     Called peaks in bed format where <Parameters>=W<optimized_window_size>_T<THR>_J<DCLUST>_S<SIZEMIN>
<name>_ALL_<Parameters>.bed                 Called peaks in bed format, centered so that the mid point coincides with the maximum of the Z-score
Properties_<name>.dat                       Table in text format with the input genomic and epigenomic properties of every called peak, if these properties are provided

Optional output: 

If a reference set of peaks is input as PREDICTION, the program outputs the list of peaks that overlap with the reference (<name>_OLD_<Parameters>.bed),
that do not overlap with it (<name>_NEW_<Parameters>.bed) and the reference peaks that are not found in the current experiment (<name>_notfound_<Parameters>.bed)

