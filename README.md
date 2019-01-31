# ZPeaks
The program ZPeaks calls peaks from a genomic experiment, comparing experiment and control over optimized windows.
Author: Ugo Bastolla <ubastolla@cbm.csic.es>
Centro de Biologia Molecular Severo Ochoa (CSIC-UAM), Madrid Spain

ZPeaks calls peaks from genomic experiments. It computes the Z score of experimental reads with respect to control reads (if the control is absent, averaged experiment reads are used) over windows. The window size is optimized in such a way to maximize the discrimination power of the method.

Installation (linux):
>mv ZPeaks.zip <dir-name>
>cd <dir-name>
>unzip ZPeaks.zip
>make
>cp ZPeaks <your-path-directory>

Run:
>ZPeaks Input_ZPeaks
where Input_ZPeaks is a configuration file structured as the sample file Input_ZPeaks_10days_F4.in provided in the package

The configuration file must specify the files that contain the experimental data in wig format (EXPER: either a single file or one file for each chromosome, mandatory), the data used for control (CONTROL: optional; if the control is absent the program will use averaged eperimental reads as a control), an optional bed file containing the peaks that you want to compare with (PREDICTION) and profiles of genomic or epigenomic data measured across the chromosomes (PROF <name>).
All types of data can be contained in a single file or in a list of files, one for each chromosome. Both individual and groups of files must be terminated with the record END at the beginning of a line.
The name of the output file may be specified (NAME:)

The program can use either computes Z-scores of the ratio between the normalized experiment and the normalized control over windows of fixed size. The size of the window is optimized by maximizing the number of detected peaks that are above a given threshold.
THR=3.0	    ! Threshold
DCLUST=160  ! Fragments above threshold are joined for d<DCLUST bases
SIZEMIN=160 ! Fragments smaller than this are ignored
DTOL=50	    ! Tolerance for identifying called peaks and reference peaks
PEAKSIZE=150 ! The property of called peaks are computed at +/-PEAKSIZE around the midpoint.

Output: 

<name>_score.wig
(Zscore in wig format)
<name>_ALL_NotCentered_<Parameters>.bed
(called peaks in bed format)
where <Parameters>=W<optimized_window_size>_T<THR>_J<DCLUST>_S<SIZEMIN>
<name>_ALL_<Parameters>.bed
(called peaks in bed format, centered so that the mid point coincides with the maximum of the Z-score)
Properties_<name>.dat
(Table in text format with the input genomic and epigenomic properties of every called peak, if these properties are provided)
If a reference set of peaks is input as PREDICTION, the program outputs the list of peaks that overlap with the reference (<name>_OLD_<Parameters>.bed), that do not overlap with it (<name>_NEW_<Parameters>.bed) and the reference peaks that are not found in the current experiment (<name>_notfound_<Parameters>.bed)

FORMAT of the configuration file:

EXPER:
DIR=<path of exper files>
<exper file> (either one file or one file for each chromosome, MANDATORY)
END
CONTROL:
DIR=<path of control files>
<control file> (either one file or one file for each chromosome, optional)
END
PREDICTION:
DIR=<path of prediction file>
<prediction file> (only one, optional)
END
PROF <profile name>:
DIR=<path of prof results>
<prof file> (one for each chromosome, optional)
END
(Same for all experimental profiles)
THR=<Threshold for positives> (default 2.00)
SIZE_MIN=<Minimum size for calling a peak>
DCLUST=<Distance threshold for joining fragments>
DTOL=<Tolerance for comparison>

