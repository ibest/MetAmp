
#This script extracts variable regions from 16S rRNA. 

##Usage

```
extract.py -i <input file path> -o <output file prefix> [options]

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --input INFILE
                        input file
  -o OUTPREFIX, --output OUTPREFIX
                        output prefix
  -fp FORWARD_PRIMER, --forward_primer FORWARD_PRIMER
                        forward primer sequence
  -rp REVERSE_PRIMER, --reverse_primer REVERSE_PRIMER
                        reverse primer sequence
  -m MATCH_AWARD, --match_award MATCH_AWARD
                        match award
  -p MISMATCH_PENALTY, --mismatch_penalty MISMATCH_PENALTY
                        mismatch penalty
  -g GAP_PENALTY, --gap GAP_PENALTY
                        gap penalty
  -pos, --positional    positional trimming
  -spos SPOS, --start_position SPOS
                        positional trimming: from
  -epos EPOS, --end_position EPOS
                        positional trimming: to
  -nr NUM_READS, --num_reads NUM_READS
                        number of reads to extract from an input file
  -v, --verbose         verbosing output
  -ver, --version       print version and exit
```

##Examples

###If you have primers:

`python extract.py -i gold.fa -o test -fp CCTACGGGAGGCAGCAG -rp CCGTCAATTCMTTTRAGN`

Here: 
gold.fa - an input file that contains 16S sequences, must be provided in FASTA format;
test - an output prefix for file that contains extracted variable regions, will be in FASTA format (i.e. test.fasta);
CCTACGGGAGGCAGCAG - forward primer
CCGTCAATTCMTTTRAGN - reverse primer


### Some primer sequences that can be used

#### V1-3 primers:

```
forward_primer = "NGAGTTTGATCCTGGCTCAG" # -m=2 -p=-5 g=-8
reverse_primer = "ATTACCGCGGCTGCTGG"
```
#### V3-5 primers:

```
forward_primer = "CCTACGGGAGGCAGCAG" # -m=2 -p=-5 g=-8
reverse_primer = "CCGTCAATTCMTTTRAGN"
```

#### V6-9 primers:

```
forward_primer = "GAATTGACGGGGRCCC" # -m=2 -p=-5 -g=-1
reverse_primer = "TACGGYTACCTTGTTAYGACTT"
```

###If you have positions (positional trimming):

`python extract.py -i gold.fa -o test -pos -spos 0 -epos 100`

###How to use Shiny version of Extracter:

* Open ui.R with for example, RStudio and then click on "Run App" button.
* Provide your 16S records, set up parameters if necessary and then click on "Start computation" button.

Table on the right will show the final statistics. 

##General suggestions:

- Change the alignment parameters (match award, mismatch penalty, open gap penalty) if average length is too small (less then 50 bp or even zero).

Last change was made on 2016/08/10
Author: Ilya Y. Zhbannikov
