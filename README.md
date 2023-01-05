# telescope
Convert [Themisto](https://github.com/algbio/Themisto)
pseudoalignments to [kallisto](https://github.com/pachterlab/kallisto)
pseudoalignments.

Telescope uses the [BitMagic](https://github.com/tlk00/BitMagic)
library.

# Installation
Either use the precompiled binary
* [Linux 64-bit binary](https://github.com/tmaklin/telescope/releases/download/v0.2.1/telescope-v0.2.1_linux_x86-64.tar.gz)
* [macOS 64-bit binary](https://github.com/tmaklin/telescope/releases/download/v0.2.1/telescope-v0.2.1_macOS_x86-64.tar.gz)

or follow the instructions for compiling telescope and libtelescope from source.

## Compiling from source
### Requirements
- C++11 compliant compiler.
- cmake
- zlib

### How-to
- Clone the repository, enter the directory and run
```
> mkdir build
> cd build
> cmake ..
> make
```
- This will compile the telescope executable in build/bin/ and the libtelescope library in build/lib/.

# Usage
## Themisto to kallisto
Convert a single pseudoalignment against 10 reference sequences to kallisto format
```
telescope --n-refs 10 -r pseudos_1.txt -o kallisto_out_folder
```
... paired-end reads (pseudoaligned separately) considering a read aligned when __either__ of the reads align
```
telescope --n-refs 10 -r pseudos_1.txt,pseudos_2.txt -o kallisto_out_folder --mode union
```
... considering a read aligned when __both__ of the reads align
```
telescope --n-refs 10 -r pseudos_1.txt,pseudos_2.txt -o kallisto_out_folder --mode intersection
```
... considering the reads separately (read ids in the second file will
be incremented by the number of reads in the first file)
```
telescope --n-refs 10 -r pseudos_1.txt,pseudos_2.txt -o kallisto_out_folder --mode unpaired
```
... go crazy and do the same with multiple files (the order of the
reads must match for the union and intersection options to make sense)
```
telescope --n-refs 10 -r pseudos_1.txt,pseudos_2.txt,pseudos_3.txt,pseudos_4.txt -o kallisto_out_folder --mode union
```

## Merge Themisto paired alignment files
Convert two pseudoalignments from paired-end reads to a single `pseudos.aln` file by intersecting the pseudoalignments
```
telescope --n-refs 10 -r pseudos_1.txt,pseudos_2.txt -o pseudos --format themisto
```

## Accepted options
telescope accepts the following flags
```
telescope -r <strand_1>,<strand_2> -o <output prefix> --n-refs <number of pseudoalignment targets>
-r              Themisto pseudoalignment(s)
-o	            Output file directory.
--n-refs        Number of reference sequences in the pseudoalignment.
--format	    Output format (kallisto or themisto, default: kallisto
--mode	        How to merge paired-end alignments (one of unpaired, union, intersection; default: unpaired)
--read-compact	Read alignments that have been compressed with alignment-writer (default: false).
--write-compact	Write themisto format alignments in alignment-writer compressed format (default: true).
--cin           Read the last alignment file from cin (default: false).
--silent	    Suppress status messages (default: false)
--help	        Print the help message.
```

# License
telescope is licensed under the GNU Lesser General Public License v2.1. 

## Dependencies
[BitMagic](https://github.com/tlk00/BitMagic) is licensed under the [Apache-2.0 license](https://opensource.org/licenses/Apache-2.0).
