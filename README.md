# telescope
Convert [Themisto](https://github.com/jnalanko/Themisto)
pseudoalignments to [kallisto](https://github.com/pachterlab/kallisto) pseudoalignments.

# Installation
Either use the precompiled binary
* [Linux 64-bit binary](https://github.com/tmaklin/telescope/releases/download/v0.2.0/telescope-v0.2.0_linux_x86-64.tar.gz)
* [macOS 64-bit binary](https://github.com/tmaklin/telescope/releases/download/v0.2.0/telescope-v0.2.0_macOS_x86-64.tar.gz)

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
Convert a single pseudoalignment against 'themisto-index' to kallisto format
```
telescope --index themisto-index -r pseudos_1.txt -o kallisto_out_folder
```
... paired-end reads (pseudoaligned separately) considering a read aligned when __either__ of the reads align
```
telescope --index themisto-index -r pseudos_1.txt,pseudos_2.txt -o kallisto_out_folder --mode union
```
... considering a read aligned when __both__ of the reads align
```
telescope --index themisto-index -r pseudos_1.txt,pseudos_2.txt -o kallisto_out_folder --mode intersection
```
... considering the reads separately (read ids in the second file will
be incremented by the number of reads in the first file)
```
telescope --index themisto-index -r pseudos_1.txt,pseudos_2.txt -o kallisto_out_folder --mode unpaired
```
... go crazy and do the same with multiple files (the order of the
reads must match for the union and intersection options to make sense)
```
telescope --index themisto-index -r pseudos_1.txt,pseudos_2.txt,pseudos_3.txt,pseudos_4.txt -o kallisto_out_folder --mode union
```

## Kallisto to Themisto
not yet implemented...

# License
telescope is licensed under the GNU Lesser General Public License v2.1. 
