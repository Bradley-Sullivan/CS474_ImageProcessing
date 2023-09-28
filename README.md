# CS474 Image Processing Project Repository

## Overview

This is a repository for CS474/674 Image Processing programming assignments.

### Authors

Mohammed Ibrahimm (CS674) (m.s.ibrahim@nevada.unr.edu)

Bradley Sullivan (CS474) (bradleysullivan@nevada.unr.edu)

### Description

Here you will find MATLAB and/or C implementations of common to advanced image processing techniques. Specific assignments are separated into their own folder-heirarchy, composed of subdirectories corresponding to each assignment part.
    
    /CS474_ImageProcessing
        /AS1
            /part1
                <C program driver>
                <MATLAB script>
            /part2
                ...
            /part3
            /part4
        /AS2
            /part1
            ...
        ...
        /images
            <relevant PGM images>
        /inc
            <relevant header files>
        /src
            <relevant source files>
        Makefile

## Compilation Instructions

Here we detail how to compile and/or execute each assignment part. 

Parts programmed in C will require:
- An installation of Make [https://www.gnu.org/software/make/]
- An installation of GCC [https://gcc.gnu.org/]

Parts written in MATLAB will require:
- MATLAB stuff idk...

### Parts in C

The parts of assignments programmed in C will need to be compiled from source before they can be executed.

To compile a specific part with a C-implementation, execute the following command from the root of the cloned or downloaded repository:

    make [assignment] [part]

Compiling a specific part as such, the compiled executable will be placed next to the driver program for that part. 

Specific assignments will follow the form `AS[1:n]` (i.e. `AS1`, `AS2`, etc.). 

Different parts will follow the form `part[1:n]` (i.e. `part1`, `part2`, etc.).

To execute the compiled program, navigate through the heirarchy depicted above to the desired part's folder. Inside the specific part's folder will contain the driver program (typically a `driver.c` file) and its associated executable (typically named `pgm`).

In order to run `pgm`, enter this command into your terminal from the specific part folder

    ./pgm <args>

NOTE: **Different programs for different assignments/parts may require user-supplied command-line arguments. First try executing without any arguments, if arguments are required a print statement will detail the correct program usage.**

This will run the compiled program. Any program output will be placed beside this executable. This includes any output images.

### Parts in MATLAB

TODO

