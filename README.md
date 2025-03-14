# Hückel Method Matrix Solver

A computational tool for building and diagonalizing Hückel matrices for molecular orbital calculations in quantum chemistry.

## Overview

This program creates and diagonalizes Hückel matrices based on input parameters, providing eigenvalues and eigenvectors that represent molecular orbital energies and coefficients. It can handle both linear (open) and cyclic molecular systems with customizable alpha (Coulomb) and beta (resonance) integrals.

## Features

- Support for both linear and cyclic molecular structures
- Customizable alpha (diagonal) and beta (off-diagonal) parameters
- Handling of alternating alpha and beta values for heterogeneous systems
- Efficient diagonalization using LAPACK library
- Output of eigenvalues and eigenvectors to separate files

## Requirements

- GCC compiler
- LAPACK and LAPACKE libraries
- Standard C libraries (stdio.h, stdlib.h, math.h)
- Python with pandas and matplotlib (for visualization)

Installation

Clone this repository
Ensure LAPACK and LAPACKE are installed on your system
Compile the program using the included makefile:

```bash
makefile
```

## Installation

1. Clone this repository
2. Ensure LAPACK and LAPACKE are installed on your system
3. Compile the program using the included makefile:

```bash
make
```

This will create the executable named `HUCKEL`.

## Usage

1. Create an `input.txt` file with the following format:

```
ATOM_NUMBER [integer]
ALPHA [value1] [value2]
BETA [value1] [value2]
OPEN [0 or 1]
```

Where:
- `ATOM_NUMBER`: Number of atoms/centers in the system
- `ALPHA`: Coulomb integral values (one or two values)
- `BETA`: Resonance integral values (one or two values)
- `OPEN`: System type (0 for cyclic, 1 for open/linear)

2. Run the program:

```bash
./HUCKEL
```

3. The program will generate two output files:
   - `eigenvalues_.txt`: Contains the eigenvalues (orbital energies)
   - `eigenvectors.txt`: Contains the eigenvectors (orbital coefficients)

## Input File Format

The program can handle different input formats, including:
- One or two alpha values
- One or two beta values

Example:
```
ATOM_NUMBER 1024
ALPHA 2 6
BETA -2 -5
OPEN 0
```

In this example:
- The system has 1024 atoms
- Alpha values alternate between 2 and 6
- Beta values alternate between -2 and -5
- The system is cyclic (OPEN = 0)

## Technical Details

### OUTPUT FOLDERS:

- Each folder in this repository represents the different permutations and combinations possible for the alpha and beta values which are scientifically relevant.
- Each of the folders also contain a `fileopen.py` which is basically the python program used to generate the plots of the eigenvalues to visualize the Fermi Levels.
- The images are saved to the corresponding folder as a `.png` file.

### Matrix Construction

- For linear systems (OPEN = 1): Only neighboring atoms are connected
- For cyclic systems (OPEN = 0): The first and last atoms are also connected
- When two alpha values are provided, they alternate along the diagonal
- When two beta values are provided, they alternate along the off-diagonals

### Diagonalization

The program uses LAPACK's `dgeev` function to calculate the eigenvalues and eigenvectors of the Hückel matrix. This provides an efficient solution for large matrices.

### Memory Management

The code includes proper memory management:
- Dynamic allocation for matrices and arrays [Not the best implementation but works decent]
- Custom 2D array allocation and deallocation functions [Taken from the Advanced Computational Techniques Course]
- Reference: https://github.com/scemama/tccm-homeworks/tree/master/project3
- Proper cleanup of all allocated resources

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Author

Copyright (c) 2025 Arya Advaith Satisha
