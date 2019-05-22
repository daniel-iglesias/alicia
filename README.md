# ALICIA
Augmented Lagrangian Implicit Constrained Inverse Analysis

## Description
ALICIA is a code for calculating heat flux density at tokamak PFCs. It is mainly developed for increased accuracy in the measurement of the SOL normal heat flux density at irregular divertor tiles using tempertures acquired by fast IR cameras.

ALICIA uses MkniX library for generating and solving the 2D numerical FEM models, and is distributed under the GPL v2 license (see License section below).

## Installation

ALICIA uses the cmake build system for cross-platform compatibility (version 2.8 or greater).

Before building ALICIA, the source code and shared library of MkniX shall be available locally. Please use Mknix devel branch.

Build ALICIA using (from the ALICIA root directory):

```
cmake -B<build directory> -H. -DCMAKE_BUILD_TYPE=<Debug|Release> -DMKNIX_LINK_DIR=<path/to/dir/containing-libmknix.so> -DMKNIX_SOURCE_DIR=<path/to/dir/containing/MkniX>
```

The use of cmake-gui is recommended. Alternatively, an example of the cmake one-liner command is provided (note that paths might differ for each user):

```
cmake -B../ALICIA-build -H. -DCMAKE_BUILD_TYPE=Release -DMKNIX_LINK_DIR=../MkniX-devel-build/src/ -DMKNIX_SOURCE_DIR=../MkniX
```

## Examples

ALICIA is provided with some examples within the folder "ALICIA-examples" so it can be tested in the deployment platform.

To list all example directories:
```
cd <build directory>/ALICIA-examples
ls
```
## Usage

Once navigated to an specific example or working directory, ALICIA can be run using two command line arguments: (1) the MkniX input file, and (2) the IR temperatures file.

e.g.
```
cd <build directory>/ALICIA-examples/90271-TPT6-augmented-5deg-alpha1E3-res1E1/
../../src/alicia input_tile6_fem.mknix temps_IR.txt 
```

## License

Copyright (C) 2016 by Daniel Iglesias

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

