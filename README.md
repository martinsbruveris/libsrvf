libsrvf
=======

What is it?
-----------

A library for shape analysis of elastic curves, using the square root velocity framework. It incorporates several algorithms to find the optimal matching between curves and functional data, including dynamic programming, the exact algorithm for functional data described in the PhD thesis of Daniel Robinson and the exact matching algorithm from Lahiri, Robinson and Klassen.


Dependecies
-----------

If you want to use libsrvf for rotational alignment, define the symbol USE_GSL=1 during compilation. libsrvf depends on the GNU Scientific Library (GSL) for its SVD routine. See www.gnu.org/software/gsl/ for instructions on getting this set up.

If you want to run the unit tests, then you will need the Boost libraries. See www.boost.org for instructions on installing these.

To build the library from source, you will need CMake, available at http://www.cmake.org/

Installation
------------

To compile the library for use by matlab, do the following.

```bash
cd libsrvf
mkdir build
cd build
cmake ..
cd ../matlab
make
```

Usage
-----

See the files in  `matlab/demo` for examples of how to use the code.

Attribution
-----------

Copyright 2018   Martins Bruveris (martins.bruveris@brunel.ac.uk)

Copyright 2012   FSU Statistical Shape Analysis and Modeling Group

This library is a derivative work of the code from the libsrvf library, which is licenced GPLv3. This code therefore is also licenced under the terms of GPLv3.

Licence
-------

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

