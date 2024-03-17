# FORTRAN to Dart

## References
 - LAPACK constants: https://netlib.org/lapack/explore-html/d4/d69/namespacela__constants.html
 - Format specifiers: https://docs.oracle.com/cd/E19957-01/805-4939/z40007437a2e/index.html
 - Intrinsics
    - https://gcc.gnu.org/onlinedocs/gfortran/Intrinsic-Procedures.html
    - https://www.ibm.com/docs/en/xl-fortran-aix/16.1.0
 - Online compiler: https://onecompiler.com/fortran/



## Regex Helpers

### Quick
 - Cleanup -> ` // \d+$`
 - To Matrix -> `\b([XY])\(([^())]+),([^()]+)\)` -> `$1[$2][$3]`
 - To Array -> `\b([XY])\(([^())]+)\)` -> `$1[$2]`
 - Goto? -> `^\s+\} // \d+`

#### WARN:
 - Fix `Box.value` in strings -> `'[^']*\.value.*'`

### Line continuations
```
^\s     (.+)(\n\s    [^\s]\s*(.*))$
      $1 $3
```

### Format continuation lines
```
^([\s\d][\s\d][\s\d][\s\d][\s\d]\s+.+)(\n\s    [^\s]\s*(.*))$
$1 $3
```

### DOUBLE PRECISION -> double
```
^\s(.*)DOUBLE PRECISION(.*)$
 $1double          $2;
```

### INTEGER -> int
```
^\s(.*)\bINTEGER\b(.*)$
 $1int    $2;
```

### LOGICAL -> bool
```
^\s(.*)\bLOGICAL\b(.*)$
 $1bool   $2;
```

### CHARACTER*1 -> String
```
^\s(.*)\bCHARACTER\b\*[0-9]\s(.*)$
 $1String      $2;
```

### CHARACTER*(*) -> List<String>
```
^\s(.*)\bCHARACTER\b\s*\*\s*\(\s*\*\s*\)(.*)$
 $1List<String>   $2;
```

### CHARACTER -> String
```
^\s(.*)\bCHARACTER\b(.*)$
 $1String   $2;
```

### Comment out INTRINSIC / COMMON
```
^\s     (\bINTRINSIC\b)\s*(.*)$
      // $1 $2
```

### Comment out EXTERNAL
```
^     ( \bEXTERNAL\b.*)$
//    $1
```

### PARAMETER (...) -> const ...
```
^\s     \bPARAMETER\b(\s*)\(\s*(.*)\)\s*$
      const    $1$2;
```

### SUBROUTINE -> void
```
^\s     SUBROUTINE\s+(\w+)\s*\((.*)$
      void $1($2 {
```

### Comments -> aligned //
```sh
find . -type f -regex '.*\.[fF]' | xargs -I % sed -Ei '' 's/^[^ \t](     [ \t]*)([^ \t].*)$/\1 \/\/ \2/g' %
```

```
^( +)(//.*)\n\1( +)//
$1$2\n$1//$3
```

### Remove empty comment lines
```sh
find . -type f -regex '.*\.[fF]' | xargs -I % sed -Ei '' 's/^[^ \t]$//g' %
```

### Comment out the remaining comments
```
^[^\s\d#]\s*(.*)$
// $1
```

### Comment out the preprocessor directives
```
^(#\s*.*)$
// $1
```

### END -> }
```sh
find . -type f -regex '.*\.[fF]' | xargs -I % sed -Ei '' 's/^      ([ ]*)END$/      \1}/g' %
```

### Comma
```
^(\s+[^/\s}](.(?!(\{|\}|;)))*)$
$1;
```

Remove artifact
```
^(\s+\});+$
$1
```

```
^(\s+\} // \d+);$
$1
```

### ELSE -> } else {
```
^(\s+)ELSE$
$1} else {
```

### END IF -> }
```
^(\s+)END\s?IF$
$1}
```

### IF ... THEN -> if ... {
```
^(\s+)(IF)\s*(.*?)\s*THEN$
$1if $3 {
```

### ELSE IF ... THEN -> } else if ... {
```
^(\s+)(ELSE\s?IF)\s*(.*?)\s*THEN$
$1} else if $3 {
```

### No arguments CALL
```
^(\s+)CALL\s+(\w+)$
$1\L$2();
```

### CALL
```
^(\s+)CALL\s+(\w+)\s*\(\s*(.*)\s*\)$
$1\L$2($3);
```
### Matrix args
```
^(\s+\w+\s+\w+\((\s*\w+\s*,)+)\s*(\w+)\s*,\s*LD\3\s*,\s*(.+)\)\s*\{$
$1 final Matrix<double> $3, final int LD$3, $4) {
```

### Array 2x args (WORK, LWORK)
```
^(\s+\w+\s+\w+\((\s*((final\s+)\w+(<\w+>)?\s+)?\w+\s*,)+)\s*IWORK\s*,\s*LIWORK\s*,\s*(.+)\)\s*\{$
$1 final Array<double> WORK, final int LWORK, $6) {
```

### Array arg
```
^(\s+\w+\s+\w+\((\s*((final\s+)\w+(<\w+>)?\s+)?\w+\s*,)+)\s*DOTYPE\s*,\s*(.+)\)\s*\{$
$1 final Array<bool> DOTYPE, $6) {
```

### Last arg
```
^(\s+\w+\s+\w+\((\s*((final\s+)\w+(<\w+>)?\s+)?\w+\s*,)+)\s*\wWORK\s*\)\s*\{$
$1 final Array<double> WORK) {
```

### First arg
```
^(\s+\w+\s+\w+\(\s*)(\w+)\s*((,\s*((final\s+)\w+(<\w+>)?\s+)?\w+\s*)*\)\s*\{)$
$1final int $2$3
```

### Unknown args
```
^(\s+\w+\s+\w+\((\s*((final\s+)\w+(<\w+>)?\s+)?\w+\s*,)+)\s*(\w+)\s*,\s*(.+)\)\s*\{$
$1 final int $6, $7) {
```

### Trailing args comma
```
^(\s+\w+\s+\w+\(\s*((final\s+)\w+(<\w+>)?\s+)?\w+\s*(,\s*((final\s+)\w+(<\w+>)?\s+)?\w+\s*?)*)\s*\)\s*\{$
$1,) {
```

### DO
Simple
```
^(\s+)DO\s+(\d*)\s*(\w+)\s*=\s*(\w+)\s*,\s*(\w*)$
$1for ($3 = $4; $3 <= $5; $3++) { // $2
```

```
\) \{ //$
) {
```

Improved
```
^(\s+)DO\s+(\d*)\s*(\w+)\s*=\s*([^,]+)\s*,\s*([^,]*)$
$1for ($3 = $4; $3 <= $5; $3++) { // $2
```

With functions LR
```
^(\s+)DO\s+(\d*)\s*(\w+)\s*=\s*(\w+\s*\([^(]*\)),\s*(\w+\s*\([^(]*\));?$
$1for ($3 = $4; $3 <= $5; $3++) { // $2
```

With function L
```
^(\s+)DO\s+(\d*)\s*(\w+)\s*=\s*(\w+\s*\([^(]*\)),\s*([^(,]+?);?$
$1for ($3 = $4; $3 <= $5; $3++) { // $2
```

Step simple (+) == 1
```
^(\s+)DO\s+(\d*)\s*(\w+)\s*=\s*(\w+)\s*,\s*(\w*)\s*,\s*1;?$
$1for ($3 = $4; $3 <= $5; $3++) { // $2
```

Step simple (+) > 1
```
^(\s+)DO\s+(\d*)\s*(\w+)\s*=\s*(\w+)\s*,\s*(\w*)\s*,\s*(\d+);?$
$1for ($3 = $4; $3 <= $5; $3 += $6) { // $2
```

Step simple (-) == -1
```
^(\s+)DO\s+(\d*)\s*(\w+)\s*=\s*(\w+)\s*,\s*(\w*)\s*,\s*-1;?$
$1for ($3 = $4; $3 >= $5; $3--) { // $2
```

Step simple (-) < -1
```
^(\s+)DO\s+(\d*)\s*(\w+)\s*=\s*(\w+)\s*,\s*(\w*)\s*,\s*-(\d+);?$
$1for ($3 = $4; $3 >= $5; $3 -= $6) { // $2
```

Step unknown sign (+-)
```
^(\s+)DO\s+(\d*)\s*(\w+)\s*=\s*(\w+)\s*,\s*(\w*)\s*,\s*([-]?\w+);?$
$1for ($3 = $4; $6 < 0 ? $3 >= $5 : $3 <= $5; $3 += $6) { // $2
```


### CONTINUE (Warning: GOTO statements affected)
```
^(\s+)(\d)(\s+)CONTINUE$
$1 $3} // $2
```

### END DO
```
^(\s+)END\s?DO$
$1}
```

### Simple IF
```
^(\s+)IF\s*\(\s*([^()]*?)\s*\)\s*([^{]*)$
$1if ($2) $3;
```

### Logical operators
```
\.EQ\.
 ==
```

### TRUE/FALSE
Repeat multiple times
```
\s*\.TRUE\.\s*(.+)$
 true $1
```

```
\s*\.TRUE\.$
 true;
```

### Floating point numbers
```
(\d)[DE]0([^\d])
$1$2
```

```
(\d)[DE]0$
$1;
```

```
(\d\.)D0([^\d])
$10$2
```

```
(\d\.)[DE]0$
$10;
```

```
(\.\d+)[DE][+-]0([^\d])
$1$2
```

```
(\.\d+)[DE][+-]0$
$1;
```

```
(\d\.)[DE][+-]0$
$10;
```

```
(\d\.)[DE][+-]0([^\d])
$10$2
```

```
(\d)\.[DE]([+-]?\d+)
$1e$2
```

```
(\d\.\d+)[DE]([+-]?\d+)
$1e$2
```

### Remove IMPLICIT NONE
```
^(\s+)(IMPLICIT NONE);?$
$1// $2
```

### RETURN -> return;
```
^(\s+)RETURN;?$
$1return;
```

### DATA
Repeaters
```
^(\s+DATA\s+.*?/.*)\b6\b\s*\*\s*\b(\d+)\b\s*[, ]?(.*/);?$
$1$2, $2, $2, $2, $2, \2$3;
```

Scalar (non-string)
```
^(\s+)DATA\s+(\w+)\s*/\s*([^',/]+?)\s*/;?$
$1const $2 = $3;
```

Scalar (string)
```
^(\s+)DATA\s+(\w+)\s*/\s*('[^']+?')\s*/;?$
$1const $2 = $3;
```

Single arrays
```
^(\s+)DATA\s+(\w+)\s*/(.*)/;?$
$1const $2 = [$3];
```

Tuples
```
^(\s+)DATA\s+(((\w+)\s*,\s*)+(\w+))\s*/\s*([^/]+)\s*/;?$
$1final ($2) = ($6);
```

Multiple declarations (1st step)
```
^(\s+DATA.+)\s+(\w+)\s*\/\s*(.*?)\s*\/;?$
$1 $2 = $3;
```

Multiple declarations (inner items - n stemp)
```
^(\s+DATA.+\/\s*,)\s+(\w+)\s*\/\s*(.*?)\s*\/\s*,
$1 $2 = $3,
```

Multiple declarations (last step)
```
^(\s+)DATA.+\s+(\w+)\s*\/\s*(.*?)\s*\/\s*,
$1const $2 = $3,
```

Array
```
^(\s+)DATA\s+(\w+)\s*/\s*(.*?)\s*/;?$
$1const $2 = [$3];
```

### IF CALL
```
^(\s+if\s*\(.*\))\s*CALL\s+((\w+)\s*(\(.*\)));?$
$1 \L$3$4;
```

### IF RETURN
```
^((\s+).*)RETURN;?$
$1return;
```

### Overlaped lines
```
^(\s+)if(.*?)      (\s+.*)$
$1if$2;\n     $3
```

### MIN/MAX/SQRT -> min/max/sqrt
```
^(\s+.*)\b(MAX|MIN)\b\s*(\(.*)$
$1\L$2$3
```

### SUBROUTINE
```
^(\s+)SUBROUTINE\s+(\w+)\s*\(\s*(.*)\s*\);?$
$1void \L$2($3) {
```

No params
```
^(\s+)SUBROUTINE\s+(\w+);?$
$1void \L$2() {
```

### FUNCTION
```
^(\s+\w+)\s+FUNCTION\s+(\w+)\s*\(\s*(.*)\s*\);?$
$1 \L$2($3) {
```

### ABS -> abs
```
^(\s+[^/\s].*[^.])\bABS\s*(\([^(]+?\))
$1$2.abs()
```

+ one level
```
^(\s+[^/\s].*[^.])\bABS\s*(\([^(]+?\([^(]+\)\s*\))
$1$2.abs()
```

### Comment out external functions
1st step
```
^(      //.*External Functions.*)\n      ([^/\s].*)$
$1\n      //- $2
```

nth step
```
^(      //.*External Functions.*)\n      (//- .*)\n      ([^/\s].*)$
$1\n      $2\n      //- $3
```
### MOD(x,y) -> (x % y)
Simple
```
\bMOD\(\s*([^(,]+?)\s*,\s*([^(,]+?)\s*\)
($1 % $2)
```

.* on the L
```
\bMOD\(\s*(.*?)\s*,\s*([^(,]+?)\s*\)
($1 % $2)
```

((x),(y)) -> (x % y)
```
\bMOD\(\s*\(\s*([^()]+?)\s*\)\s*,\s*\(\s*([^()]+?)\s*\)\s*\)
($1 % $2)
```

### Array assignment
```
^(\s+\w+)\(\s*([^=]+?)\s*\)\s*=\s+
$1[$2] =
```

### DBLE -> toDouble
A word
```
^(\s+.*)\bDBLE\s*\(\s*(\w+?)\s*\)
$1$2.toDouble()
```

Simple
```
^(\s+.*)\bDBLE\s*\(\s*([^()]+?)\s*\)
$1($2).toDouble()
```

1 level
```
^(\s+.*)\bDBLE\s*\(\s*([^()]*\([^()]+?\)+[^()]*?)\s*\)
$1($2).toDouble()
```

### ++ --
```
(\b\w[^\s]*)\s*=\s*\1\s*([+-])\s1\b
$1$2$2
```

### += -= *= /=
Matrices
```
(\b\w[^\s]*\[[^\[\]]*?\]\[[^\[\]]*?\])\s*=\s*\1\s*([\/*+-])
$1 $2=
```
Vars
```
(\b\w[^\s]*\s*)=\s*\1\s*([*\/+-])
$1$2=
```

### Functios
```
^\s*\b(void|String|int|double|Complex|bool)\s+\w+\((((\s*\n?\s*final\s+\w+(<\w+>)?\s+\w+),\s*)*(\4\n?,?)?)\)\s*\{
^( *\b(void|String|int|double|Complex|bool)\s+\w+\(((\s*final\s+\w+(<\w+>)?\s+\w+)\s*,\s*)*)(\s*final\s+\w+(<\w+>)?\s+\w+\s*,?)(((\s*\n?\s*final\s+\w+(<\w+>)?\s+\w+)\s*,\s*)*\s*\)\s*\{)
^( *\b(void|String|int|double|Complex|bool)\s+\w+\(((\s*final\s+\w+(<\w+>)?\s+\w+)\s*,\s*)*?)((\s*final\s+\w+<\w+>\s+)(\w+)(\s*,?))(((\s*\n?\s*final\s+\w+(<\w+>)?\s+\w+)\s*,\s*)*\s*\)\s*\{)
```
Add trailing `_` to Matrix/Array args
```sh
find . -type f -regex '.*\.dart' | xargs -I % perl -0777 -i'' -pe 's/(\n? *\b(void|String|int|double|Complex|bool)\s+\w+\((([\s\n]*final\s+\w+(<\w+>)?\s+\w+)\s*,[\s\n]*)*)(([\s\n]*final\s+((Matrix|Array)+<\w+>)\s+)(\b(\w*)[^_]\b)([\s\n]*,?))((([\s\n]*final\s+\w+(<\w+>)?\s+\w+)[\s\n]*,?\s*)*\)\s*\{[^\n]*(\n\/\/[^\n]*)*)/$1$7$10_$12$13\\n  final $10 = $10_.dim();/g' %
```

# LAPACK

[![Build Status](https://travis-ci.org/Reference-LAPACK/lapack.svg?branch=master)](https://travis-ci.org/Reference-LAPACK/lapack)
![CMake](https://github.com/Reference-LAPACK/lapack/actions/workflows/cmake.yml/badge.svg)
![Makefile](https://github.com/Reference-LAPACK/lapack/actions/workflows/makefile.yml/badge.svg)
[![Appveyor](https://ci.appveyor.com/api/projects/status/bh38iin398msrbtr?svg=true)](https://ci.appveyor.com/project/langou/lapack/)
[![codecov](https://codecov.io/gh/Reference-LAPACK/lapack/branch/master/graph/badge.svg)](https://codecov.io/gh/Reference-LAPACK/lapack)
[![Packaging status](https://repology.org/badge/tiny-repos/lapack.svg)](https://repology.org/metapackage/lapack/versions)
[![OpenSSF Scorecard](https://api.securityscorecards.dev/projects/github.com/Reference-LAPACK/lapack/badge)](https://securityscorecards.dev/viewer/?uri=github.com/Reference-LAPACK/lapack)

* VERSION 1.0   :  February 29, 1992
* VERSION 1.0a  :  June 30, 1992
* VERSION 1.0b  :  October 31, 1992
* VERSION 1.1   :  March 31, 1993
* VERSION 2.0   :  September 30, 1994
* VERSION 3.0   :  June 30, 1999
* VERSION 3.0 + update :  October 31, 1999
* VERSION 3.0 + update :  May 31, 2000
* VERSION 3.1   : November 2006
* VERSION 3.1.1 : February 2007
* VERSION 3.2   : November 2008
* VERSION 3.2.1 : April 2009
* VERSION 3.2.2 : June 2010
* VERSION 3.3.0 : November 2010
* VERSION 3.3.1 : April 2011
* VERSION 3.4.0 : November 2011
* VERSION 3.4.1 : April 2012
* VERSION 3.4.2 : September 2012
* VERSION 3.5.0 : November 2013
* VERSION 3.6.0 : November 2015
* VERSION 3.6.1 : June 2016
* VERSION 3.7.0 : December 2016
* VERSION 3.7.1 : June 2017
* VERSION 3.8.0 : November 2017
* VERSION 3.9.0 : November 2019
* VERSION 3.9.1 : April 2021
* VERSION 3.10.0 : June 2021
* VERSION 3.10.1 : April 2022
* VERSION 3.11.0 : November 2022
* VERSION 3.12.0 : November 2023

LAPACK is a library of Fortran subroutines for solving the most commonly
occurring problems in numerical linear algebra.

LAPACK is a freely-available software package. It can be included in commercial
software packages (and has been). We only ask that that proper credit be given
to the authors, for example by citing the LAPACK Users' Guide. The license used
for the software is the [modified BSD license](https://github.com/Reference-LAPACK/lapack/blob/master/LICENSE).

Like all software, it is copyrighted. It is not trademarked, but we do ask the
following: if you modify the source for these routines we ask that you change
the name of the routine and comment the changes made to the original.

We will gladly answer any questions regarding the software. If a modification
is done, however, it is the responsibility of the person who modified the
routine to provide support.

LAPACK is [available from GitHub](https://github.com/Reference-LAPACK/lapack).
LAPACK releases are also [available on netlib](http://www.netlib.org/lapack/).

The distribution contains (1) the Fortran source for LAPACK, and (2) its
testing programs.  It also contains (3) the Fortran reference implementation of
the Basic Linear Algebra Subprograms (the Level 1, 2, and 3 BLAS) needed by
LAPACK.  However this code is intended for use only if there is no other
implementation of the BLAS already available on your machine; the efficiency of
LAPACK depends very much on the efficiency of the BLAS.  It also contains (4)
CBLAS, a C interface to the BLAS, and (5) LAPACKE, a C interface to LAPACK.

## Installation

 - LAPACK can be installed with `make`. The configuration must be set in the
   `make.inc` file. A `make.inc.example` for a Linux machine running GNU compilers
   is given in the main directory. Some specific `make.inc` are also available in
   the `INSTALL` directory.
 - LAPACK includes also the [CMake](https://cmake.org/) build.  You will need
   to have CMake installed on your machine.  CMake will allow an easy
   installation on a Windows Machine.  An example CMake build to install the
   LAPACK library under `$HOME/.local/lapack/` is:
   ```sh
   mkdir build
   cd build
   cmake -DCMAKE_INSTALL_LIBDIR=$HOME/.local/lapack ..
   cmake --build . -j --target install
   ```
 - LAPACK can be built and installed using [vcpkg](https://github.com/Microsoft/vcpkg/) dependency manager:
   ```sh
   git clone https://github.com/Microsoft/vcpkg.git
   cd vcpkg
   ./bootstrap-vcpkg.sh  # ./bootstrap-vcpkg.bat for Windows
   ./vcpkg integrate install
   ./vcpkg install lapack
   ```
   The lapack port in vcpkg is kept up to date by Microsoft team members and community contributors. If the version is out of date, please [create an issue or pull request](https://github.com/Microsoft/vcpkg) on the vcpkg repository.

## User Support

LAPACK has been thoroughly tested, on many different types of computers. The
LAPACK project supports the package in the sense that reports of errors or poor
performance will gain immediate attention from the developers. Such reports,
descriptions of interesting applications, and other comments should be sent by
email to [the LAPACK team](mailto:lapack@icl.utk.edu).

A list of known problems, bugs, and compiler errors for LAPACK is
[maintained on netlib](http://www.netlib.org/lapack/release_notes.html).
Please see as well the [GitHub issue tracker](https://github.com/Reference-LAPACK/lapack/issues).

For further information on LAPACK please read our [FAQ](http://www.netlib.org/lapack/faq.html)
and [Users' Guide](http://www.netlib.org/lapack/lug/lapack_lug.html).
A [user forum](http://icl.cs.utk.edu/lapack-forum/) and specific information for
[running LAPACK under Windows](http://icl.cs.utk.edu/lapack-for-windows/lapack/).
is also available to help you with the LAPACK library.


## Testing

LAPACK includes a thorough test suite. We recommend that, after compilation,
you run the test suite.

For complete information on the LAPACK Testing please consult LAPACK Working
Note 41 "Installation Guide for LAPACK".


## LAPACKE

LAPACK now includes the [LAPACKE](http://www.netlib.org/lapack/lapacke.html)
package.  LAPACKE is a Standard C language API for LAPACK that was born from a
collaboration of the LAPACK and INTEL Math Kernel Library teams.
