# GGCount

### A C++ program for fast computation of the number of paths in a grid graph.

This program computes the number of corner-to-corner simple paths and the
number of simple cycles in a grid graph.
It have extended integer sequences [A007764](http://oeis.org/A007764) and
[A140517](http://oeis.org/A140517) in November 2013.
Technical details are shown in the following report: 

> H. Iwashita, Y. Nakazawa, J. Kawahara, T. Uno, and S. Minato,
  "Efficient Computation of the Number of Paths in a Grid Graph with Minimal
  Perfect Hash Functions",
  Technical Report TCS-TR-A-13-64, Division of Computer Science, Graduate
  School of Information Science and Technology, Hokkaido University, 2013.
  ([pdf](http://www-alg.ist.hokudai.ac.jp/~thomas/TCSTR/tcstr_13_64/tcstr_13_64.pdf))

## Requirements

Modern Linux or Linux-like environment including:

* GCC with C++11 and OpenMP features
* GNU Make

This program is tested only on 64-bit Linux.

## Compilation

Just type `make` in the source directory to get both single-thread version
(`ggcount_single`) and multi-thread version (`ggcount_multi`) of the
executable program.  Please check the `Makefile` if necessary.

## Usage

> `$` `./ggcount_single` [ `-c` ] [ `-h` ] _columns_ [ _rows_ ] [ `%`_modulus_ ]

> `$` `./ggcount_multi` [ `-c` ] [ `-h` ] _columns_ [ _rows_ ] [ `%`_modulus_ ]

> `$` `./ggcount.pl` [ `-c` ] [ `-h` ] _size_

Arguments _columns_ and _rows_ specify the horizontal and vertical sizes of a
grid graph.  The latter one can be ommited if they are the same.

When option `-c` is given, it computes the number of cycles in the graph;
otherwise it computes the number of paths connecting opposite corners of the
graph.  Option `-h` is a switch to count Hamiltonian paths/cycles.

Optional argument `%`_modulus_ makes the program use modular arithmetic rather
than bignum arithmetic.  The _modulus_ value must not be greater than 64-bit
unsigned integer.

`ggcount.pl` is the Perl utility that repeats executution of `./ggcount_multi`
with coprime `%`_modulus_ values and computes the final answer from those
results using Chinese remainder theorem.
It only supports square grid graphs.

## See also

* [Time with class! Let's count!](http://www.youtube.com/watch?v=Q4gTV4r0zRs):
  the Youtube-animation that have motivated us to tackle this problem.

* [Graphillion](http://graphillion.org): a Python library for those who want
  to compute something more than just counting.

