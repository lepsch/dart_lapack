// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:io';

import 'package:dart_lapack/src/install/slamch.dart';
import 'package:dart_lapack/src/nio.dart';

void main() {
  final EPS = slamch('Epsilon');
  final SFMIN = slamch('Safe minimum');
  final BASE = slamch('Base');
  final PREC = slamch('Precision');
  final T = slamch('Number of digits in mantissa');
  final RND = slamch('Rounding mode');
  final EMIN = slamch('Minimum exponent');
  final RMIN = slamch('Underflow threshold');
  final EMAX = slamch('Largest exponent');
  final RMAX = slamch('Overflow threshold');

  final nout = Nout(stdout);

  nout.println(' Epsilon                      = $EPS');
  nout.println(' Safe minimum                 = $SFMIN');
  nout.println(' Base                         = $BASE');
  nout.println(' Precision                    = $PREC');
  nout.println(' Number of digits in mantissa = $T');
  nout.println(' Rounding mode                = $RND');
  nout.println(' Minimum exponent             = $EMIN');
  nout.println(' Underflow threshold          = $RMIN');
  nout.println(' Largest exponent             = $EMAX');
  nout.println(' Overflow threshold           = $RMAX');
  nout.println(' Reciprocal of safe minimum   = ${1 / SFMIN}');
}
