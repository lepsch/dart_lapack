// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:io';

import 'package:lapack/src/install/dlamch.dart';
import 'package:lapack/src/nio.dart';

void main() {
  final EPS = dlamch('Epsilon');
  final SFMIN = dlamch('Safe minimum');
  final BASE = dlamch('Base');
  final PREC = dlamch('Precision');
  final T = dlamch('Number of digits in mantissa');
  final RND = dlamch('Rounding mode');
  final EMIN = dlamch('Minimum exponent');
  final RMIN = dlamch('Underflow threshold');
  final EMAX = dlamch('Largest exponent');
  final RMAX = dlamch('Overflow threshold');

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
