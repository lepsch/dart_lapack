// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/ilaenv.dart';

void main() {
  int IEEEOK;

  print('We are about to check whether infinity arithmetic');
  print('can be trusted.  If this test hangs, set');
  print('ILAENV = 0 for ISPEC = 11 in LAPACK/SRC/ilaenv.f');

  IEEEOK = ilaenv(11, 'ILAENV', 'N', 1, 2, 3, 4);
  print('');

  if (IEEEOK == 0) {
    print('Infinity arithmetic did not perform per the ieee spec');
  } else {
    print('Infinity arithmetic performed as per the ieee spec.');
    print('However, this is not an exhaustive test and does not');
    print('guarantee that infinity arithmetic meets the ieee spec.');
  }

  print('');
  // ilaenv( 10, ...) checks both infinity and NaN arithmetic
  // infinity has already been checked so checking NaN now
  print('We are about to check whether NaN arithmetic');
  print('can be trusted.  If this test hangs, set');
  print('ILAENV = 0 for ISPEC = 10 in LAPACK/SRC/ilaenv.f');
  IEEEOK = ilaenv(10, 'ILAENV', 'N', 1, 2, 3, 4);

  print('');
  if (IEEEOK == 0) {
    print('NaN arithmetic did not perform per the ieee spec');
  } else {
    print('NaN arithmetic performed as per the ieee spec.');
    print('However, this is not an exhaustive test and does not');
    print('guarantee that NaN arithmetic meets the ieee spec.');
  }
  print('');
}
