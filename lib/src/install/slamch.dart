// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/install/lsame.dart';
import 'package:dart_lapack/src/intrinsics/digits.dart';
import 'package:dart_lapack/src/intrinsics/epsilon.dart';
import 'package:dart_lapack/src/intrinsics/huge.dart';
import 'package:dart_lapack/src/intrinsics/maxexponent.dart';
import 'package:dart_lapack/src/intrinsics/minexponent.dart';
import 'package:dart_lapack/src/intrinsics/radix.dart';
import 'package:dart_lapack/src/intrinsics/tiny.dart';

double slamch(final String CMACH) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  const ONE = 1.0, ZERO = 0.0;
  double RND, EPS, SFMIN, SMALL;

  // Assume rounding, not chopping. Always.

  RND = ONE;

  if (ONE == RND) {
    EPS = epsilon32(ZERO) * 0.5;
  } else {
    EPS = epsilon32(ZERO);
  }

  if (lsame(CMACH, 'E')) {
    return EPS;
  } else if (lsame(CMACH, 'S')) {
    SFMIN = tiny32(ZERO);
    SMALL = ONE / huge32(ZERO);
    if (SMALL >= SFMIN) {
      // Use SMALL plus a bit, to avoid the possibility of rounding
      // causing overflow when computing  1/sfmin.
      SFMIN = SMALL * (ONE + EPS);
    }
    return SFMIN;
  } else if (lsame(CMACH, 'B')) {
    return radix32(ZERO);
  } else if (lsame(CMACH, 'P')) {
    return EPS * radix32(ZERO);
  } else if (lsame(CMACH, 'N')) {
    return digits32(ZERO).toDouble();
  } else if (lsame(CMACH, 'R')) {
    return RND;
  } else if (lsame(CMACH, 'M')) {
    return minexponent32(ZERO).toDouble();
  } else if (lsame(CMACH, 'U')) {
    return tiny32(ZERO);
  } else if (lsame(CMACH, 'L')) {
    return maxexponent32(ZERO).toDouble();
  } else if (lsame(CMACH, 'O')) {
    return huge32(ZERO);
  } else {
    return ZERO;
  }
}

// ***********************************************************************
// > \brief \b SLAMC3
// > \details
// > \b Purpose:
// > \verbatim
// > SLAMC3  is intended to force  A  and  B  to be stored prior to doing
// > the addition of  A  and  B ,  for use in situations where optimizers
// > might hold one of these in a register.
// > \endverbatim
// > \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
// >
// > \ingroup lamc3
// >
// > \param[in] A
// > \verbatim
// > \endverbatim
// >
// > \param[in] B
// > \verbatim
// >          The values A and B.
// > \endverbatim
// >
double slamc3(final double A, final double B) {
// -- LAPACK auxiliary routine --
// Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  return A + B;
}
