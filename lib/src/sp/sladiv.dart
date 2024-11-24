// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

void sladiv(A, B, C, D, P, Q) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  double A, B, C, D, P, Q;
  // ..

// =====================================================================

  // .. Parameters ..
  double BS;
  const BS = 2.0;
  double HALF;
  const HALF = 0.5;
  double TWO;
  const TWO = 2.0;

  // .. Local Scalars ..
  double AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS;
  // ..
  // .. External Functions ..
  //- REAL               SLAMCH;
  // EXTERNAL SLAMCH
  // ..
  // .. External Subroutines ..
  // EXTERNAL SLADIV1
  // ..
  // .. Intrinsic Functions ..
  // INTRINSIC ABS, MAX
  // ..

  AA = A;
  BB = B;
  CC = C;
  DD = D;
  AB = max((A).abs(), (B).abs());
  CD = max((C).abs(), (D).abs());
  S = 1.0;

  OV = SLAMCH('Overflow threshold');
  UN = SLAMCH('Safe minimum');
  EPS = SLAMCH('Epsilon');
  BE = BS / (EPS * EPS);

  if (AB >= HALF * OV) {
    AA = HALF * AA;
    BB = HALF * BB;
    S = TWO * S;
  }
  if (CD >= HALF * OV) {
    CC = HALF * CC;
    DD = HALF * DD;
    S = HALF * S;
  }
  if (AB <= UN * BS / EPS) {
    AA = AA * BE;
    BB = BB * BE;
    S = S / BE;
  }
  if (CD <= UN * BS / EPS) {
    CC = CC * BE;
    DD = DD * BE;
    S = S * BE;
  }
  if ((D).abs() <= (C).abs()) {
    sladiv1(AA, BB, CC, DD, P, Q);
  } else {
    sladiv1(BB, AA, DD, CC, P, Q);
    Q = -Q;
  }
  P = P * S;
  Q = Q * S;

  return;
}

// > \ingroup ladiv

void sladiv1(A, B, C, D, P, Q) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  double A, B, C, D, P, Q;
  // ..

// =====================================================================

  // .. Parameters ..
  double ONE;
  const ONE = 1.0;

  // .. Local Scalars ..
  double R, T;
  // ..
  // .. External Functions ..
  //- REAL               SLADIV2;
  // EXTERNAL SLADIV2
  // ..

  R = D / C;
  T = ONE / (C + D * R);
  P = SLADIV2(A, B, C, D, R, T);
  A = -A;
  Q = SLADIV2(B, A, C, D, R, T);

  return;
}

// > \ingroup ladiv

double sladiv2(A, B, C, D, R, T) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  double A, B, C, D, R, T;
  // ..

// =====================================================================

  // .. Parameters ..
  double ZERO;
  const ZERO = 0.0;

  // .. Local Scalars ..
  double BR;
  // ..

  if (R != ZERO) {
    BR = B * R;
    if (BR != ZERO) {
      SLADIV2 = (A + BR) * T;
    } else {
      SLADIV2 = A * T + (B * T) * R;
    }
  } else {
    SLADIV2 = (A + D * (B / C)) * T;
  }
}
