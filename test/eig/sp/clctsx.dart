// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

bool clctsx(ALPHA, BETA) {
// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

  // .. Scalar Arguments ..
  Complex ALPHA, BETA;
  // ..

// =====================================================================

  // .. Parameters ..
  // REAL               ZERO
  // PARAMETER          ( ZERO = 0.0 )
  // COMPLEX            CZERO
  // PARAMETER          ( CZERO = ( 0.0, 0.0 ) )
  // ..
  // .. Scalars in Common ..
  bool FS;
  int I, M, MPLUSN, N;
  // ..
  // .. Common blocks ..
  // COMMON / MN / M, N, MPLUSN, I, FS
  // ..
  // .. Save statement ..
  SAVE;
  // ..

  if (FS) {
    I = I + 1;
    if (I <= M) {
      CLCTSX = false;
    } else {
      CLCTSX = true;
    }
    if (I == MPLUSN) {
      FS = false;
      I = 0;
    }
  } else {
    I = I + 1;
    if (I <= N) {
      CLCTSX = true;
    } else {
      CLCTSX = false;
    }
    if (I == MPLUSN) {
      FS = true;
      I = 0;
    }
  }

  // IF( BETA == CZERO ) THEN
  // CLCTSX = ( REAL( ALPHA ) > ZERO )
  // ELSE
  // CLCTSX = ( REAL( ALPHA/BETA ) > ZERO )
  // END IF
}
