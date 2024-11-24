// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      double second() {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

// =====================================================================

      // .. Local Scalars ..
      double               T1;
      double               TARRAY( 2 );
      // ..
      // .. Intrinsic Functions ..
      double               ETIME;
      // INTRINSIC ETIME

      T1 = ETIME( TARRAY );
      SECOND = TARRAY( 1 );
      }
