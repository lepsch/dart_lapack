// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      bool sisnan(final int SIN,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double, INTENT(IN) :: SIN;
      // ..

// =====================================================================

// .. External Functions ..
      bool    SLAISNAN;
      // EXTERNAL SLAISNAN
// ..
// .. Executable Statements ..
      SISNAN = SLAISNAN(SIN,SIN);
      }
