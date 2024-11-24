// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      bool slaisnan(final int SIN1, final int SIN2,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double, INTENT(IN) :: SIN1, SIN2;
      // ..

// =====================================================================

      SLAISNAN = (SIN1 != SIN2);
      }
