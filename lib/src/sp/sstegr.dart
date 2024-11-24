// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void sstegr(final int JOBZ, final int RANGE, final int N, final int D, final int E, final int VL, final int VU, final int IL, final int IU, final int ABSTOL, final int M, final int W, final Matrix<double> Z_, final int LDZ, final int ISUPPZ, final Array<double> WORK_, final int LWORK, final Array<int> IWORK_, final int LIWORK, final Box<int> INFO,) {
  final Z = Z_.dim();
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, RANGE;
      int                IL, INFO, IU, LDZ, LIWORK, LWORK, M, N;
      double             ABSTOL, VL, VU;
      int                ISUPPZ( * ), IWORK( * );
      double               D( * ), E( * ), W( * ), WORK( * );
      double               Z( LDZ, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool    TRYRAC;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSTEMR
      INFO = 0;
      TRYRAC = false;
       sstemr(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK, LWORK, IWORK, LIWORK, INFO );
      }
