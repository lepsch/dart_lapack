// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void sla_lin_berr(final int N, final int NZ, final int NRHS, final int RES, final int AYB, final int BERR,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                N, NZ, NRHS;
      double               AYB( N, NRHS ), BERR( NRHS );
      double               RES( N, NRHS );
      // ..

// =====================================================================

      // .. Local Scalars ..
      double               TMP;
      int                I, J;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. External Functions ..
      // EXTERNAL SLAMCH
      double               SLAMCH;
      double               SAFE1;

      // Adding SAFE1 to the numerator guards against spuriously zero
      // residuals.  A similar safeguard is in the SLA_yyAMV routine used
      // to compute AYB.

      SAFE1 = SLAMCH( 'Safe minimum' );
      SAFE1 = (NZ+1)*SAFE1;

      for (J = 1; J <= NRHS; J++) {
         BERR[J] = 0.0;
         for (I = 1; I <= N; I++) {
            if (AYB(I,J) != 0.0) {
               TMP = (SAFE1+(RES(I,J) ).abs() )/AYB(I,J);
               BERR[J] = max( BERR(J), TMP );
            }

      // If AYB is exactly 0.0 (and if computed by SLA_yyAMV), then we know
      // the true residual also must be exactly 0.0.

         }
      }
      }
