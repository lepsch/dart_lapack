// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cbdt02(final int M, final int N, final Matrix<double> B_, final int LDB, final Matrix<double> C_, final int LDC, final Matrix<double> U_, final int LDU, final Array<double> _WORK_, final Array<double> RWORK_, final int RESID,) {
  final B = B_.dim();
  final C = C_.dim();
  final U = U_.dim();
  final _WORK = _WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDB, LDC, LDU, M, N;
      double               RESID;
      double               RWORK( * );
      Complex            B( LDB, * ), C( LDC, * ), U( LDU, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double               BNORM, EPS, REALMN;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SCASUM, SLAMCH;
      // EXTERNAL CLANGE, SCASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL

      // Quick return if possible

      RESID = ZERO;
      if (M <= 0 || N <= 0) return;
      REALMN = double( max( M, N ) );
      EPS = SLAMCH( 'Precision' );

      // Compute norm(B - U * C)

      for (J = 1; J <= N; J++) { // 10
         ccopy(M, B( 1, J ), 1, WORK, 1 );
         cgemv('No transpose', M, M, -CMPLX( ONE ), U, LDU, C( 1, J ), 1, CMPLX( ONE ), WORK, 1 );
         RESID = max( RESID, SCASUM( M, WORK, 1 ) );
      } // 10

      // Compute norm of B.

      BNORM = CLANGE( '1', M, N, B, LDB, RWORK );

      if ( BNORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( BNORM >= RESID ) {
            RESID = ( RESID / BNORM ) / ( REALMN*EPS );
         } else {
            if ( BNORM < ONE ) {
               RESID = ( min( RESID, REALMN*BNORM ) / BNORM ) / ( REALMN*EPS );
            } else {
               RESID = min( RESID / BNORM, REALMN ) / ( REALMN*EPS );
            }
         }
      }
      }
