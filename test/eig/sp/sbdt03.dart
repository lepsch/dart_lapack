// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void sbdt03(final int UPLO, final int N, final int KD, final int D, final int E, final Matrix<double> U_, final int LDU, final int S, final Matrix<double> VT_, final int LDVT, final Array<double> _WORK_, final int RESID,) {
  final U = U_.dim();
  final VT = VT_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                KD, LDU, LDVT, N;
      double               RESID;
      double               D( * ), E( * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * );
      // ..

// ======================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, J;
      double               BNORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ISAMAX;
      //- REAL               SASUM, SLAMCH;
      // EXTERNAL lsame, ISAMAX, SASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL

      // Quick return if possible

      RESID = ZERO;
      if (N <= 0) return;

      // Compute B - U * S * V' one column at a time.

      BNORM = ZERO;
      if ( KD >= 1 ) {

         // B is bidiagonal.

         if ( lsame( UPLO, 'U' ) ) {

            // B is upper bidiagonal.

            for (J = 1; J <= N; J++) { // 20
               for (I = 1; I <= N; I++) { // 10
                  WORK[N+I] = S( I )*VT( I, J );
               } // 10
               sgemv('No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 );
               WORK[J] = WORK( J ) + D( J );
               if ( J > 1 ) {
                  WORK[J-1] = WORK( J-1 ) + E( J-1 );
                  BNORM = max( BNORM, ( D( J ) ).abs()+( E( J-1 ) ).abs() );
               } else {
                  BNORM = max( BNORM, ( D( J ) ).abs() );
               }
               RESID = max( RESID, SASUM( N, WORK, 1 ) );
            } // 20
         } else {

            // B is lower bidiagonal.

            for (J = 1; J <= N; J++) { // 40
               for (I = 1; I <= N; I++) { // 30
                  WORK[N+I] = S( I )*VT( I, J );
               } // 30
               sgemv('No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 );
               WORK[J] = WORK( J ) + D( J );
               if ( J < N ) {
                  WORK[J+1] = WORK( J+1 ) + E( J );
                  BNORM = max( BNORM, ( D( J ) ).abs()+( E( J ) ).abs() );
               } else {
                  BNORM = max( BNORM, ( D( J ) ).abs() );
               }
               RESID = max( RESID, SASUM( N, WORK, 1 ) );
            } // 40
         }
      } else {

         // B is diagonal.

         for (J = 1; J <= N; J++) { // 60
            for (I = 1; I <= N; I++) { // 50
               WORK[N+I] = S( I )*VT( I, J );
            } // 50
            sgemv('No transpose', N, N, -ONE, U, LDU, WORK( N+1 ), 1, ZERO, WORK, 1 );
            WORK[J] = WORK( J ) + D( J );
            RESID = max( RESID, SASUM( N, WORK, 1 ) );
         } // 60
         J = ISAMAX( N, D, 1 );
         BNORM = ( D( J ) ).abs();
      }

      // Compute norm(B - U * S * V') / ( n * norm(B) * EPS )

      EPS = SLAMCH( 'Precision' );

      if ( BNORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         if ( BNORM >= RESID ) {
            RESID = ( RESID / BNORM ) / ( REAL( N )*EPS );
         } else {
            if ( BNORM < ONE ) {
               RESID = ( min( RESID, double( N )*BNORM ) / BNORM ) / ( REAL( N )*EPS );
            } else {
               RESID = min( RESID / BNORM, REAL( N ) ) / ( REAL( N )*EPS );
            }
         }
      }

      }
