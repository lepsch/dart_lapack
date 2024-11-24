// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void csyt03(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> AINV_, final int LDAINV, final Matrix<double> WORK_, final int LDWORK, final Array<double> RWORK_, final int RCOND, final int RESID,) {
  final A = A_.dim();
  final AINV = AINV_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, LDAINV, LDWORK, N;
      double               RCOND, RESID;
      double               RWORK( * );
      Complex            A( LDA, * ), AINV( LDAINV, * ), WORK( LDWORK, * );
      // ..

// =====================================================================


      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I, J;
      double               AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANGE, CLANSY, SLAMCH;
      // EXTERNAL lsame, CLANGE, CLANSY, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSYMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL

      // Quick exit if N = 0

      if ( N <= 0 ) {
         RCOND = ONE;
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANSY( '1', UPLO, N, A, LDA, RWORK );
      AINVNM = CLANSY( '1', UPLO, N, AINV, LDAINV, RWORK );
      if ( ANORM <= ZERO || AINVNM <= ZERO ) {
         RCOND = ZERO;
         RESID = ONE / EPS;
         return;
      }
      RCOND = ( ONE/ANORM ) / AINVNM;

      // Expand AINV into a full matrix and call CSYMM to multiply
      // AINV on the left by A (store the result in WORK).

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J - 1; I++) { // 10
               AINV[J][I] = AINV( I, J );
            } // 10
         } // 20
      } else {
         for (J = 1; J <= N; J++) { // 40
            for (I = J + 1; I <= N; I++) { // 30
               AINV[J][I] = AINV( I, J );
            } // 30
         } // 40
      }
      csymm('Left', UPLO, N, N, -CONE, A, LDA, AINV, LDAINV, CZERO, WORK, LDWORK );

      // Add the identity matrix to WORK .

      for (I = 1; I <= N; I++) { // 50
         WORK[I][I] = WORK( I, I ) + CONE;
      } // 50

      // Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)

      RESID = CLANGE( '1', N, N, WORK, LDWORK, RWORK );

      RESID = ( ( RESID*RCOND )/EPS ) / REAL( N );

      }
