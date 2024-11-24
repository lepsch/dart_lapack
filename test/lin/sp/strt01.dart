// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void strt01(final int UPLO, final int DIAG, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> AINV_, final int LDAINV, final int RCOND, final Array<double> _WORK_, final int RESID,) {
  final A = A_.dim();
  final AINV = AINV_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, UPLO;
      int                LDA, LDAINV, N;
      double               RCOND, RESID;
      double               A( LDA, * ), AINV( LDAINV, * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double               AINVNM, ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, SLANTR;
      // EXTERNAL lsame, SLAMCH, SLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL STRMV
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
      ANORM = SLANTR( '1', UPLO, DIAG, N, N, A, LDA, WORK );
      AINVNM = SLANTR( '1', UPLO, DIAG, N, N, AINV, LDAINV, WORK );
      if ( ANORM <= ZERO || AINVNM <= ZERO ) {
         RCOND = ZERO;
         RESID = ONE / EPS;
         return;
      }
      RCOND = ( ONE / ANORM ) / AINVNM;

      // Set the diagonal of AINV to 1 if AINV has unit diagonal.

      if ( lsame( DIAG, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 10
            AINV[J][J] = ONE;
         } // 10
      }

      // Compute A * AINV, overwriting AINV.

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 20
            strmv('Upper', 'No transpose', DIAG, J, A, LDA, AINV( 1, J ), 1 );
         } // 20
      } else {
         for (J = 1; J <= N; J++) { // 30
            strmv('Lower', 'No transpose', DIAG, N-J+1, A( J, J ), LDA, AINV( J, J ), 1 );
         } // 30
      }

      // Subtract 1 from each diagonal element to form A*AINV - I.

      for (J = 1; J <= N; J++) { // 40
         AINV[J][J] = AINV( J, J ) - ONE;
      } // 40

      // Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)

      RESID = SLANTR( '1', UPLO, 'Non-unit', N, N, AINV, LDAINV, WORK );

      RESID = ( ( RESID*RCOND ) / REAL( N ) ) / EPS;

      }
