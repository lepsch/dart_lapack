// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void strt02(final int UPLO, final int TRANS, final int DIAG, final int N, final int NRHS, final Matrix<double> A_, final int LDA, final Matrix<double> X_, final int LDX, final Matrix<double> B_, final int LDB, final Array<double> _WORK_, final int RESID,) {
  final A = A_.dim();
  final X = X_.dim();
  final B = B_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      double               RESID;
      double               A( LDA, * ), B( LDB, * ), WORK( * ), X( LDX, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                J;
      double               ANORM, BNORM, EPS, XNORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SASUM, SLAMCH, SLANTR;
      // EXTERNAL lsame, SASUM, SLAMCH, SLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, STRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Quick exit if N = 0 or NRHS = 0

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Compute the 1-norm of op(A).

      if ( lsame( TRANS, 'N' ) ) {
         ANORM = SLANTR( '1', UPLO, DIAG, N, N, A, LDA, WORK );
      } else {
         ANORM = SLANTR( 'I', UPLO, DIAG, N, N, A, LDA, WORK );
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the maximum over the number of right hand sides of
      //    norm(op(A)*X - B) / ( norm(op(A)) * norm(X) * EPS )

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 10
         scopy(N, X( 1, J ), 1, WORK, 1 );
         strmv(UPLO, TRANS, DIAG, N, A, LDA, WORK, 1 );
         saxpy(N, -ONE, B( 1, J ), 1, WORK, 1 );
         BNORM = SASUM( N, WORK, 1 );
         XNORM = SASUM( N, X( 1, J ), 1 );
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS );
         }
      } // 10

      }
