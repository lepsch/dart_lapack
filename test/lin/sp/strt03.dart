// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void strt03(final int UPLO, final int TRANS, final int DIAG, final int N, final int NRHS, final Matrix<double> A_, final int LDA, final int SCALE, final int CNORM, final int TSCAL, final Matrix<double> X_, final int LDX, final Matrix<double> B_, final int LDB, final Array<double> _WORK_, final int RESID,) {
  final A = A_.dim();
  final X = X_.dim();
  final B = B_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIAG, TRANS, UPLO;
      int                LDA, LDB, LDX, N, NRHS;
      double               RESID, SCALE, TSCAL;
      double               A( LDA, * ), B( LDB, * ), CNORM( * ), WORK( * ), X( LDX, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                IX, J;
      double               BIGNUM, EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ISAMAX;
      //- REAL               SLAMCH;
      // EXTERNAL lsame, ISAMAX, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SSCAL, STRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, REAL

      // Quick exit if N = 0

      if ( N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }
      EPS = SLAMCH( 'Epsilon' );
      SMLNUM = SLAMCH( 'Safe minimum' );
      BIGNUM = ONE / SMLNUM;

      // Compute the norm of the triangular matrix A using the column
      // norms already computed by SLATRS.

      TNORM = ZERO;
      if ( lsame( DIAG, 'N' ) ) {
         for (J = 1; J <= N; J++) { // 10
            TNORM = max( TNORM, TSCAL*( A( J, J ) ).abs()+CNORM( J ) );
         } // 10
      } else {
         for (J = 1; J <= N; J++) { // 20
            TNORM = max( TNORM, TSCAL+CNORM( J ) );
         } // 20
      }

      // Compute the maximum over the number of right hand sides of
      //    norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 30
         scopy(N, X( 1, J ), 1, WORK, 1 );
         IX = ISAMAX( N, WORK, 1 );
         XNORM = max( ONE, ( X( IX, J ) ).abs() );
         XSCAL = ( ONE / XNORM ) / REAL( N );
         sscal(N, XSCAL, WORK, 1 );
         strmv(UPLO, TRANS, DIAG, N, A, LDA, WORK, 1 );
         saxpy(N, -SCALE*XSCAL, B( 1, J ), 1, WORK, 1 );
         IX = ISAMAX( N, WORK, 1 );
         ERR = TSCAL*( WORK( IX ) ).abs();
         IX = ISAMAX( N, X( 1, J ), 1 );
         XNORM = ( X( IX, J ) ).abs();
         if ( ERR*SMLNUM <= XNORM ) {
            if (XNORM > ZERO) ERR = ERR / XNORM;
         } else {
            if (ERR > ZERO) ERR = ONE / EPS;
         }
         if ( ERR*SMLNUM <= TNORM ) {
            if (TNORM > ZERO) ERR = ERR / TNORM;
         } else {
            if (ERR > ZERO) ERR = ONE / EPS;
         }
         RESID = max( RESID, ERR );
      } // 30

      }
