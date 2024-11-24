// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void cptsvx(final int FACT, final int N, final int NRHS, final int D, final int E, final int DF, final int EF, final Matrix<double> B_, final int LDB, final Matrix<double> X_, final int LDX, final int RCOND, final int FERR, final int BERR, final Array<double> _WORK_, final Array<double> RWORK_, final Box<int> INFO,) {
  final B = B_.dim();
  final X = X_.dim();
  final _WORK = _WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             FACT;
      int                INFO, LDB, LDX, N, NRHS;
      double               RCOND;
      double               BERR( * ), D( * ), DF( * ), FERR( * ), RWORK( * )       Complex            B( LDB, * ), E( * ), EF( * ), WORK( * ), X( LDX, * );
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      bool               NOFACT;
      double               ANORM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANHT, SLAMCH;
      // EXTERNAL lsame, CLANHT, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CLACPY, CPTCON, CPTRFS, CPTTRF, CPTTRS, SCOPY, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      NOFACT = lsame( FACT, 'N' );
      if ( !NOFACT && !lsame( FACT, 'F' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( LDX < max( 1, N ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('CPTSVX', -INFO );
         return;
      }

      if ( NOFACT ) {

         // Compute the L*D*L**H (or U**H*D*U) factorization of A.

         scopy(N, D, 1, DF, 1 );
         if (N > 1) ccopy( N-1, E, 1, EF, 1 );
         cpttrf(N, DF, EF, INFO );

         // Return if INFO is non-zero.

         if ( INFO > 0 ) {
            RCOND = ZERO;
            return;
         }
      }

      // Compute the norm of the matrix A.

      ANORM = CLANHT( '1', N, D, E );

      // Compute the reciprocal of the condition number of A.

      cptcon(N, DF, EF, ANORM, RCOND, RWORK, INFO );

      // Compute the solution vectors X.

      clacpy('Full', N, NRHS, B, LDB, X, LDX );
      cpttrs('Lower', N, NRHS, DF, EF, X, LDX, INFO );

      // Use iterative refinement to improve the computed solutions and
      // compute error bounds and backward error estimates for them.

      cptrfs('Lower', N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

      // Set INFO = N+1 if the matrix is singular to working precision.

      if( RCOND < SLAMCH( 'Epsilon' ) ) INFO = N + 1;

      }
