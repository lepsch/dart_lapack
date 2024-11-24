// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      void sget54(final int N, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final Matrix<double> S_, final int LDS, final Matrix<double> T_, final int LDT, final Matrix<double> U_, final int LDU, final Matrix<double> V_, final int LDV, final Array<double> _WORK_, final int RESULT,) {
  final A = A_.dim();
  final B = B_.dim();
  final S = S_.dim();
  final T = T_.dim();
  final U = U_.dim();
  final V = V_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LDS, LDT, LDU, LDV, N;
      double               RESULT;
      double               A( LDA, * ), B( LDB, * ), S( LDS, * ), T( LDT, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               ABNORM, ULP, UNFL, WNORM;
      double               DUM( 1 );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE;
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL

      RESULT = ZERO;
      if (N <= 0) return;

      // Constants

      UNFL = SLAMCH( 'Safe minimum' );
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );

      // compute the norm of (A,B)

      slacpy('Full', N, N, A, LDA, WORK, N );
      slacpy('Full', N, N, B, LDB, WORK( N*N+1 ), N );
      ABNORM = max( SLANGE( '1', N, 2*N, WORK, N, DUM ), UNFL );

      // Compute W1 = A - U*S*V', and put in the array WORK(1:N*N)

      slacpy(' ', N, N, A, LDA, WORK, N );
      sgemm('N', 'N', N, N, N, ONE, U, LDU, S, LDS, ZERO, WORK( N*N+1 ), N );

      sgemm('N', 'C', N, N, N, -ONE, WORK( N*N+1 ), N, V, LDV, ONE, WORK, N );

      // Compute W2 = B - U*T*V', and put in the workarray W(N*N+1:2*N*N)

      slacpy(' ', N, N, B, LDB, WORK( N*N+1 ), N );
      sgemm('N', 'N', N, N, N, ONE, U, LDU, T, LDT, ZERO, WORK( 2*N*N+1 ), N );

      sgemm('N', 'C', N, N, N, -ONE, WORK( 2*N*N+1 ), N, V, LDV, ONE, WORK( N*N+1 ), N );

      // Compute norm(W)/ ( ulp*norm((A,B)) )

      WNORM = SLANGE( '1', N, 2*N, WORK, N, DUM );

      if ( ABNORM > WNORM ) {
         RESULT = ( WNORM / ABNORM ) / ( 2*N*ULP );
      } else {
         if ( ABNORM < ONE ) {
            RESULT = ( min( WNORM, 2*N*ABNORM ) / ABNORM ) / ( 2*N*ULP );
         } else {
            RESULT = min( WNORM / ABNORM, REAL( 2*N ) ) / ( 2*N*ULP );
         }
      }

      }
