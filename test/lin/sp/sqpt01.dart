// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

      double sqpt01(final int M, final int N, final int K, final int A, final int AF, final int LDA, final int TAU, final int JPVT, final Array<double> WORK_, final int LWORK,) {
  final WORK = WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      int                JPVT( * );
      double               A( LDA, * ), AF( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO, J;
      double               NORMA;
      double               RWORK( 1 );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE;
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SORMQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL

      SQPT01 = ZERO;

      // Test if there is enough workspace

      if ( LWORK < M*N+N ) {
         xerbla('SQPT01', 10 );
         return;
      }

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      NORMA = SLANGE( 'One-norm', M, N, A, LDA, RWORK );

      for (J = 1; J <= K; J++) {
         for (I = 1; I <= min( J, M ); I++) {
            WORK[( J-1 )*M+I] = AF( I, J );
         }
         for (I = J + 1; I <= M; I++) {
            WORK[( J-1 )*M+I] = ZERO;
         }
      }
      for (J = K + 1; J <= N; J++) {
         scopy(M, AF( 1, J ), 1, WORK( ( J-1 )*M+1 ), 1 );
      }

      sormqr('Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO );

      for (J = 1; J <= N; J++) {

         // Compare i-th column of QR and jpvt(i)-th column of A

         saxpy(M, -ONE, A( 1, JPVT( J ) ), 1, WORK( ( J-1 )*M+1 ), 1 );
      }

      SQPT01 = SLANGE( 'One-norm', M, N, WORK, M, RWORK ) / ( REAL( max( M, N ) )*SLAMCH( 'Epsilon' ) )       IF( NORMA != ZERO ) SQPT01 = SQPT01 / NORMA;

      }
