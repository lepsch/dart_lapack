      double dqpt01(M, N, K, A, AF, LDA, TAU, JPVT, WORK, LWORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      int                JPVT( * );
      double             A( LDA, * ), AF( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO, J;
      double             NORMA;
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DORMQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN

      DQPT01 = ZERO;

      // Test if there is enough workspace

      if ( LWORK < M*N+N ) {
         xerbla('DQPT01', 10 );
         return;
      }

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      NORMA = dlange( 'One-norm', M, N, A, LDA, RWORK );

      for (J = 1; J <= K; J++) {

         // Copy the upper triangular part of the factor R stored
         // in AF(1:K,1:K) into the work array WORK.

         for (I = 1; I <= min( J, M ); I++) {
            WORK[( J-1 )*M+I] = AF( I, J );
         }

         // Zero out the elements below the diagonal in the work array.

         for (I = J + 1; I <= M; I++) {
            WORK[( J-1 )*M+I] = ZERO;
         }
      }

      // Copy columns (K+1,N) from AF into the work array WORK.
      // AF(1:K,K+1:N) contains the rectangular block of the upper trapezoidal
      // factor R, AF(K+1:M,K+1:N) contains the partially updated residual
      // matrix of R.

      for (J = K + 1; J <= N; J++) {
         dcopy(M, AF( 1, J ), 1, WORK( ( J-1 )*M+1 ), 1 );
      }

      dormqr('Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO );

      for (J = 1; J <= N; J++) {

         // Compare J-th column of QR and JPVT(J)-th column of A.

         daxpy(M, -ONE, A( 1, JPVT( J ) ), 1, WORK( ( J-1 )*M+1 ), 1 );
      }

      DQPT01 = dlange( 'One-norm', M, N, WORK, M, RWORK ) / ( (max( M, N )).toDouble()*dlamch( 'Epsilon' ) )       IF( NORMA != ZERO ) DQPT01 = DQPT01 / NORMA;

      }
