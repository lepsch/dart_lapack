      double cqpt01(final int M, final int N, final int K, final int A, final int AF, final int LDA, final int TAU, final int JPVT, final Array<double> WORK, final int LWORK,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      int                JPVT( * );
      Complex            A( LDA, * ), AF( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO, J;
      double               NORMA;
      double               RWORK( 1 );
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CUNMQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL

      CQPT01 = ZERO;

      // Test if there is enough workspace

      if ( LWORK < M*N+N ) {
         xerbla('CQPT01', 10 );
         return;
      }

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      NORMA = CLANGE( 'One-norm', M, N, A, LDA, RWORK );

      for (J = 1; J <= K; J++) {
         for (I = 1; I <= min( J, M ); I++) {
            WORK[( J-1 )*M+I] = AF( I, J );
         }
         for (I = J + 1; I <= M; I++) {
            WORK[( J-1 )*M+I] = ZERO;
         }
      }
      for (J = K + 1; J <= N; J++) {
         ccopy(M, AF( 1, J ), 1, WORK( ( J-1 )*M+1 ), 1 );
      }

      cunmqr('Left', 'No transpose', M, N, K, AF, LDA, TAU, WORK, M, WORK( M*N+1 ), LWORK-M*N, INFO );

      for (J = 1; J <= N; J++) {

         // Compare i-th column of QR and jpvt(i)-th column of A

         caxpy(M, CMPLX( -ONE ), A( 1, JPVT( J ) ), 1, WORK( ( J-1 )*M+1 ), 1 );
      }

      CQPT01 = CLANGE( 'One-norm', M, N, WORK, M, RWORK ) / ( REAL( max( M, N ) )*SLAMCH( 'Epsilon' ) )       IF( NORMA != ZERO ) CQPT01 = CQPT01 / NORMA;

      }
