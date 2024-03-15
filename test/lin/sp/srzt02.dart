      double srzt02(final int M, final int N, final int AF, final int LDA, final int TAU, final Array<double> WORK_, final int LWORK,) {
  final WORK = WORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LWORK, M, N;
      double               AF( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO;
      double               RWORK( 1 );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE;
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASET, SORMRZ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL

      SRZT02 = ZERO;

      if ( LWORK < N*N+N ) {
         xerbla('SRZT02', 7 );
         return;
      }

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      // Q := I

      slaset('Full', N, N, ZERO, ONE, WORK, N );

      // Q := P(1) * ... * P(m) * Q

      sormrz('Left', 'No transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO );

      // Q := P(m) * ... * P(1) * Q

      sormrz('Left', 'Transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO );

      // Q := Q - I

      for (I = 1; I <= N; I++) { // 10
         WORK[( I-1 )*N+I] = WORK( ( I-1 )*N+I ) - ONE;
      } // 10

      SRZT02 = SLANGE( 'One-norm', N, N, WORK, N, RWORK ) / ( SLAMCH( 'Epsilon' )*REAL( max( M, N ) ) );
      }