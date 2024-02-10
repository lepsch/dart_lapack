      double drzt02(M, N, AF, LDA, TAU, final Array<double> WORK, final int LWORK) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LWORK, M, N;
      double             AF( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO;
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET, DORMRZ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX

      DRZT02 = ZERO;

      if ( LWORK < N*N+N ) {
         xerbla('DRZT02', 7 );
         return;
      }

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      // Q := I

      dlaset('Full', N, N, ZERO, ONE, WORK, N );

      // Q := P(1) * ... * P(m) * Q

      dormrz('Left', 'No transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO );

      // Q := P(m) * ... * P(1) * Q

      dormrz('Left', 'Transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO );

      // Q := Q - I

      for (I = 1; I <= N; I++) { // 10
         WORK[( I-1 )*N+I] = WORK( ( I-1 )*N+I ) - ONE;
      } // 10

      DRZT02 = dlange( 'One-norm', N, N, WORK, N, RWORK ) / ( dlamch( 'Epsilon' )*(max( M, N )).toDouble() );
      }
