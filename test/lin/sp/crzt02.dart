      double crzt02(M, N, AF, LDA, TAU, WORK, LWORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LWORK, M, N;
      Complex            AF( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO;
      double               RWORK( 1 );
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SLAMCH;
      // EXTERNAL CLANGE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASET, CUNMRZ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL

      CRZT02 = ZERO;

      if ( LWORK < N*N+N ) {
         xerbla('CRZT02', 7 );
         return;
      }

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      // Q := I

      claset('Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), WORK, N );

      // Q := P(1) * ... * P(m) * Q

      cunmrz('Left', 'No transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO );

      // Q := P(m)' * ... * P(1)' * Q

      cunmrz('Left', 'Conjugate transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO );

      // Q := Q - I

      for (I = 1; I <= N; I++) { // 10
         WORK[( I-1 )*N+I] = WORK( ( I-1 )*N+I ) - ONE;
      } // 10

      CRZT02 = CLANGE( 'One-norm', N, N, WORK, N, RWORK ) / ( SLAMCH( 'Epsilon' )*REAL( max( M, N ) ) );
      return;
      }
