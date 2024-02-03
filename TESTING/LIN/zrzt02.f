      double zrzt02(M, N, AF, LDA, TAU, WORK, LWORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      Complex         AF( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO;
      // ..
      // .. Local Arrays ..
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE;
      // EXTERNAL DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLASET, ZUNMRZ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX
      // ..
      // .. Executable Statements ..

      ZRZT02 = ZERO;

      if ( LWORK < N*N+N ) {
         xerbla('ZRZT02', 7 );
         return;
      }

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      // Q := I

      zlaset('Full', N, N, DCMPLX( ZERO ), DCMPLX( ONE ), WORK, N );

      // Q := P(1) * ... * P(m) * Q

      zunmrz('Left', 'No transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO );

      // Q := P(m)' * ... * P(1)' * Q

      zunmrz('Left', 'Conjugate transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO );

      // Q := Q - I

      for (I = 1; I <= N; I++) { // 10
         WORK( ( I-1 )*N+I ) = WORK( ( I-1 )*N+I ) - ONE;
      } // 10

      ZRZT02 = ZLANGE( 'One-norm', N, N, WORK, N, RWORK ) / ( DLAMCH( 'Epsilon' )*DBLE( max( M, N ) ) );
      return;
      }
