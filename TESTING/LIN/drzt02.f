      double           FUNCTION DRZT02( M, N, AF, LDA, TAU, WORK, LWORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             AF( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO;
      // ..
      // .. Local Arrays ..
      double             RWORK( 1 );
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE;
      // EXTERNAL DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET, DORMRZ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX
      // ..
      // .. Executable Statements ..

      DRZT02 = ZERO

      if ( LWORK.LT.N*N+N ) {
         xerbla('DRZT02', 7 );
         RETURN
      }

      // Quick return if possible

      IF( M.LE.0 .OR. N.LE.0 ) RETURN

      // Q := I

      dlaset('Full', N, N, ZERO, ONE, WORK, N );

      // Q := P(1) * ... * P(m) * Q

      dormrz('Left', 'No transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO );

      // Q := P(m) * ... * P(1) * Q

      dormrz('Left', 'Transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO );

      // Q := Q - I

      DO 10 I = 1, N
         WORK( ( I-1 )*N+I ) = WORK( ( I-1 )*N+I ) - ONE
   10 CONTINUE

      DRZT02 = DLANGE( 'One-norm', N, N, WORK, N, RWORK ) / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N ) ) )
      RETURN

      // End of DRZT02

      }
