      REAL             FUNCTION SRZT02( M, N, AF, LDA, TAU, WORK, LWORK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               AF( LDA, * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I, INFO;
      // ..
      // .. Local Arrays ..
      REAL               RWORK( 1 )
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE
      // EXTERNAL SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASET, SORMRZ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL
      // ..
      // .. Executable Statements ..

      SRZT02 = ZERO

      IF( LWORK.LT.N*N+N ) THEN
         CALL XERBLA( 'SRZT02', 7 )
         RETURN
      END IF

      // Quick return if possible

      IF( M.LE.0 .OR. N.LE.0 ) RETURN

      // Q := I

      CALL SLASET( 'Full', N, N, ZERO, ONE, WORK, N )

      // Q := P(1) * ... * P(m) * Q

      CALL SORMRZ( 'Left', 'No transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO )

      // Q := P(m) * ... * P(1) * Q

      CALL SORMRZ( 'Left', 'Transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO )

      // Q := Q - I

      DO 10 I = 1, N
         WORK( ( I-1 )*N+I ) = WORK( ( I-1 )*N+I ) - ONE
   10 CONTINUE

      SRZT02 = SLANGE( 'One-norm', N, N, WORK, N, RWORK ) / ( SLAMCH( 'Epsilon' )*REAL( MAX( M, N ) ) )
      RETURN

      // End of SRZT02

      END
