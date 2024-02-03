      double           FUNCTION ZRZT02( M, N, AF, LDA, TAU, WORK, LWORK );
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
*     ..
*     .. Array Arguments ..
      COMPLEX*16         AF( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      int                I, INFO;
*     ..
*     .. Local Arrays ..
      double             RWORK( 1 );
*     ..
*     .. External Functions ..
      double             DLAMCH, ZLANGE;
      EXTERNAL           DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLASET, ZUNMRZ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX
*     ..
*     .. Executable Statements ..
*
      ZRZT02 = ZERO
*
      IF( LWORK.LT.N*N+N ) THEN
         CALL XERBLA( 'ZRZT02', 7 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) RETURN
*
*     Q := I
*
      CALL ZLASET( 'Full', N, N, DCMPLX( ZERO ), DCMPLX( ONE ), WORK, N )
*
*     Q := P(1) * ... * P(m) * Q
*
      CALL ZUNMRZ( 'Left', 'No transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO )
*
*     Q := P(m)' * ... * P(1)' * Q
*
      CALL ZUNMRZ( 'Left', 'Conjugate transpose', N, N, M, N-M, AF, LDA, TAU, WORK, N, WORK( N*N+1 ), LWORK-N*N, INFO )
*
*     Q := Q - I
*
      DO 10 I = 1, N
         WORK( ( I-1 )*N+I ) = WORK( ( I-1 )*N+I ) - ONE
   10 CONTINUE
*
      ZRZT02 = ZLANGE( 'One-norm', N, N, WORK, N, RWORK ) / ( DLAMCH( 'Epsilon' )*DBLE( MAX( M, N ) ) )
      RETURN
*
*     End of ZRZT02
*
      END
