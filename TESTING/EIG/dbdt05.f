      SUBROUTINE DBDT05( M, N, A, LDA, S, NS, U, LDU, VT, LDVT, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDU, LDVT, M, N, NS
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
*     ..
*
* ======================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ANORM, EPS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM, DLAMCH, DLANGE
      EXTERNAL           LSAME, IDAMAX, DASUM, DLAMCH, DLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      RESID = ZERO
      IF( MIN( M, N ).LE.0 .OR. NS.LE.0 ) RETURN
*
      EPS = DLAMCH( 'Precision' )
      ANORM = DLANGE( 'M', M, N, A, LDA, WORK )
*
*     Compute U' * A * V.
*
      CALL DGEMM( 'N', 'T', M, NS, N, ONE, A, LDA, VT, LDVT, ZERO, WORK( 1+NS*NS ), M )       CALL DGEMM( 'T', 'N', NS, NS, M, -ONE, U, LDU, WORK( 1+NS*NS ), M, ZERO, WORK, NS )
*
*     norm(S - U' * B * V)
*
      J = 0
      DO 10 I = 1, NS
         WORK( J+I ) =  WORK( J+I ) + S( I )
         RESID = MAX( RESID, DASUM( NS, WORK( J+1 ), 1 ) )
         J = J + NS
   10 CONTINUE
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         IF( ANORM.GE.RESID ) THEN
            RESID = ( RESID / ANORM ) / ( DBLE( N )*EPS )
         ELSE
            IF( ANORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, DBLE( N )*ANORM ) / ANORM ) / ( DBLE( N )*EPS )
            ELSE
               RESID = MIN( RESID / ANORM, DBLE( N ) ) / ( DBLE( N )*EPS )
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DBDT05
*
      END
