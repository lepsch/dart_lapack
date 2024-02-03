      SUBROUTINE CBDT05( M, N, A, LDA, S, NS, U, LDU, VT, LDVT, WORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDU, LDVT, M, N, NS;
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               S( * )
      COMPLEX            A( LDA, * ), U( * ), VT( LDVT, * ), WORK( * )
*     ..
*
* ======================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      int                I, J;
      REAL               ANORM, EPS
*     ..
*     .. Local Arrays ..
      REAL               DUM( 1 )
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ISAMAX;
      REAL               SASUM, SCASUM, SLAMCH, CLANGE
      // EXTERNAL LSAME, ISAMAX, SASUM, SCASUM, SLAMCH, CLANGE
*     ..
*     .. External Subroutines ..
      // EXTERNAL CGEMM
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      RESID = ZERO
      IF( MIN( M, N ).LE.0 .OR. NS.LE.0 ) RETURN
*
      EPS = SLAMCH( 'Precision' )
      ANORM = CLANGE( 'M', M, N, A, LDA, DUM )
*
*     Compute U' * A * V.
*
      CALL CGEMM( 'N', 'C', M, NS, N, CONE, A, LDA, VT, LDVT, CZERO, WORK( 1+NS*NS ), M )       CALL CGEMM( 'C', 'N', NS, NS, M, -CONE, U, LDU, WORK( 1+NS*NS ), M, CZERO, WORK, NS )
*
*     norm(S - U' * B * V)
*
      J = 0
      DO 10 I = 1, NS
         WORK( J+I ) =  WORK( J+I ) + CMPLX( S( I ), ZERO )
         RESID = MAX( RESID, SCASUM( NS, WORK( J+1 ), 1 ) )
         J = J + NS
   10 CONTINUE
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         IF( ANORM.GE.RESID ) THEN
            RESID = ( RESID / ANORM ) / ( REAL( N )*EPS )
         ELSE
            IF( ANORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, REAL( N )*ANORM ) / ANORM ) / ( REAL( N )*EPS )
            ELSE
               RESID = MIN( RESID / ANORM, REAL( N ) ) / ( REAL( N )*EPS )
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of CBDT05
*
      END
