      SUBROUTINE ZBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RWORK,
     $                   RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDB, LDC, LDU, M, N
      DOUBLE PRECISION   RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         B( LDB, * ), C( LDC, * ), U( LDU, * ),
     $                   WORK( * )
*     ..
*
* ======================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J
      DOUBLE PRECISION   BNORM, EPS, REALMN
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DZASUM, ZLANGE
      EXTERNAL           DLAMCH, DZASUM, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      RESID = ZERO
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
      REALMN = DBLE( MAX( M, N ) )
      EPS = DLAMCH( 'Precision' )
*
*     Compute norm(B - U * C)
*
      DO 10 J = 1, N
         CALL ZCOPY( M, B( 1, J ), 1, WORK, 1 )
         CALL ZGEMV( 'No transpose', M, M, -DCMPLX( ONE ), U, LDU,
     $               C( 1, J ), 1, DCMPLX( ONE ), WORK, 1 )
         RESID = MAX( RESID, DZASUM( M, WORK, 1 ) )
   10 CONTINUE
*
*     Compute norm of B.
*
      BNORM = ZLANGE( '1', M, N, B, LDB, RWORK )
*
      IF( BNORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO )
     $      RESID = ONE / EPS
      ELSE
         IF( BNORM.GE.RESID ) THEN
            RESID = ( RESID / BNORM ) / ( REALMN*EPS )
         ELSE
            IF( BNORM.LT.ONE ) THEN
               RESID = ( MIN( RESID, REALMN*BNORM ) / BNORM ) /
     $                 ( REALMN*EPS )
            ELSE
               RESID = MIN( RESID / BNORM, REALMN ) / ( REALMN*EPS )
            END IF
         END IF
      END IF
      RETURN
*
*     End of ZBDT02
*
      END