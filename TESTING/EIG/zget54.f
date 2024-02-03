      SUBROUTINE ZGET54( N, A, LDA, B, LDB, S, LDS, T, LDT, U, LDU, V, LDV, WORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDB, LDS, LDT, LDU, LDV, N
      DOUBLE PRECISION   RESULT
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), S( LDS, * ), T( LDT, * ), U( LDU, * ), V( LDV, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   ABNORM, ULP, UNFL, WNORM
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMM, ZLACPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      RESULT = ZERO
      IF( N.LE.0 ) RETURN
*
*     Constants
*
      UNFL = DLAMCH( 'Safe minimum' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
*
*     compute the norm of (A,B)
*
      CALL ZLACPY( 'Full', N, N, A, LDA, WORK, N )
      CALL ZLACPY( 'Full', N, N, B, LDB, WORK( N*N+1 ), N )
      ABNORM = MAX( ZLANGE( '1', N, 2*N, WORK, N, DUM ), UNFL )
*
*     Compute W1 = A - U*S*V', and put in the array WORK(1:N*N)
*
      CALL ZLACPY( ' ', N, N, A, LDA, WORK, N )
      CALL ZGEMM( 'N', 'N', N, N, N, CONE, U, LDU, S, LDS, CZERO, WORK( N*N+1 ), N )
*
      CALL ZGEMM( 'N', 'C', N, N, N, -CONE, WORK( N*N+1 ), N, V, LDV, CONE, WORK, N )
*
*     Compute W2 = B - U*T*V', and put in the workarray W(N*N+1:2*N*N)
*
      CALL ZLACPY( ' ', N, N, B, LDB, WORK( N*N+1 ), N )
      CALL ZGEMM( 'N', 'N', N, N, N, CONE, U, LDU, T, LDT, CZERO, WORK( 2*N*N+1 ), N )
*
      CALL ZGEMM( 'N', 'C', N, N, N, -CONE, WORK( 2*N*N+1 ), N, V, LDV, CONE, WORK( N*N+1 ), N )
*
*     Compute norm(W)/ ( ulp*norm((A,B)) )
*
      WNORM = ZLANGE( '1', N, 2*N, WORK, N, DUM )
*
      IF( ABNORM.GT.WNORM ) THEN
         RESULT = ( WNORM / ABNORM ) / ( 2*N*ULP )
      ELSE
         IF( ABNORM.LT.ONE ) THEN
            RESULT = ( MIN( WNORM, 2*N*ABNORM ) / ABNORM ) / ( 2*N*ULP )
         ELSE
            RESULT = MIN( WNORM / ABNORM, DBLE( 2*N ) ) / ( 2*N*ULP )
         END IF
      END IF
*
      RETURN
*
*     End of ZGET54
*
      END
