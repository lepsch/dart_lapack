      SUBROUTINE CGET10( M, N, A, LDA, B, LDB, WORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDB, M, N
      REAL               RESULT
*     ..
*     .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      int                J
      REAL               ANORM, EPS, UNFL, WNORM
*     ..
*     .. External Functions ..
      REAL               SCASUM, SLAMCH, CLANGE
      EXTERNAL           SCASUM, SLAMCH, CLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CCOPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         RESULT = ZERO
         RETURN
      END IF
*
      UNFL = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
*
      WNORM = ZERO
      DO 10 J = 1, N
         CALL CCOPY( M, A( 1, J ), 1, WORK, 1 )
         CALL CAXPY( M, CMPLX( -ONE ), B( 1, J ), 1, WORK, 1 )
         WNORM = MAX( WNORM, SCASUM( N, WORK, 1 ) )
   10 CONTINUE
*
      ANORM = MAX( CLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
*
      IF( ANORM.GT.WNORM ) THEN
         RESULT = ( WNORM / ANORM ) / ( M*EPS )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*EPS )
         ELSE
            RESULT = MIN( WNORM / ANORM, REAL( M ) ) / ( M*EPS )
         END IF
      END IF
*
      RETURN
*
*     End of CGET10
*
      END
