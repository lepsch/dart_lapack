      SUBROUTINE ZGET10( M, N, A, LDA, B, LDB, WORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDB, M, N;
      double             RESULT;
*     ..
*     .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      int                J;
      double             ANORM, EPS, UNFL, WNORM;
*     ..
*     .. External Functions ..
      double             DLAMCH, DZASUM, ZLANGE;
      // EXTERNAL DLAMCH, DZASUM, ZLANGE
*     ..
*     .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
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
      UNFL = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
*
      WNORM = ZERO
      DO 10 J = 1, N
         CALL ZCOPY( M, A( 1, J ), 1, WORK, 1 )
         CALL ZAXPY( M, DCMPLX( -ONE ), B( 1, J ), 1, WORK, 1 )
         WNORM = MAX( WNORM, DZASUM( N, WORK, 1 ) )
   10 CONTINUE
*
      ANORM = MAX( ZLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
*
      IF( ANORM.GT.WNORM ) THEN
         RESULT = ( WNORM / ANORM ) / ( M*EPS )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*EPS )
         ELSE
            RESULT = MIN( WNORM / ANORM, DBLE( M ) ) / ( M*EPS )
         END IF
      END IF
*
      RETURN
*
*     End of ZGET10
*
      END
