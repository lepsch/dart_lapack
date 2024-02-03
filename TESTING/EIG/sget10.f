      SUBROUTINE SGET10( M, N, A, LDA, B, LDB, WORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, M, N;
      REAL               RESULT
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                J;
      REAL               ANORM, EPS, UNFL, WNORM
      // ..
      // .. External Functions ..
      REAL               SASUM, SLAMCH, SLANGE
      // EXTERNAL SASUM, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( M.LE.0 .OR. N.LE.0 ) THEN
         RESULT = ZERO
         RETURN
      END IF

      UNFL = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )

      WNORM = ZERO
      DO 10 J = 1, N
         CALL SCOPY( M, A( 1, J ), 1, WORK, 1 )
         CALL SAXPY( M, -ONE, B( 1, J ), 1, WORK, 1 )
         WNORM = MAX( WNORM, SASUM( N, WORK, 1 ) )
   10 CONTINUE

      ANORM = MAX( SLANGE( '1', M, N, A, LDA, WORK ), UNFL )

      IF( ANORM.GT.WNORM ) THEN
         RESULT = ( WNORM / ANORM ) / ( M*EPS )
      ELSE
         IF( ANORM.LT.ONE ) THEN
            RESULT = ( MIN( WNORM, M*ANORM ) / ANORM ) / ( M*EPS )
         ELSE
            RESULT = MIN( WNORM / ANORM, REAL( M ) ) / ( M*EPS )
         END IF
      END IF

      RETURN

      // End of SGET10

      }
