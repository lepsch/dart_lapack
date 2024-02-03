      SUBROUTINE SSPT01( UPLO, N, A, AFAC, IPIV, C, LDC, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                LDC, N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( * ), AFAC( * ), C( LDC, * ), RWORK( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, JC;
      REAL               ANORM, EPS
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANSP, SLANSY
      // EXTERNAL LSAME, SLAMCH, SLANSP, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAVSP, SLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..
*
      // Quick exit if N = 0.
*
      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
      // Determine EPS and the norm of A.
*
      EPS = SLAMCH( 'Epsilon' )
      ANORM = SLANSP( '1', UPLO, N, A, RWORK )
*
      // Initialize C to the identity matrix.
*
      CALL SLASET( 'Full', N, N, ZERO, ONE, C, LDC )
*
      // Call SLAVSP to form the product D * U' (or D * L' ).
*
      CALL SLAVSP( UPLO, 'Transpose', 'Non-unit', N, N, AFAC, IPIV, C, LDC, INFO )
*
      // Call SLAVSP again to multiply by U ( or L ).
*
      CALL SLAVSP( UPLO, 'No transpose', 'Unit', N, N, AFAC, IPIV, C, LDC, INFO )
*
      // Compute the difference  C - A .
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         JC = 0
         DO 20 J = 1, N
            DO 10 I = 1, J
               C( I, J ) = C( I, J ) - A( JC+I )
   10       CONTINUE
            JC = JC + J
   20    CONTINUE
      ELSE
         JC = 1
         DO 40 J = 1, N
            DO 30 I = J, N
               C( I, J ) = C( I, J ) - A( JC+I-J )
   30       CONTINUE
            JC = JC + N - J + 1
   40    CONTINUE
      END IF
*
      // Compute norm( C - A ) / ( N * norm(A) * EPS )
*
      RESID = SLANSY( '1', UPLO, N, C, LDC, RWORK )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
      END IF
*
      RETURN
*
      // End of SSPT01
*
      END
