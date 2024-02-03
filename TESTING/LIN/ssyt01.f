      SUBROUTINE SSYT01( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), RWORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      REAL               ANORM, EPS
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANSY
      // EXTERNAL LSAME, SLAMCH, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASET, SLAVSY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF

      // Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK )

      // Initialize C to the identity matrix.

      CALL SLASET( 'Full', N, N, ZERO, ONE, C, LDC )

      // Call SLAVSY to form the product D * U' (or D * L' ).

      CALL SLAVSY( UPLO, 'Transpose', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO )

      // Call SLAVSY again to multiply by U (or L ).

      CALL SLAVSY( UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO )

      // Compute the difference  C - A .

      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, J
               C( I, J ) = C( I, J ) - A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE
         DO 40 J = 1, N
            DO 30 I = J, N
               C( I, J ) = C( I, J ) - A( I, J )
   30       CONTINUE
   40    CONTINUE
      END IF

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = SLANSY( '1', UPLO, N, C, LDC, RWORK )

      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
      END IF

      RETURN

      // End of SSYT01

      END
