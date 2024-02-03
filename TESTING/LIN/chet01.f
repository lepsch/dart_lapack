      SUBROUTINE CHET01( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      REAL               ANORM, EPS
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHE, SLAMCH
      // EXTERNAL LSAME, CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLAVHE, CLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, REAL
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
      ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK )
*
      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.
*
      DO 10 J = 1, N
         IF( AIMAG( AFAC( J, J ) ).NE.ZERO ) THEN
            RESID = ONE / EPS
            RETURN
         END IF
   10 CONTINUE
*
      // Initialize C to the identity matrix.
*
      CALL CLASET( 'Full', N, N, CZERO, CONE, C, LDC )
*
      // Call CLAVHE to form the product D * U' (or D * L' ).
*
      CALL CLAVHE( UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO )
*
      // Call CLAVHE again to multiply by U (or L ).
*
      CALL CLAVHE( UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO )
*
      // Compute the difference  C - A .
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 30 J = 1, N
            DO 20 I = 1, J - 1
               C( I, J ) = C( I, J ) - A( I, J )
   20       CONTINUE
            C( J, J ) = C( J, J ) - REAL( A( J, J ) )
   30    CONTINUE
      ELSE
         DO 50 J = 1, N
            C( J, J ) = C( J, J ) - REAL( A( J, J ) )
            DO 40 I = J + 1, N
               C( I, J ) = C( I, J ) - A( I, J )
   40       CONTINUE
   50    CONTINUE
      END IF
*
      // Compute norm( C - A ) / ( N * norm(A) * EPS )
*
      RESID = CLANHE( '1', UPLO, N, C, LDC, RWORK )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
      END IF
*
      RETURN
*
      // End of CHET01
*
      END
