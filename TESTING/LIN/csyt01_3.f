      SUBROUTINE CSYT01_3( UPLO, N, A, LDA, AFAC, LDAFAC, E, IPIV, C, LDC, RWORK, RESID )
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
      COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), E( * )
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
      REAL               SLAMCH, CLANSY
      // EXTERNAL LSAME, SLAMCH, CLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASET, CLAVSY_ROOK, CSYCONVF_ROOK
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
      // a) Revert to multipliers of L
*
      CALL CSYCONVF_ROOK( UPLO, 'R', N, AFAC, LDAFAC, E, IPIV, INFO )
*
      // 1) Determine EPS and the norm of A.
*
      EPS = SLAMCH( 'Epsilon' )
      ANORM = CLANSY( '1', UPLO, N, A, LDA, RWORK )
*
      // 2) Initialize C to the identity matrix.
*
      CALL CLASET( 'Full', N, N, CZERO, CONE, C, LDC )
*
      // 3) Call ZLAVSY_ROOK to form the product D * U' (or D * L' ).
*
      CALL CLAVSY_ROOK( UPLO, 'Transpose', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO )
*
      // 4) Call ZLAVSY_ROOK again to multiply by U (or L ).
*
      CALL CLAVSY_ROOK( UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO )
*
      // 5) Compute the difference  C - A .
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 1, N
            DO I = 1, J
               C( I, J ) = C( I, J ) - A( I, J )
            END DO
         END DO
      ELSE
         DO J = 1, N
            DO I = J, N
               C( I, J ) = C( I, J ) - A( I, J )
            END DO
         END DO
      END IF
*
      // 6) Compute norm( C - A ) / ( N * norm(A) * EPS )
*
      RESID = CLANSY( '1', UPLO, N, C, LDC, RWORK )
*
      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
      END IF

*
      // b) Convert to factor of L (or U)
*
      CALL CSYCONVF_ROOK( UPLO, 'C', N, AFAC, LDAFAC, E, IPIV, INFO )
*
      RETURN
*
      // End of CSYT01_3
*
      END
