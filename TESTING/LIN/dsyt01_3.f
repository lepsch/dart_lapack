      SUBROUTINE DSYT01_3( UPLO, N, A, LDA, AFAC, LDAFAC, E, IPIV, C, LDC, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), E( * ), RWORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANSY;
      // EXTERNAL LSAME, DLAMCH, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET, DLAVSY_ROOK, DSYCONVF_ROOK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      IF( N.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF

      // a) Revert to multipliers of L

      CALL DSYCONVF_ROOK( UPLO, 'R', N, AFAC, LDAFAC, E, IPIV, INFO )

      // 1) Determine EPS and the norm of A.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = DLANSY( '1', UPLO, N, A, LDA, RWORK )

      // 2) Initialize C to the identity matrix.

      CALL DLASET( 'Full', N, N, ZERO, ONE, C, LDC )

      // 3) Call DLAVSY_ROOK to form the product D * U' (or D * L' ).

      CALL DLAVSY_ROOK( UPLO, 'Transpose', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO )

      // 4) Call DLAVSY_ROOK again to multiply by U (or L ).

      CALL DLAVSY_ROOK( UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO )

      // 5) Compute the difference  C - A.

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

      // 6) Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = DLANSY( '1', UPLO, N, C, LDC, RWORK )

      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      ELSE
         RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
      END IF


      // b) Convert to factor of L (or U)

      CALL DSYCONVF_ROOK( UPLO, 'C', N, AFAC, LDAFAC, E, IPIV, INFO )

      RETURN

      // End of DSYT01_3

      END
