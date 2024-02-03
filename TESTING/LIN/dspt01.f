      SUBROUTINE DSPT01( UPLO, N, A, AFAC, IPIV, C, LDC, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDC, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             A( * ), AFAC( * ), C( LDC, * ), RWORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, JC;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANSP, DLANSY;
      // EXTERNAL LSAME, DLAMCH, DLANSP, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET, DLAVSP
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

      // Determine EPS and the norm of A.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = DLANSP( '1', UPLO, N, A, RWORK )

      // Initialize C to the identity matrix.

      CALL DLASET( 'Full', N, N, ZERO, ONE, C, LDC )

      // Call DLAVSP to form the product D * U' (or D * L' ).

      CALL DLAVSP( UPLO, 'Transpose', 'Non-unit', N, N, AFAC, IPIV, C, LDC, INFO )

      // Call DLAVSP again to multiply by U ( or L ).

      CALL DLAVSP( UPLO, 'No transpose', 'Unit', N, N, AFAC, IPIV, C, LDC, INFO )

      // Compute the difference  C - A .

      IF( LSAME( UPLO, 'U' ) ) THEN
         JC = 0
         DO 20 J = 1, N
            DO 10 I = 1, J
               C( I, J ) = C( I, J ) - A( JC+I )
   10       CONTINUE
            JC = JC + J
   20    CONTINUE
      } else {
         JC = 1
         DO 40 J = 1, N
            DO 30 I = J, N
               C( I, J ) = C( I, J ) - A( JC+I-J )
   30       CONTINUE
            JC = JC + N - J + 1
   40    CONTINUE
      END IF

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = DLANSY( '1', UPLO, N, C, LDC, RWORK )

      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
      END IF

      RETURN

      // End of DSPT01

      }
