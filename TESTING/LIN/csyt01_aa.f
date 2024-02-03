      SUBROUTINE CSYT01_AA( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID )

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
      COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * )
      REAL               RWORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE  = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               ANORM, EPS
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, CLANSY
      // EXTERNAL LSAME, SLAMCH, CLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASET, CLAVSY
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
      ANORM = CLANSY( '1', UPLO, N, A, LDA, RWORK )

      // Initialize C to the tridiagonal matrix T.

      CALL CLASET( 'Full', N, N, CZERO, CZERO, C, LDC )
      CALL CLACPY( 'F', 1, N, AFAC( 1, 1 ), LDAFAC+1, C( 1, 1 ), LDC+1 )
      IF( N.GT.1 ) THEN
         IF( LSAME( UPLO, 'U' ) ) THEN
            CALL CLACPY( 'F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 1, 2 ), LDC+1 )             CALL CLACPY( 'F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 2, 1 ), LDC+1 )
         } else {
            CALL CLACPY( 'F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 1, 2 ), LDC+1 )             CALL CLACPY( 'F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 2, 1 ), LDC+1 )
         ENDIF

         // Call CTRMM to form the product U' * D (or L * D ).

         IF( LSAME( UPLO, 'U' ) ) THEN
            CALL CTRMM( 'Left', UPLO, 'Transpose', 'Unit', N-1, N, CONE, AFAC( 1, 2 ), LDAFAC, C( 2, 1 ), LDC )
         } else {
            CALL CTRMM( 'Left', UPLO, 'No transpose', 'Unit', N-1, N, CONE, AFAC( 2, 1 ), LDAFAC, C( 2, 1 ), LDC )
         END IF

         // Call CTRMM again to multiply by U (or L ).

         IF( LSAME( UPLO, 'U' ) ) THEN
            CALL CTRMM( 'Right', UPLO, 'No transpose', 'Unit', N, N-1, CONE, AFAC( 1, 2 ), LDAFAC, C( 1, 2 ), LDC )
         } else {
            CALL CTRMM( 'Right', UPLO, 'Transpose', 'Unit', N, N-1, CONE, AFAC( 2, 1 ), LDAFAC, C( 1, 2 ), LDC )
         END IF
      ENDIF

      // Apply symmetric pivots

      DO J = N, 1, -1
         I = IPIV( J )
         IF( I.NE.J ) CALL CSWAP( N, C( J, 1 ), LDC, C( I, 1 ), LDC )
      END DO
      DO J = N, 1, -1
         I = IPIV( J )
         IF( I.NE.J ) CALL CSWAP( N, C( 1, J ), 1, C( 1, I ), 1 )
      END DO


      // Compute the difference  C - A .

      IF( LSAME( UPLO, 'U' ) ) THEN
         DO J = 1, N
            DO I = 1, J
               C( I, J ) = C( I, J ) - A( I, J )
            END DO
         END DO
      } else {
         DO J = 1, N
            DO I = J, N
               C( I, J ) = C( I, J ) - A( I, J )
            END DO
         END DO
      END IF

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = CLANSY( '1', UPLO, N, C, LDC, RWORK )

      IF( ANORM.LE.ZERO ) THEN
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
      END IF

      RETURN

      // End of CSYT01_AA

      }
