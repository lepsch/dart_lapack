      SUBROUTINE SSYT01_ROOK( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID )

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
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
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
      // EXTERNAL SLASET, SLAVSY_ROOK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK )

      // Initialize C to the identity matrix.

      slaset('Full', N, N, ZERO, ONE, C, LDC );

      // Call SLAVSY_ROOK to form the product D * U' (or D * L' ).

      slavsy_rook(UPLO, 'Transpose', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Call SLAVSY_ROOK again to multiply by U (or L ).

      slavsy_rook(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Compute the difference  C - A .

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               C( I, J ) = C( I, J ) - A( I, J )
            } // 10
         } // 20
      } else {
         for (J = 1; J <= N; J++) { // 40
            for (I = J; I <= N; I++) { // 30
               C( I, J ) = C( I, J ) - A( I, J )
            } // 30
         } // 40
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = SLANSY( '1', UPLO, N, C, LDC, RWORK )

      if ( ANORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS
      }

      RETURN

      // End of SSYT01_ROOK

      }
