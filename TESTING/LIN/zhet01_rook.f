      SUBROUTINE ZHET01_ROOK( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID )

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
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             ZLANHE, DLAMCH;
      // EXTERNAL LSAME, ZLANHE, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASET, ZLAVHE_ROOK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DIMAG, DBLE
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Determine EPS and the norm of A.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      for (J = 1; J <= N; J++) { // 10
         if ( DIMAG( AFAC( J, J ) ).NE.ZERO ) {
            RESID = ONE / EPS
            RETURN
         }
      } // 10

      // Initialize C to the identity matrix.

      zlaset('Full', N, N, CZERO, CONE, C, LDC );

      // Call ZLAVHE_ROOK to form the product D * U' (or D * L' ).

      zlavhe_rook(UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Call ZLAVHE_ROOK again to multiply by U (or L ).

      zlavhe_rook(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Compute the difference  C - A .

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 30
            for (I = 1; I <= J - 1; I++) { // 20
               C( I, J ) = C( I, J ) - A( I, J )
            } // 20
            C( J, J ) = C( J, J ) - DBLE( A( J, J ) )
         } // 30
      } else {
         for (J = 1; J <= N; J++) { // 50
            C( J, J ) = C( J, J ) - DBLE( A( J, J ) )
            for (I = J + 1; I <= N; I++) { // 40
               C( I, J ) = C( I, J ) - A( I, J )
            } // 40
         } // 50
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = ZLANHE( '1', UPLO, N, C, LDC, RWORK )

      if ( ANORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         RESID = ( ( RESID/DBLE( N ) )/ANORM ) / EPS
      }

      RETURN

      // End of ZHET01_ROOK

      }
