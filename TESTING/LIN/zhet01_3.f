      SUBROUTINE ZHET01_3( UPLO, N, A, LDA, AFAC, LDAFAC, E, IPIV, C, LDC, RWORK, RESID )

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
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), E( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
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
      // EXTERNAL ZLASET, ZLAVHE_ROOK, ZSYCONVF_ROOK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DIMAG, DBLE
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO
         RETURN
      }

      // a) Revert to multipliers of L

      zsyconvf_rook(UPLO, 'R', N, AFAC, LDAFAC, E, IPIV, INFO );

      // 1) Determine EPS and the norm of A.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      for (J = 1; J <= N; J++) {
         if ( DIMAG( AFAC( J, J ) ) != ZERO ) {
            RESID = ONE / EPS
            RETURN
         }
      }

      // 2) Initialize C to the identity matrix.

      zlaset('Full', N, N, CZERO, CONE, C, LDC );

      // 3) Call ZLAVHE_ROOK to form the product D * U' (or D * L' ).

      zlavhe_rook(UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // 4) Call ZLAVHE_RK again to multiply by U (or L ).

      zlavhe_rook(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // 5) Compute the difference  C - A .

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J - 1; I++) {
               C( I, J ) = C( I, J ) - A( I, J )
            }
            C( J, J ) = C( J, J ) - DBLE( A( J, J ) )
         }
      } else {
         for (J = 1; J <= N; J++) {
            C( J, J ) = C( J, J ) - DBLE( A( J, J ) )
            for (I = J + 1; I <= N; I++) {
               C( I, J ) = C( I, J ) - A( I, J )
            }
         }
      }

      // 6) Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = ZLANHE( '1', UPLO, N, C, LDC, RWORK )

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID/DBLE( N ) )/ANORM ) / EPS
      }

      // b) Convert to factor of L (or U)

      zsyconvf_rook(UPLO, 'C', N, AFAC, LDAFAC, E, IPIV, INFO );

      RETURN

      // End of ZHET01_3

      }
