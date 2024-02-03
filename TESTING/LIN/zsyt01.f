      void zsyt01(UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX*16         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANSY;
      // EXTERNAL LSAME, DLAMCH, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASET, ZLAVSY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Determine EPS and the norm of A.

      EPS = DLAMCH( 'Epsilon' );
      ANORM = ZLANSY( '1', UPLO, N, A, LDA, RWORK );

      // Initialize C to the identity matrix.

      zlaset('Full', N, N, CZERO, CONE, C, LDC );

      // Call ZLAVSY to form the product D * U' (or D * L' ).

      zlavsy(UPLO, 'Transpose', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Call ZLAVSY again to multiply by U (or L ).

      zlavsy(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Compute the difference  C - A .

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               C( I, J ) = C( I, J ) - A( I, J );
            } // 10
         } // 20
      } else {
         for (J = 1; J <= N; J++) { // 40
            for (I = J; I <= N; I++) { // 30
               C( I, J ) = C( I, J ) - A( I, J );
            } // 30
         } // 40
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = ZLANSY( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS;
      }

      return;

      // End of ZSYT01

      }
