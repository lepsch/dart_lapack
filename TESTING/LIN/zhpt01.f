      SUBROUTINE ZHPT01( UPLO, N, A, AFAC, IPIV, C, LDC, RWORK, RESID )

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
      double             RWORK( * );
      COMPLEX*16         A( * ), AFAC( * ), C( LDC, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, JC;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANHE, ZLANHP;
      // EXTERNAL LSAME, DLAMCH, ZLANHE, ZLANHP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASET, ZLAVHP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DIMAG
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Determine EPS and the norm of A.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZLANHP( '1', UPLO, N, A, RWORK )

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      JC = 1
      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 10
            if ( DIMAG( AFAC( JC ) ).NE.ZERO ) {
               RESID = ONE / EPS
               RETURN
            }
            JC = JC + J + 1
   10    CONTINUE
      } else {
         for (J = 1; J <= N; J++) { // 20
            if ( DIMAG( AFAC( JC ) ).NE.ZERO ) {
               RESID = ONE / EPS
               RETURN
            }
            JC = JC + N - J + 1
   20    CONTINUE
      }

      // Initialize C to the identity matrix.

      zlaset('Full', N, N, CZERO, CONE, C, LDC );

      // Call ZLAVHP to form the product D * U' (or D * L' ).

      zlavhp(UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, IPIV, C, LDC, INFO );

      // Call ZLAVHP again to multiply by U ( or L ).

      zlavhp(UPLO, 'No transpose', 'Unit', N, N, AFAC, IPIV, C, LDC, INFO );

      // Compute the difference  C - A .

      if ( LSAME( UPLO, 'U' ) ) {
         JC = 0
         for (J = 1; J <= N; J++) { // 40
            DO 30 I = 1, J - 1
               C( I, J ) = C( I, J ) - A( JC+I )
   30       CONTINUE
            C( J, J ) = C( J, J ) - DBLE( A( JC+J ) )
            JC = JC + J
   40    CONTINUE
      } else {
         JC = 1
         for (J = 1; J <= N; J++) { // 60
            C( J, J ) = C( J, J ) - DBLE( A( JC ) )
            DO 50 I = J + 1, N
               C( I, J ) = C( I, J ) - A( JC+I-J )
   50       CONTINUE
            JC = JC + N - J + 1
   60    CONTINUE
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = ZLANHE( '1', UPLO, N, C, LDC, RWORK )

      if ( ANORM.LE.ZERO ) {
         IF( RESID.NE.ZERO ) RESID = ONE / EPS
      } else {
         RESID = ( ( RESID / DBLE( N ) ) / ANORM ) / EPS
      }

      RETURN

      // End of ZHPT01

      }
