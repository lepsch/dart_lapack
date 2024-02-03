      SUBROUTINE ZGBT02( TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                KL, KU, LDA, LDB, LDX, M, N, NRHS;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I1, I2, J, KD, N1;
      double             ANORM, BNORM, EPS, TEMP, XNORM;
      COMPLEX*16         ZDUM
      // ..
      // .. External Functions ..
      bool               DISNAN, LSAME;
      double             DLAMCH, DZASUM;
      // EXTERNAL DISNAN, DLAMCH, DZASUM, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGBMV
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, MIN
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Quick return if N = 0 pr NRHS = 0

      if ( M.LE.0 .OR. N.LE.0 .OR. NRHS.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' )
      ANORM = ZERO
      if ( LSAME( TRANS, 'N' ) ) {

         // Find norm1(A).

         KD = KU + 1
         for (J = 1; J <= N; J++) { // 10
            I1 = MAX( KD+1-J, 1 )
            I2 = MIN( KD+M-J, KL+KD )
            if ( I2.GE.I1 ) {
               TEMP = DZASUM( I2-I1+1, A( I1, J ), 1 )
               IF( ANORM.LT.TEMP .OR. DISNAN( TEMP ) ) ANORM = TEMP
            }
         } // 10
      } else {

         // Find normI(A).

         for (I1 = 1; I1 <= M; I1++) { // 12
            RWORK( I1 ) = ZERO
         } // 12
         for (J = 1; J <= N; J++) { // 16
            KD = KU + 1 - J
            DO 14 I1 = MAX( 1, J-KU ), MIN( M, J+KL )
               RWORK( I1 ) = RWORK( I1 ) + CABS1( A( KD+I1, J ) )
            } // 14
         } // 16
         for (I1 = 1; I1 <= M; I1++) { // 18
            TEMP = RWORK( I1 )
            IF( ANORM.LT.TEMP .OR. DISNAN( TEMP ) ) ANORM = TEMP
         } // 18
      }
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      if ( LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' ) ) {
         N1 = N
      } else {
         N1 = M
      }

      // Compute B - op(A)*X

      for (J = 1; J <= NRHS; J++) { // 20
         zgbmv(TRANS, M, N, KL, KU, -CONE, A, LDA, X( 1, J ), 1, CONE, B( 1, J ), 1 );
      } // 20

      // Compute the maximum over the number of right hand sides of
         // norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).

      RESID = ZERO
      for (J = 1; J <= NRHS; J++) { // 30
         BNORM = DZASUM( N1, B( 1, J ), 1 )
         XNORM = DZASUM( N1, X( 1, J ), 1 )
         if ( XNORM.LE.ZERO ) {
            RESID = ONE / EPS
         } else {
            RESID = MAX( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS )
         }
      } // 30

      RETURN

      // End of ZGBT02

      }
