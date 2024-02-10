      void zgbt02(TRANS, M, N, KL, KU, NRHS, final Matrix<double> A, final int LDA, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, final Array<double> RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                KL, KU, LDA, LDB, LDX, M, N, NRHS;
      double             RESID;
      double             RWORK( * );
      Complex         A( LDA, * ), B( LDB, * ), X( LDX, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      int                I1, I2, J, KD, N1;
      double             ANORM, BNORM, EPS, TEMP, XNORM;
      Complex         ZDUM;
      // ..
      // .. External Functions ..
      //- bool               DISNAN, lsame;
      //- double             DLAMCH, DZASUM;
      // EXTERNAL DISNAN, DLAMCH, DZASUM, lsame
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
      CABS1[ZDUM] = ( ZDUM.toDouble() ).abs() + ( DIMAG( ZDUM ) ).abs();

      // Quick return if N = 0 pr NRHS = 0

      if ( M <= 0 || N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = dlamch( 'Epsilon' );
      ANORM = ZERO;
      if ( lsame( TRANS, 'N' ) ) {

         // Find norm1(A).

         KD = KU + 1;
         for (J = 1; J <= N; J++) { // 10
            I1 = max( KD+1-J, 1 );
            I2 = min( KD+M-J, KL+KD );
            if ( I2 >= I1 ) {
               TEMP = DZASUM( I2-I1+1, A( I1, J ), 1 );
               if( ANORM < TEMP || disnan( TEMP ) ) ANORM = TEMP;
            }
         } // 10
      } else {

         // Find normI(A).

         for (I1 = 1; I1 <= M; I1++) { // 12
            RWORK[I1] = ZERO;
         } // 12
         for (J = 1; J <= N; J++) { // 16
            KD = KU + 1 - J;
            for (I1 = max( 1, J-KU ); I1 <= min( M, J+KL ); I1++) { // 14
               RWORK[I1] = RWORK( I1 ) + CABS1( A( KD+I1, J ) );
            } // 14
         } // 16
         for (I1 = 1; I1 <= M; I1++) { // 18
            TEMP = RWORK( I1 );
            if( ANORM < TEMP || disnan( TEMP ) ) ANORM = TEMP;
         } // 18
      }
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      if ( lsame( TRANS, 'T' ) || lsame( TRANS, 'C' ) ) {
         N1 = N;
      } else {
         N1 = M;
      }

      // Compute B - op(A)*X

      for (J = 1; J <= NRHS; J++) { // 20
         zgbmv(TRANS, M, N, KL, KU, -CONE, A, LDA, X( 1, J ), 1, CONE, B( 1, J ), 1 );
      } // 20

      // Compute the maximum over the number of right hand sides of
      //    norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 30
         BNORM = DZASUM( N1, B( 1, J ), 1 );
         XNORM = DZASUM( N1, X( 1, J ), 1 );
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM / ANORM ) / XNORM ) / EPS );
         }
      } // 30

      }
