      void cgbt02(final int TRANS, final int M, final int N, final int KL, final int KU, final int NRHS, final Matrix<double> A_, final int LDA, final Matrix<double> X_, final int LDX, final Matrix<double> B_, final int LDB, final Array<double> RWORK_, final int RESID,) {
  final A = A_.dim();
  final X = X_.dim();
  final B = B_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                KL, KU, LDA, LDB, LDX, M, N, NRHS;
      double               RESID;
      double               RWORK( * );
      Complex            A( LDA, * ), B( LDB, * ), X( LDX, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      int                I1, I2, J, KD, N1;
      double               ANORM, BNORM, EPS, TEMP, XNORM;
      Complex            CDUM;
      // ..
      // .. External Functions ..
      //- bool               lsame, SISNAN;
      //- REAL               SCASUM, SLAMCH;
      // EXTERNAL lsame, SCASUM, SISNAN, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGBMV
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
      // ..
      // .. Statement Function definitions ..
      CABS1[CDUM] = ( double( CDUM ) ).abs() + ( AIMAG( CDUM ) ).abs();

      // Quick return if N = 0 pr NRHS = 0

      if ( M <= 0 || N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = ZERO;
      if ( lsame( TRANS, 'N' ) ) {

         // Find norm1(A).

         KD = KU + 1;
         for (J = 1; J <= N; J++) { // 10
            I1 = max( KD+1-J, 1 );
            I2 = min( KD+M-J, KL+KD );
            if ( I2 >= I1 ) {
               TEMP = SCASUM( I2-I1+1, A( I1, J ), 1 );
               if( ANORM < TEMP || SISNAN( TEMP ) ) ANORM = TEMP;
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
            if( ANORM < TEMP || SISNAN( TEMP ) ) ANORM = TEMP;
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
         cgbmv(TRANS, M, N, KL, KU, -CONE, A, LDA, X( 1, J ), 1, CONE, B( 1, J ), 1 );
      } // 20

      // Compute the maximum over the number of right hand sides of
      //    norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).

      RESID = ZERO;
      for (J = 1; J <= NRHS; J++) { // 30
         BNORM = SCASUM( N1, B( 1, J ), 1 );
         XNORM = SCASUM( N1, X( 1, J ), 1 );
         if ( XNORM <= ZERO ) {
            RESID = ONE / EPS;
         } else {
            RESID = max( RESID, ( ( BNORM/ANORM )/XNORM )/EPS );
         }
      } // 30

      }
