      void cgbt02(TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDX, B, LDB, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                KL, KU, LDA, LDB, LDX, M, N, NRHS;
      REAL               RESID;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * );
      COMPLEX            A( LDA, * ), B( LDB, * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I1, I2, J, KD, N1;
      REAL               ANORM, BNORM, EPS, TEMP, XNORM;
      COMPLEX            CDUM;
      // ..
      // .. External Functions ..
      //- bool               LSAME, SISNAN;
      //- REAL               SCASUM, SLAMCH;
      // EXTERNAL LSAME, SCASUM, SISNAN, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGBMV
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
      // ..
      // .. Statement Function definitions ..
      CABS1( CDUM ) = ( REAL( CDUM ) ).abs() + ( AIMAG( CDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      // Quick return if N = 0 pr NRHS = 0

      if ( M <= 0 || N <= 0 || NRHS <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = ZERO;
      if ( LSAME( TRANS, 'N' ) ) {

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
            RWORK( I1 ) = ZERO;
         } // 12
         for (J = 1; J <= N; J++) { // 16
            KD = KU + 1 - J;
            for (I1 = max( 1, J-KU ); I1 <= min( M, J+KL ); I1++) { // 14
               RWORK( I1 ) = RWORK( I1 ) + CABS1( A( KD+I1, J ) );
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

      if ( LSAME( TRANS, 'T' ) || LSAME( TRANS, 'C' ) ) {
         N1 = N;
      } else {
         N1 = M;
      }

      // Compute B - op(A)*X

      for (J = 1; J <= NRHS; J++) { // 20
         cgbmv(TRANS, M, N, KL, KU, -CONE, A, LDA, X( 1, J ), 1, CONE, B( 1, J ), 1 );
      } // 20

      // Compute the maximum over the number of right hand sides of
         // norm(B - op(A)*X) / ( norm(op(A)) * norm(X) * EPS ).

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

      return;
      }
