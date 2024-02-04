      void cgeequ(M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      double               AMAX, COLCND, ROWCND;
      // ..
      // .. Array Arguments ..
      double               C( * ), R( * );
      Complex            A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double               BIGNUM, RCMAX, RCMIN, SMLNUM;
      Complex            ZDUM;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
      // ..
      // .. Statement Functions ..
      double               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( double( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('CGEEQU', -INFO );
         return;
      }

      // Quick return if possible

      if ( M == 0 || N == 0 ) {
         ROWCND = ONE;
         COLCND = ONE;
         AMAX = ZERO;
         return;
      }

      // Get machine constants.

      SMLNUM = SLAMCH( 'S' );
      BIGNUM = ONE / SMLNUM;

      // Compute row scale factors.

      for (I = 1; I <= M; I++) { // 10
         R[I] = ZERO;
      } // 10

      // Find the maximum element in each row.

      for (J = 1; J <= N; J++) { // 30
         for (I = 1; I <= M; I++) { // 20
            R[I] = max( R( I ), CABS1( A( I, J ) ) );
         } // 20
      } // 30

      // Find the maximum and minimum scale factors.

      RCMIN = BIGNUM;
      RCMAX = ZERO;
      for (I = 1; I <= M; I++) { // 40
         RCMAX = max( RCMAX, R( I ) );
         RCMIN = min( RCMIN, R( I ) );
      } // 40
      AMAX = RCMAX;

      if ( RCMIN == ZERO ) {

         // Find the first zero scale factor and return an error code.

         for (I = 1; I <= M; I++) { // 50
            if ( R( I ) == ZERO ) {
               INFO = I;
               return;
            }
         } // 50
      } else {

         // Invert the scale factors.

         for (I = 1; I <= M; I++) { // 60
            R[I] = ONE / min( max( R( I ), SMLNUM ), BIGNUM );
         } // 60

         // Compute ROWCND = min(R(I)) / max(R(I))

         ROWCND = max( RCMIN, SMLNUM ) / min( RCMAX, BIGNUM );
      }

      // Compute column scale factors

      for (J = 1; J <= N; J++) { // 70
         C[J] = ZERO;
      } // 70

      // Find the maximum element in each column,
      // assuming the row scaling computed above.

      for (J = 1; J <= N; J++) { // 90
         for (I = 1; I <= M; I++) { // 80
            C[J] = max( C( J ), CABS1( A( I, J ) )*R( I ) );
         } // 80
      } // 90

      // Find the maximum and minimum scale factors.

      RCMIN = BIGNUM;
      RCMAX = ZERO;
      for (J = 1; J <= N; J++) { // 100
         RCMIN = min( RCMIN, C( J ) );
         RCMAX = max( RCMAX, C( J ) );
      } // 100

      if ( RCMIN == ZERO ) {

         // Find the first zero scale factor and return an error code.

         for (J = 1; J <= N; J++) { // 110
            if ( C( J ) == ZERO ) {
               INFO = M + J;
               return;
            }
         } // 110
      } else {

         // Invert the scale factors.

         for (J = 1; J <= N; J++) { // 120
            C[J] = ONE / min( max( C( J ), SMLNUM ), BIGNUM );
         } // 120

         // Compute COLCND = min(C(J)) / max(C(J))

         COLCND = max( RCMIN, SMLNUM ) / min( RCMAX, BIGNUM );
      }

      return;
      }