      SUBROUTINE CGEEQUB( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      REAL               AMAX, COLCND, ROWCND
      // ..
      // .. Array Arguments ..
      REAL               C( * ), R( * )
      COMPLEX            A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               BIGNUM, RCMAX, RCMIN, SMLNUM, RADIX, LOGRDX
      COMPLEX            ZDUM
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, LOG, REAL, AIMAG
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }
      if ( INFO != 0 ) {
         xerbla('CGEEQUB', -INFO );
         RETURN
      }

      // Quick return if possible.

      if ( M == 0 || N == 0 ) {
         ROWCND = ONE
         COLCND = ONE
         AMAX = ZERO
         RETURN
      }

      // Get machine constants.  Assume SMLNUM is a power of the radix.

      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      RADIX = SLAMCH( 'B' )
      LOGRDX = LOG( RADIX )

      // Compute row scale factors.

      for (I = 1; I <= M; I++) { // 10
         R( I ) = ZERO
      } // 10

      // Find the maximum element in each row.

      for (J = 1; J <= N; J++) { // 30
         for (I = 1; I <= M; I++) { // 20
            R( I ) = MAX( R( I ), CABS1( A( I, J ) ) )
         } // 20
      } // 30
      for (I = 1; I <= M; I++) {
         if ( R( I ).GT.ZERO ) {
            R( I ) = RADIX**INT( LOG(R( I ) ) / LOGRDX )
         }
      }

      // Find the maximum and minimum scale factors.

      RCMIN = BIGNUM
      RCMAX = ZERO
      for (I = 1; I <= M; I++) { // 40
         RCMAX = MAX( RCMAX, R( I ) )
         RCMIN = MIN( RCMIN, R( I ) )
      } // 40
      AMAX = RCMAX

      if ( RCMIN == ZERO ) {

         // Find the first zero scale factor and return an error code.

         for (I = 1; I <= M; I++) { // 50
            if ( R( I ) == ZERO ) {
               INFO = I
               RETURN
            }
         } // 50
      } else {

         // Invert the scale factors.

         for (I = 1; I <= M; I++) { // 60
            R( I ) = ONE / MIN( MAX( R( I ), SMLNUM ), BIGNUM )
         } // 60

         // Compute ROWCND = min(R(I)) / max(R(I)).

         ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      }

      // Compute column scale factors.

      for (J = 1; J <= N; J++) { // 70
         C( J ) = ZERO
      } // 70

      // Find the maximum element in each column,
      // assuming the row scaling computed above.

      for (J = 1; J <= N; J++) { // 90
         for (I = 1; I <= M; I++) { // 80
            C( J ) = MAX( C( J ), CABS1( A( I, J ) )*R( I ) )
         } // 80
         if ( C( J ).GT.ZERO ) {
            C( J ) = RADIX**INT( LOG( C( J ) ) / LOGRDX )
         }
      } // 90

      // Find the maximum and minimum scale factors.

      RCMIN = BIGNUM
      RCMAX = ZERO
      for (J = 1; J <= N; J++) { // 100
         RCMIN = MIN( RCMIN, C( J ) )
         RCMAX = MAX( RCMAX, C( J ) )
      } // 100

      if ( RCMIN == ZERO ) {

         // Find the first zero scale factor and return an error code.

         for (J = 1; J <= N; J++) { // 110
            if ( C( J ) == ZERO ) {
               INFO = M + J
               RETURN
            }
         } // 110
      } else {

         // Invert the scale factors.

         for (J = 1; J <= N; J++) { // 120
            C( J ) = ONE / MIN( MAX( C( J ), SMLNUM ), BIGNUM )
         } // 120

         // Compute COLCND = min(C(J)) / max(C(J)).

         COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      }

      RETURN

      // End of CGEEQUB

      }
