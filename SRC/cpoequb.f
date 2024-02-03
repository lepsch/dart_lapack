      SUBROUTINE CPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, N;
      REAL               AMAX, SCOND;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * );
      REAL               S( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               SMIN, BASE, TMP;
      // ..
      // .. External Functions ..
      REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT, LOG, INT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      // Positive definite only performs 1 pass of equilibration.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -3;
      }
      if ( INFO != 0 ) {
         xerbla('CPOEQUB', -INFO );
         return;
      }

      // Quick return if possible.

      if ( N == 0 ) {
         SCOND = ONE;
         AMAX = ZERO;
         return;
      }

      BASE = SLAMCH( 'B' );
      TMP = -0.5 / LOG ( BASE );

      // Find the minimum and maximum diagonal elements.

      S( 1 ) = REAL( A( 1, 1 ) );
      SMIN = S( 1 );
      AMAX = S( 1 );
      for (I = 2; I <= N; I++) { // 10
         S( I ) = REAL( A( I, I ) );
         SMIN = MIN( SMIN, S( I ) );
         AMAX = MAX( AMAX, S( I ) );
      } // 10

      if ( SMIN <= ZERO ) {

         // Find the first non-positive diagonal element and return.

         for (I = 1; I <= N; I++) { // 20
            if ( S( I ) <= ZERO ) {
               INFO = I;
               return;
            }
         } // 20
      } else {

         // Set the scale factors to the reciprocals
         // of the diagonal elements.

         for (I = 1; I <= N; I++) { // 30
            S( I ) = BASE ** INT( TMP * LOG( S( I ) ) );
         } // 30

         // Compute SCOND = min(S(I)) / max(S(I)).

         SCOND = SQRT( SMIN ) / SQRT( AMAX );
      }

      return;

      // End of CPOEQUB

      }
