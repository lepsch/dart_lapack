      SUBROUTINE SPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, N;
      REAL               AMAX, SCOND
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), S( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               SMIN, BASE, TMP
      // ..
      // .. External Functions ..
      REAL               SLAMCH
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

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -1
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -3
      }
      if ( INFO != 0 ) {
         xerbla('SPOEQUB', -INFO );
         RETURN
      }

      // Quick return if possible.

      if ( N == 0 ) {
         SCOND = ONE
         AMAX = ZERO
         RETURN
      }

      BASE = SLAMCH( 'B' )
      TMP = -0.5 / LOG ( BASE )

      // Find the minimum and maximum diagonal elements.

      S( 1 ) = A( 1, 1 )
      SMIN = S( 1 )
      AMAX = S( 1 )
      for (I = 2; I <= N; I++) { // 10
         S( I ) = A( I, I )
         SMIN = MIN( SMIN, S( I ) )
         AMAX = MAX( AMAX, S( I ) )
      } // 10

      if ( SMIN.LE.ZERO ) {

         // Find the first non-positive diagonal element and return.

         for (I = 1; I <= N; I++) { // 20
            if ( S( I ).LE.ZERO ) {
               INFO = I
               RETURN
            }
         } // 20
      } else {

         // Set the scale factors to the reciprocals
         // of the diagonal elements.

         for (I = 1; I <= N; I++) { // 30
            S( I ) = BASE ** INT( TMP * LOG( S( I ) ) )
         } // 30

         // Compute SCOND = min(S(I)) / max(S(I)).

         SCOND = SQRT( SMIN ) / SQRT( AMAX )
      }

      RETURN

      // End of SPOEQUB

      }
