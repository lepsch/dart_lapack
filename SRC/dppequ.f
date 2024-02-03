      SUBROUTINE DPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      double             AMAX, SCOND;
      // ..
      // .. Array Arguments ..
      double             AP( * ), S( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, JJ;
      double             SMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      }
      if ( INFO != 0 ) {
         xerbla('DPPEQU', -INFO );
         RETURN
      }

      // Quick return if possible

      if ( N == 0 ) {
         SCOND = ONE
         AMAX = ZERO
         RETURN
      }

      // Initialize SMIN and AMAX.

      S( 1 ) = AP( 1 )
      SMIN = S( 1 )
      AMAX = S( 1 )

      if ( UPPER ) {

         // UPLO = 'U':  Upper triangle of A is stored.
         // Find the minimum and maximum diagonal elements.

         JJ = 1
         for (I = 2; I <= N; I++) { // 10
            JJ = JJ + I
            S( I ) = AP( JJ )
            SMIN = MIN( SMIN, S( I ) )
            AMAX = MAX( AMAX, S( I ) )
         } // 10

      } else {

         // UPLO = 'L':  Lower triangle of A is stored.
         // Find the minimum and maximum diagonal elements.

         JJ = 1
         for (I = 2; I <= N; I++) { // 20
            JJ = JJ + N - I + 2
            S( I ) = AP( JJ )
            SMIN = MIN( SMIN, S( I ) )
            AMAX = MAX( AMAX, S( I ) )
         } // 20
      }

      if ( SMIN.LE.ZERO ) {

         // Find the first non-positive diagonal element and return.

         for (I = 1; I <= N; I++) { // 30
            if ( S( I ).LE.ZERO ) {
               INFO = I
               RETURN
            }
         } // 30
      } else {

         // Set the scale factors to the reciprocals
         // of the diagonal elements.

         for (I = 1; I <= N; I++) { // 40
            S( I ) = ONE / SQRT( S( I ) )
         } // 40

         // Compute SCOND = min(S(I)) / max(S(I))

         SCOND = SQRT( SMIN ) / SQRT( AMAX )
      }
      RETURN

      // End of DPPEQU

      }
