      void sppequ(final int UPLO, final int N, final int AP, final int S, final int SCOND, final int AMAX, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, N;
      double               AMAX, SCOND;
      double               AP( * ), S( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               UPPER;
      int                I, JJ;
      double               SMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('SPPEQU', -INFO );
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         SCOND = ONE;
         AMAX = ZERO;
         return;
      }

      // Initialize SMIN and AMAX.

      S[1] = AP( 1 );
      SMIN = S( 1 );
      AMAX = S( 1 );

      if ( UPPER ) {

         // UPLO = 'U':  Upper triangle of A is stored.
         // Find the minimum and maximum diagonal elements.

         JJ = 1;
         for (I = 2; I <= N; I++) { // 10
            JJ = JJ + I;
            S[I] = AP( JJ );
            SMIN = min( SMIN, S( I ) );
            AMAX = max( AMAX, S( I ) );
         } // 10

      } else {

         // UPLO = 'L':  Lower triangle of A is stored.
         // Find the minimum and maximum diagonal elements.

         JJ = 1;
         for (I = 2; I <= N; I++) { // 20
            JJ = JJ + N - I + 2;
            S[I] = AP( JJ );
            SMIN = min( SMIN, S( I ) );
            AMAX = max( AMAX, S( I ) );
         } // 20
      }

      if ( SMIN <= ZERO ) {

         // Find the first non-positive diagonal element and return.

         for (I = 1; I <= N; I++) { // 30
            if ( S( I ) <= ZERO ) {
               INFO = I;
               return;
            }
         } // 30
      } else {

         // Set the scale factors to the reciprocals
         // of the diagonal elements.

         for (I = 1; I <= N; I++) { // 40
            S[I] = ONE / sqrt( S( I ) );
         } // 40

         // Compute SCOND = min(S(I)) / max(S(I))

         SCOND = sqrt( SMIN ) / sqrt( AMAX );
      }
      }
