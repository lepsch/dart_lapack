      void clarnv(IDIST, ISEED, N, X ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IDIST, N;
      int                ISEED( 4 );
      Complex            X( * );
      // ..

      double               ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      int                LV;
      const              LV = 128 ;
      double               TWOPI;
      const      TWOPI = 6.28318530717958647692528676655900576839 ;
      int                I, IL, IV;
      double               U( LV );
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, EXP, LOG, MIN, SQRT
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARUV

      for (IV = 1; LV / 2 < 0 ? IV >= N : IV <= N; IV += LV / 2) { // 60
         IL = min( LV / 2, N-IV+1 );

         // Call SLARUV to generate 2*IL real numbers from a uniform (0,1)
         // distribution (2*IL <= LV)

         slaruv(ISEED, 2*IL, U );

         if ( IDIST == 1 ) {

            // Copy generated numbers

            for (I = 1; I <= IL; I++) { // 10
               X[IV+I-1] = CMPLX( U( 2*I-1 ), U( 2*I ) );
            } // 10
         } else if ( IDIST == 2 ) {

            // Convert generated numbers to uniform (-1,1) distribution

            for (I = 1; I <= IL; I++) { // 20
               X[IV+I-1] = CMPLX( TWO*U( 2*I-1 )-ONE, TWO*U( 2*I )-ONE );
            } // 20
         } else if ( IDIST == 3 ) {

            // Convert generated numbers to normal (0,1) distribution

            for (I = 1; I <= IL; I++) { // 30
               X[IV+I-1] = sqrt( -TWO*LOG( U( 2*I-1 ) ) )* EXP( CMPLX( ZERO, TWOPI*U( 2*I ) ) );
            } // 30
         } else if ( IDIST == 4 ) {

            // Convert generated numbers to complex numbers uniformly
            // distributed on the unit disk

            for (I = 1; I <= IL; I++) { // 40
               X[IV+I-1] = sqrt( U( 2*I-1 ) )* EXP( CMPLX( ZERO, TWOPI*U( 2*I ) ) );
            } // 40
         } else if ( IDIST == 5 ) {

            // Convert generated numbers to complex numbers uniformly
            // distributed on the unit circle

            for (I = 1; I <= IL; I++) { // 50
               X[IV+I-1] = EXP( CMPLX( ZERO, TWOPI*U( 2*I ) ) );
            } // 50
         }
      } // 60
      return;
      }
