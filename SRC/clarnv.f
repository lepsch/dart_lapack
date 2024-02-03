      SUBROUTINE CLARNV( IDIST, ISEED, N, X )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IDIST, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      COMPLEX            X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 ;
      int                LV;
      const              LV = 128 ;
      REAL               TWOPI
      const      TWOPI = 6.28318530717958647692528676655900576839E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, IL, IV;
      // ..
      // .. Local Arrays ..
      REAL               U( LV )
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, EXP, LOG, MIN, SQRT
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARUV
      // ..
      // .. Executable Statements ..

      DO 60 IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )

         // Call SLARUV to generate 2*IL real numbers from a uniform (0,1)
         // distribution (2*IL <= LV)

         slaruv(ISEED, 2*IL, U );

         if ( IDIST == 1 ) {

            // Copy generated numbers

            for (I = 1; I <= IL; I++) { // 10
               X( IV+I-1 ) = CMPLX( U( 2*I-1 ), U( 2*I ) )
            } // 10
         } else if ( IDIST == 2 ) {

            // Convert generated numbers to uniform (-1,1) distribution

            for (I = 1; I <= IL; I++) { // 20
               X( IV+I-1 ) = CMPLX( TWO*U( 2*I-1 )-ONE, TWO*U( 2*I )-ONE )
            } // 20
         } else if ( IDIST == 3 ) {

            // Convert generated numbers to normal (0,1) distribution

            for (I = 1; I <= IL; I++) { // 30
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )* EXP( CMPLX( ZERO, TWOPI*U( 2*I ) ) )
            } // 30
         } else if ( IDIST == 4 ) {

            // Convert generated numbers to complex numbers uniformly
            // distributed on the unit disk

            for (I = 1; I <= IL; I++) { // 40
               X( IV+I-1 ) = SQRT( U( 2*I-1 ) )* EXP( CMPLX( ZERO, TWOPI*U( 2*I ) ) )
            } // 40
         } else if ( IDIST == 5 ) {

            // Convert generated numbers to complex numbers uniformly
            // distributed on the unit circle

            for (I = 1; I <= IL; I++) { // 50
               X( IV+I-1 ) = EXP( CMPLX( ZERO, TWOPI*U( 2*I ) ) )
            } // 50
         }
      } // 60
      RETURN

      // End of CLARNV

      }
