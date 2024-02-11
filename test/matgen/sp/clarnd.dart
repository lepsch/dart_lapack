      Complex clarnd(final int IDIST, final int ISEED,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IDIST;
      int                ISEED( 4 );
      // ..

      double               ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      double               TWOPI;
      const      TWOPI = 6.28318530717958647692528676655900576839 ;
      double               T1, T2;
      // ..
      // .. External Functions ..
      //- REAL               SLARAN;
      // EXTERNAL SLARAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, EXP, LOG, SQRT

      // Generate a pair of real random numbers from a uniform (0,1)
      // distribution

      T1 = SLARAN( ISEED );
      T2 = SLARAN( ISEED );

      if ( IDIST == 1 ) {

         // real and imaginary parts each uniform (0,1)

         CLARND = CMPLX( T1, T2 );
      } else if ( IDIST == 2 ) {

         // real and imaginary parts each uniform (-1,1)

         CLARND = CMPLX( TWO*T1-ONE, TWO*T2-ONE );
      } else if ( IDIST == 3 ) {

         // real and imaginary parts each normal (0,1)

         CLARND = sqrt( -TWO*LOG( T1 ) )*EXP( CMPLX( ZERO, TWOPI*T2 ) );
      } else if ( IDIST == 4 ) {

         // uniform distribution on the unit disc abs(z) <= 1

         CLARND = sqrt( T1 )*EXP( CMPLX( ZERO, TWOPI*T2 ) );
      } else if ( IDIST == 5 ) {

         // uniform distribution on the unit circle abs(z) = 1

         CLARND = EXP( CMPLX( ZERO, TWOPI*T2 ) );
      }
      }
