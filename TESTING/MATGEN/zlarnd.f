      Complex   FUNCTION ZLARND( IDIST, ISEED );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IDIST;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      double             TWOPI;
      const      TWOPI = 6.28318530717958647692528676655900576839 ;
      // ..
      // .. Local Scalars ..
      double             T1, T2;
      // ..
      // .. External Functions ..
      double             DLARAN;
      // EXTERNAL DLARAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, EXP, LOG, SQRT
      // ..
      // .. Executable Statements ..

      // Generate a pair of real random numbers from a uniform (0,1)
      // distribution

      T1 = DLARAN( ISEED );
      T2 = DLARAN( ISEED );

      if ( IDIST == 1 ) {

         // real and imaginary parts each uniform (0,1)

         ZLARND = DCMPLX( T1, T2 );
      } else if ( IDIST == 2 ) {

         // real and imaginary parts each uniform (-1,1)

         ZLARND = DCMPLX( TWO*T1-ONE, TWO*T2-ONE );
      } else if ( IDIST == 3 ) {

         // real and imaginary parts each normal (0,1)

         ZLARND = sqrt( -TWO*LOG( T1 ) )*EXP( DCMPLX( ZERO, TWOPI*T2 ) );
      } else if ( IDIST == 4 ) {

         // uniform distribution on the unit disc abs(z) <= 1

         ZLARND = sqrt( T1 )*EXP( DCMPLX( ZERO, TWOPI*T2 ) );
      } else if ( IDIST == 5 ) {

         // uniform distribution on the unit circle abs(z) = 1

         ZLARND = EXP( DCMPLX( ZERO, TWOPI*T2 ) );
      }
      return;
      }
