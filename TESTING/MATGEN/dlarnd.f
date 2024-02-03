      double           FUNCTION DLARND( IDIST, ISEED );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IDIST;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, TWO;
      const              ONE = 1.0, TWO = 2.0 ;
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
      // INTRINSIC COS, LOG, SQRT
      // ..
      // .. Executable Statements ..

      // Generate a real random number from a uniform (0,1) distribution

      T1 = DLARAN( ISEED );

      if ( IDIST == 1 ) {

         // uniform (0,1)

         DLARND = T1;
      } else if ( IDIST == 2 ) {

         // uniform (-1,1)

         DLARND = TWO*T1 - ONE;
      } else if ( IDIST == 3 ) {

         // normal (0,1)

         T2 = DLARAN( ISEED );
         DLARND = SQRT( -TWO*LOG( T1 ) )*COS( TWOPI*T2 );
      }
      RETURN;

      // End of DLARND

      }
