      REAL             FUNCTION SLARND( IDIST, ISEED );

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
      REAL               ONE, TWO;
      const              ONE = 1.0, TWO = 2.0 ;
      REAL               TWOPI;
      const      TWOPI = 6.28318530717958647692528676655900576839 ;
      // ..
      // .. Local Scalars ..
      REAL               T1, T2;
      // ..
      // .. External Functions ..
      REAL               SLARAN;
      // EXTERNAL SLARAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC COS, LOG, SQRT
      // ..
      // .. Executable Statements ..

      // Generate a real random number from a uniform (0,1) distribution

      T1 = SLARAN( ISEED );

      if ( IDIST == 1 ) {

         // uniform (0,1)

         SLARND = T1;
      } else if ( IDIST == 2 ) {

         // uniform (-1,1)

         SLARND = TWO*T1 - ONE;
      } else if ( IDIST == 3 ) {

         // normal (0,1)

         T2 = SLARAN( ISEED );
         SLARND = SQRT( -TWO*LOG( T1 ) )*COS( TWOPI*T2 );
      }
      return;

      // End of SLARND

      }
