      double dlaran(ISEED ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Array Arguments ..
      int                ISEED( 4 );
      // ..

// =====================================================================

      // .. Parameters ..
      int                M1, M2, M3, M4;
      const              M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 ;
      double             ONE;
      const              ONE = 1.0 ;
      int                IPW2;
      double             R;
      const              IPW2 = 4096, R = ONE / IPW2 ;
      // ..
      // .. Local Scalars ..
      int                IT1, IT2, IT3, IT4;
      double             RNDOUT;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MOD
      // ..
      // .. Executable Statements ..
      } // 10

      // multiply the seed by the multiplier modulo 2**48

      IT4 = ISEED( 4 )*M4;
      IT3 = IT4 / IPW2;
      IT4 = IT4 - IPW2*IT3;
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3;
      IT2 = IT3 / IPW2;
      IT3 = IT3 - IPW2*IT2;
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2;
      IT1 = IT2 / IPW2;
      IT2 = IT2 - IPW2*IT1;
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 + ISEED( 4 )*M1;
      IT1 = (IT1 % IPW2);

      // return updated seed

      ISEED( 1 ) = IT1;
      ISEED( 2 ) = IT2;
      ISEED( 3 ) = IT3;
      ISEED( 4 ) = IT4;

      // convert 48-bit integer to a real number in the interval (0,1)

      RNDOUT = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )+R* ( DBLE( IT4 ) ) ) ) );

      if (RNDOUT == 1.0) {
         // If a real number has n bits of precision, and the first
         // n bits of the 48-bit integer above happen to be all 1 (which
         // will occur about once every 2**n calls), then DLARAN will
         // be rounded to exactly 1.0.
         // Since DLARAN is not supposed to return exactly 0.0 or 1.0
         // (and some callers of DLARAN, such as CLARND, depend on that),
         // the statistically correct thing to do in this situation is
         // simply to iterate again.
         // N.B. the case DLARAN = 0.0 should not be possible.

         GOTO 10;
      }

      DLARAN = RNDOUT;
      return;
      }
