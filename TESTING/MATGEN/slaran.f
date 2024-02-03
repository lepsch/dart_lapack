      REAL FUNCTION SLARAN( ISEED )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Array Arguments ..
      int                ISEED( 4 );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                M1, M2, M3, M4;
      PARAMETER          ( M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      int                IPW2;
      REAL               R
      PARAMETER          ( IPW2 = 4096, R = ONE / IPW2 )
      // ..
      // .. Local Scalars ..
      int                IT1, IT2, IT3, IT4;
      REAL               RNDOUT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD, REAL
      // ..
      // .. Executable Statements ..
  10  CONTINUE

      // multiply the seed by the multiplier modulo 2**48

      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 + ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )

      // return updated seed

      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4

      // convert 48-bit integer to a real number in the interval (0,1)

      RNDOUT = R*( REAL( IT1 )+R*( REAL( IT2 )+R*( REAL( IT3 )+R* ( REAL( IT4 ) ) ) ) )

      IF (RNDOUT.EQ.1.0) THEN
         // If a real number has n bits of precision, and the first
         // n bits of the 48-bit integer above happen to be all 1 (which
         // will occur about once every 2**n calls), then SLARAN will
         // be rounded to exactly 1.0. In IEEE single precision arithmetic,
        t // his will happen relatively often since n = 24.
         // Since SLARAN is not supposed to return exactly 0.0 or 1.0
         // (and some callers of SLARAN, such as CLARND, depend on that),
        t // he statistically correct thing to do in this situation is
         // simply to iterate again.
         // N.B. the case SLARAN = 0.0 should not be possible.

         GOTO 10
      END IF

      SLARAN = RNDOUT
      RETURN

      // End of SLARAN

      }
