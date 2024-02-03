      REAL             FUNCTION SLARND( IDIST, ISEED )

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
      REAL               ONE, TWO
      const              ONE = 1.0E+0, TWO = 2.0E+0 ;
      REAL               TWOPI
      const      TWOPI = 6.28318530717958647692528676655900576839E+0 ;
      // ..
      // .. Local Scalars ..
      REAL               T1, T2
      // ..
      // .. External Functions ..
      REAL               SLARAN
      // EXTERNAL SLARAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC COS, LOG, SQRT
      // ..
      // .. Executable Statements ..

      // Generate a real random number from a uniform (0,1) distribution

      T1 = SLARAN( ISEED )

      IF( IDIST.EQ.1 ) THEN

         // uniform (0,1)

         SLARND = T1
      ELSE IF( IDIST.EQ.2 ) THEN

         // uniform (-1,1)

         SLARND = TWO*T1 - ONE
      ELSE IF( IDIST.EQ.3 ) THEN

         // normal (0,1)

         T2 = SLARAN( ISEED )
         SLARND = SQRT( -TWO*LOG( T1 ) )*COS( TWOPI*T2 )
      END IF
      RETURN

      // End of SLARND

      }
