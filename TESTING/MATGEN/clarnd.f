      COMPLEX FUNCTION CLARND( IDIST, ISEED )

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
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 ;
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
      // INTRINSIC CMPLX, EXP, LOG, SQRT
      // ..
      // .. Executable Statements ..

      // Generate a pair of real random numbers from a uniform (0,1)
      // distribution

      T1 = SLARAN( ISEED )
      T2 = SLARAN( ISEED )

      IF( IDIST.EQ.1 ) THEN

         // real and imaginary parts each uniform (0,1)

         CLARND = CMPLX( T1, T2 )
      ELSE IF( IDIST.EQ.2 ) THEN

         // real and imaginary parts each uniform (-1,1)

         CLARND = CMPLX( TWO*T1-ONE, TWO*T2-ONE )
      ELSE IF( IDIST.EQ.3 ) THEN

         // real and imaginary parts each normal (0,1)

         CLARND = SQRT( -TWO*LOG( T1 ) )*EXP( CMPLX( ZERO, TWOPI*T2 ) )
      ELSE IF( IDIST.EQ.4 ) THEN

         // uniform distribution on the unit disc abs(z) <= 1

         CLARND = SQRT( T1 )*EXP( CMPLX( ZERO, TWOPI*T2 ) )
      ELSE IF( IDIST.EQ.5 ) THEN

         // uniform distribution on the unit circle abs(z) = 1

         CLARND = EXP( CMPLX( ZERO, TWOPI*T2 ) )
      END IF
      RETURN

      // End of CLARND

      }
