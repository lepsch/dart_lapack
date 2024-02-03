      COMPLEX*16   FUNCTION ZLARND( IDIST, ISEED )

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
      double             ZERO, ONE, TWO;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
      double             TWOPI;
      PARAMETER  ( TWOPI = 6.28318530717958647692528676655900576839D+0 )
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

      T1 = DLARAN( ISEED )
      T2 = DLARAN( ISEED )

      IF( IDIST.EQ.1 ) THEN

         // real and imaginary parts each uniform (0,1)

         ZLARND = DCMPLX( T1, T2 )
      ELSE IF( IDIST.EQ.2 ) THEN

         // real and imaginary parts each uniform (-1,1)

         ZLARND = DCMPLX( TWO*T1-ONE, TWO*T2-ONE )
      ELSE IF( IDIST.EQ.3 ) THEN

         // real and imaginary parts each normal (0,1)

         ZLARND = SQRT( -TWO*LOG( T1 ) )*EXP( DCMPLX( ZERO, TWOPI*T2 ) )
      ELSE IF( IDIST.EQ.4 ) THEN

         // uniform distribution on the unit disc abs(z) <= 1

         ZLARND = SQRT( T1 )*EXP( DCMPLX( ZERO, TWOPI*T2 ) )
      ELSE IF( IDIST.EQ.5 ) THEN

         // uniform distribution on the unit circle abs(z) = 1

         ZLARND = EXP( DCMPLX( ZERO, TWOPI*T2 ) )
      END IF
      RETURN

      // End of ZLARND

      END
