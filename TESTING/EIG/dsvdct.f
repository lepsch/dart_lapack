      SUBROUTINE DSVDCT( N, S, E, SHIFT, NUM )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                N, NUM;
      double             SHIFT;
      // ..
      // .. Array Arguments ..
      double             E( * ), S( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      PARAMETER          ( ONE = 1.0D0 )
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D0 )
      // ..
      // .. Local Scalars ..
      int                I;
      double             M1, M2, MX, OVFL, SOV, SSHIFT, SSUN, SUN, TMP, TOM, U, UNFL;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Get machine constants

      UNFL = 2*DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL

      // Find largest entry

      MX = ABS( S( 1 ) )
      DO 10 I = 1, N - 1
         MX = MAX( MX, ABS( S( I+1 ) ), ABS( E( I ) ) )
   10 CONTINUE

      IF( MX.EQ.ZERO ) THEN
         IF( SHIFT.LT.ZERO ) THEN
            NUM = 0
         ELSE
            NUM = 2*N
         END IF
         RETURN
      END IF

      // Compute scale factors as in Kahan's report

      SUN = SQRT( UNFL )
      SSUN = SQRT( SUN )
      SOV = SQRT( OVFL )
      TOM = SSUN*SOV
      IF( MX.LE.ONE ) THEN
         M1 = ONE / MX
         M2 = TOM
      ELSE
         M1 = ONE
         M2 = TOM / MX
      END IF

      // Begin counting

      U = ONE
      NUM = 0
      SSHIFT = ( SHIFT*M1 )*M2
      U = -SSHIFT
      IF( U.LE.SUN ) THEN
         IF( U.LE.ZERO ) THEN
            NUM = NUM + 1
            IF( U.GT.-SUN ) U = -SUN
         ELSE
            U = SUN
         END IF
      END IF
      TMP = ( S( 1 )*M1 )*M2
      U = -TMP*( TMP / U ) - SSHIFT
      IF( U.LE.SUN ) THEN
         IF( U.LE.ZERO ) THEN
            NUM = NUM + 1
            IF( U.GT.-SUN ) U = -SUN
         ELSE
            U = SUN
         END IF
      END IF
      DO 20 I = 1, N - 1
         TMP = ( E( I )*M1 )*M2
         U = -TMP*( TMP / U ) - SSHIFT
         IF( U.LE.SUN ) THEN
            IF( U.LE.ZERO ) THEN
               NUM = NUM + 1
               IF( U.GT.-SUN ) U = -SUN
            ELSE
               U = SUN
            END IF
         END IF
         TMP = ( S( I+1 )*M1 )*M2
         U = -TMP*( TMP / U ) - SSHIFT
         IF( U.LE.SUN ) THEN
            IF( U.LE.ZERO ) THEN
               NUM = NUM + 1
               IF( U.GT.-SUN ) U = -SUN
            ELSE
               U = SUN
            END IF
         END IF
   20 CONTINUE
      RETURN

      // End of DSVDCT

      END
