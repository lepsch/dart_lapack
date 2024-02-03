      SUBROUTINE DSTECT( N, A, B, SHIFT, NUM )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                N, NUM;
      double             SHIFT;
      // ..
      // .. Array Arguments ..
      double             A( * ), B( * );
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE, THREE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, THREE = 3.0D0 )
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
*
      // Get machine constants
*
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = DLAMCH( 'Overflow' )
*
      // Find largest entry
*
      MX = ABS( A( 1 ) )
      DO 10 I = 1, N - 1
         MX = MAX( MX, ABS( A( I+1 ) ), ABS( B( I ) ) )
   10 CONTINUE
*
      // Handle easy cases, including zero matrix
*
      IF( SHIFT.GE.THREE*MX ) THEN
         NUM = N
         RETURN
      END IF
      IF( SHIFT.LT.-THREE*MX ) THEN
         NUM = 0
         RETURN
      END IF
*
      // Compute scale factors as in Kahan's report
      // At this point, MX .NE. 0 so we can divide by it
*
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
*
      // Begin counting
*
      NUM = 0
      SSHIFT = ( SHIFT*M1 )*M2
      U = ( A( 1 )*M1 )*M2 - SSHIFT
      IF( U.LE.SUN ) THEN
         IF( U.LE.ZERO ) THEN
            NUM = NUM + 1
            IF( U.GT.-SUN ) U = -SUN
         ELSE
            U = SUN
         END IF
      END IF
      DO 20 I = 2, N
         TMP = ( B( I-1 )*M1 )*M2
         U = ( ( A( I )*M1 )*M2-TMP*( TMP / U ) ) - SSHIFT
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
*
      // End of DSTECT
*
      END
