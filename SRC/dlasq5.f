      SUBROUTINE DLASQ5( I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, IEEE, EPS )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      bool               IEEE;
      int                I0, N0, PP;
      double             DMIN, DMIN1, DMIN2, DN, DNM1, DNM2, TAU, SIGMA, EPS;
*     ..
*     .. Array Arguments ..
      double             Z( * );
*     ..
*
*  =====================================================================
*
*     .. Parameter ..
      double             ZERO, HALF;
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5 )
*     ..
*     .. Local Scalars ..
      int                J4, J4P2;
      double             D, EMIN, TEMP, DTHRESH;
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( ( N0-I0-1 ).LE.0 ) RETURN
*
      DTHRESH = EPS*(SIGMA+TAU)
      IF( TAU.LT.DTHRESH*HALF ) TAU = ZERO
      IF( TAU.NE.ZERO ) THEN
      J4 = 4*I0 + PP - 3
      EMIN = Z( J4+4 )
      D = Z( J4 ) - TAU
      DMIN = D
      DMIN1 = -Z( J4 )
*
      IF( IEEE ) THEN
*
*        Code for IEEE arithmetic.
*
         IF( PP.EQ.0 ) THEN
            DO 10 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-2 ) = D + Z( J4-1 )
               TEMP = Z( J4+1 ) / Z( J4-2 )
               D = D*TEMP - TAU
               DMIN = MIN( DMIN, D )
               Z( J4 ) = Z( J4-1 )*TEMP
               EMIN = MIN( Z( J4 ), EMIN )
   10       CONTINUE
         ELSE
            DO 20 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-3 ) = D + Z( J4 )
               TEMP = Z( J4+2 ) / Z( J4-3 )
               D = D*TEMP - TAU
               DMIN = MIN( DMIN, D )
               Z( J4-1 ) = Z( J4 )*TEMP
               EMIN = MIN( Z( J4-1 ), EMIN )
   20       CONTINUE
         END IF
*
*        Unroll last two steps.
*
         DNM2 = D
         DMIN2 = DMIN
         J4 = 4*( N0-2 ) - PP
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM2 + Z( J4P2 )
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
         DMIN = MIN( DMIN, DNM1 )
*
         DMIN1 = DMIN
         J4 = J4 + 4
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM1 + Z( J4P2 )
         Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
         DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
         DMIN = MIN( DMIN, DN )
*
      ELSE
*
*        Code for non IEEE arithmetic.
*
         IF( PP.EQ.0 ) THEN
            DO 30 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-2 ) = D + Z( J4-1 )
               IF( D.LT.ZERO ) THEN
                  RETURN
               ELSE
                  Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
                  D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU
               END IF
               DMIN = MIN( DMIN, D )
               EMIN = MIN( EMIN, Z( J4 ) )
   30       CONTINUE
         ELSE
            DO 40 J4 = 4*I0, 4*( N0-3 ), 4
               Z( J4-3 ) = D + Z( J4 )
               IF( D.LT.ZERO ) THEN
                  RETURN
               ELSE
                  Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
                  D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU
               END IF
               DMIN = MIN( DMIN, D )
               EMIN = MIN( EMIN, Z( J4-1 ) )
   40       CONTINUE
         END IF
*
*        Unroll last two steps.
*
         DNM2 = D
         DMIN2 = DMIN
         J4 = 4*( N0-2 ) - PP
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM2 + Z( J4P2 )
         IF( DNM2.LT.ZERO ) THEN
            RETURN
         ELSE
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
         END IF
         DMIN = MIN( DMIN, DNM1 )
*
         DMIN1 = DMIN
         J4 = J4 + 4
         J4P2 = J4 + 2*PP - 1
         Z( J4-2 ) = DNM1 + Z( J4P2 )
         IF( DNM1.LT.ZERO ) THEN
            RETURN
         ELSE
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
         END IF
         DMIN = MIN( DMIN, DN )
*
      END IF
      ELSE
*     This is the version that sets d's to zero if they are small enough
         J4 = 4*I0 + PP - 3
         EMIN = Z( J4+4 )
         D = Z( J4 ) - TAU
         DMIN = D
         DMIN1 = -Z( J4 )
         IF( IEEE ) THEN
*
*     Code for IEEE arithmetic.
*
            IF( PP.EQ.0 ) THEN
               DO 50 J4 = 4*I0, 4*( N0-3 ), 4
                  Z( J4-2 ) = D + Z( J4-1 )
                  TEMP = Z( J4+1 ) / Z( J4-2 )
                  D = D*TEMP - TAU
                  IF( D.LT.DTHRESH ) D = ZERO
                  DMIN = MIN( DMIN, D )
                  Z( J4 ) = Z( J4-1 )*TEMP
                  EMIN = MIN( Z( J4 ), EMIN )
 50            CONTINUE
            ELSE
               DO 60 J4 = 4*I0, 4*( N0-3 ), 4
                  Z( J4-3 ) = D + Z( J4 )
                  TEMP = Z( J4+2 ) / Z( J4-3 )
                  D = D*TEMP - TAU
                  IF( D.LT.DTHRESH ) D = ZERO
                  DMIN = MIN( DMIN, D )
                  Z( J4-1 ) = Z( J4 )*TEMP
                  EMIN = MIN( Z( J4-1 ), EMIN )
 60            CONTINUE
            END IF
*
*     Unroll last two steps.
*
            DNM2 = D
            DMIN2 = DMIN
            J4 = 4*( N0-2 ) - PP
            J4P2 = J4 + 2*PP - 1
            Z( J4-2 ) = DNM2 + Z( J4P2 )
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
            DMIN = MIN( DMIN, DNM1 )
*
            DMIN1 = DMIN
            J4 = J4 + 4
            J4P2 = J4 + 2*PP - 1
            Z( J4-2 ) = DNM1 + Z( J4P2 )
            Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
            DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
            DMIN = MIN( DMIN, DN )
*
         ELSE
*
*     Code for non IEEE arithmetic.
*
            IF( PP.EQ.0 ) THEN
               DO 70 J4 = 4*I0, 4*( N0-3 ), 4
                  Z( J4-2 ) = D + Z( J4-1 )
                  IF( D.LT.ZERO ) THEN
                     RETURN
                  ELSE
                     Z( J4 ) = Z( J4+1 )*( Z( J4-1 ) / Z( J4-2 ) )
                     D = Z( J4+1 )*( D / Z( J4-2 ) ) - TAU
                  END IF
                  IF( D.LT.DTHRESH) D = ZERO
                  DMIN = MIN( DMIN, D )
                  EMIN = MIN( EMIN, Z( J4 ) )
 70            CONTINUE
            ELSE
               DO 80 J4 = 4*I0, 4*( N0-3 ), 4
                  Z( J4-3 ) = D + Z( J4 )
                  IF( D.LT.ZERO ) THEN
                     RETURN
                  ELSE
                     Z( J4-1 ) = Z( J4+2 )*( Z( J4 ) / Z( J4-3 ) )
                     D = Z( J4+2 )*( D / Z( J4-3 ) ) - TAU
                  END IF
                  IF( D.LT.DTHRESH) D = ZERO
                  DMIN = MIN( DMIN, D )
                  EMIN = MIN( EMIN, Z( J4-1 ) )
 80            CONTINUE
            END IF
*
*     Unroll last two steps.
*
            DNM2 = D
            DMIN2 = DMIN
            J4 = 4*( N0-2 ) - PP
            J4P2 = J4 + 2*PP - 1
            Z( J4-2 ) = DNM2 + Z( J4P2 )
            IF( DNM2.LT.ZERO ) THEN
               RETURN
            ELSE
               Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
               DNM1 = Z( J4P2+2 )*( DNM2 / Z( J4-2 ) ) - TAU
            END IF
            DMIN = MIN( DMIN, DNM1 )
*
            DMIN1 = DMIN
            J4 = J4 + 4
            J4P2 = J4 + 2*PP - 1
            Z( J4-2 ) = DNM1 + Z( J4P2 )
            IF( DNM1.LT.ZERO ) THEN
               RETURN
            ELSE
               Z( J4 ) = Z( J4P2+2 )*( Z( J4P2 ) / Z( J4-2 ) )
               DN = Z( J4P2+2 )*( DNM1 / Z( J4-2 ) ) - TAU
            END IF
            DMIN = MIN( DMIN, DN )
*
         END IF
      END IF
*
      Z( J4+2 ) = DN
      Z( 4*N0-PP ) = EMIN
      RETURN
*
*     End of DLASQ5
*
      END
