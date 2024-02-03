      SUBROUTINE DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN;
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, HALF, ONE, TWO;
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0, TWO = 2.0D0 )
      double             MULTPL;
      PARAMETER          ( MULTPL = 4.0D+0 )
      // ..
      // .. Local Scalars ..
      double             AA, BB, BCMAX, BCMIS, CC, CS1, DD, EPS, P, SAB, SAC, SCALE, SIGMA, SN1, TAU, TEMP, Z, SAFMIN, SAFMN2, SAFMX2;
      int                COUNT;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLAPY2;
      // EXTERNAL DLAMCH, DLAPY2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SIGN, SQRT
      // ..
      // .. Executable Statements ..

      SAFMIN = DLAMCH( 'S' )
      EPS = DLAMCH( 'P' )
      SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / LOG( DLAMCH( 'B' ) ) / TWO )
      SAFMX2 = ONE / SAFMN2
      IF( C.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO

      ELSE IF( B.EQ.ZERO ) THEN

         // Swap rows and columns

         CS = ZERO
         SN = ONE
         TEMP = D
         D = A
         A = TEMP
         B = -C
         C = ZERO

      ELSE IF( ( A-D ).EQ.ZERO .AND. SIGN( ONE, B ).NE.SIGN( ONE, C ) ) THEN
         CS = ONE
         SN = ZERO

      ELSE

         TEMP = A - D
         P = HALF*TEMP
         BCMAX = MAX( ABS( B ), ABS( C ) )
         BCMIS = MIN( ABS( B ), ABS( C ) )*SIGN( ONE, B )*SIGN( ONE, C )
         SCALE = MAX( ABS( P ), BCMAX )
         Z = ( P / SCALE )*P + ( BCMAX / SCALE )*BCMIS

         // If Z is of the order of the machine accuracy, postpone the
         // decision on the nature of eigenvalues

         IF( Z.GE.MULTPL*EPS ) THEN

            // Real eigenvalues. Compute A and D.

            Z = P + SIGN( SQRT( SCALE )*SQRT( Z ), P )
            A = D + Z
            D = D - ( BCMAX / Z )*BCMIS

            // Compute B and the rotation matrix

            TAU = DLAPY2( C, Z )
            CS = Z / TAU
            SN = C / TAU
            B = B - C
            C = ZERO

         ELSE

            // Complex eigenvalues, or real (almost) equal eigenvalues.
            // Make diagonal elements equal.

            COUNT = 0
            SIGMA = B + C
   10       CONTINUE
            COUNT = COUNT + 1
            SCALE = MAX( ABS(TEMP), ABS(SIGMA) )
            IF( SCALE.GE.SAFMX2 ) THEN
               SIGMA = SIGMA * SAFMN2
               TEMP = TEMP * SAFMN2
               IF (COUNT .LE. 20) GOTO 10
            END IF
            IF( SCALE.LE.SAFMN2 ) THEN
               SIGMA = SIGMA * SAFMX2
               TEMP = TEMP * SAFMX2
               IF (COUNT .LE. 20) GOTO 10
            END IF
            P = HALF*TEMP
            TAU = DLAPY2( SIGMA, TEMP )
            CS = SQRT( HALF*( ONE+ABS( SIGMA ) / TAU ) )
            SN = -( P / ( TAU*CS ) )*SIGN( ONE, SIGMA )

            // Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
                    // [ CC  DD ]   [ C  D ] [ SN  CS ]

            AA = A*CS + B*SN
            BB = -A*SN + B*CS
            CC = C*CS + D*SN
            DD = -C*SN + D*CS

            // Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
                    // [ C  D ]   [-SN  CS ] [ CC  DD ]

            A = AA*CS + CC*SN
            B = BB*CS + DD*SN
            C = -AA*SN + CC*CS
            D = -BB*SN + DD*CS

            TEMP = HALF*( A+D )
            A = TEMP
            D = TEMP

            IF( C.NE.ZERO ) THEN
               IF( B.NE.ZERO ) THEN
                  IF( SIGN( ONE, B ).EQ.SIGN( ONE, C ) ) THEN

                     // Real eigenvalues: reduce to upper triangular form

                     SAB = SQRT( ABS( B ) )
                     SAC = SQRT( ABS( C ) )
                     P = SIGN( SAB*SAC, C )
                     TAU = ONE / SQRT( ABS( B+C ) )
                     A = TEMP + P
                     D = TEMP - P
                     B = B - C
                     C = ZERO
                     CS1 = SAB*TAU
                     SN1 = SAC*TAU
                     TEMP = CS*CS1 - SN*SN1
                     SN = CS*SN1 + SN*CS1
                     CS = TEMP
                  END IF
               ELSE
                  B = -C
                  C = ZERO
                  TEMP = CS
                  CS = -SN
                  SN = TEMP
               END IF
            END IF
         END IF

      END IF

      // Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).

      RT1R = A
      RT2R = D
      IF( C.EQ.ZERO ) THEN
         RT1I = ZERO
         RT2I = ZERO
      ELSE
         RT1I = SQRT( ABS( B ) )*SQRT( ABS( C ) )
         RT2I = -RT1I
      END IF
      RETURN

      // End of DLANV2

      END
