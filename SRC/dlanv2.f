      void dlanv2(A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN;
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, HALF, ONE, TWO;
      const              ZERO = 0.0, HALF = 0.5, ONE = 1.0, TWO = 2.0 ;
      double             MULTPL;
      const              MULTPL = 4.0 ;
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

      SAFMIN = DLAMCH( 'S' );
      EPS = DLAMCH( 'P' );
      SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / LOG( DLAMCH( 'B' ) ) / TWO );
      SAFMX2 = ONE / SAFMN2;
      if ( C == ZERO ) {
         CS = ONE;
         SN = ZERO;

      } else if ( B == ZERO ) {

         // Swap rows and columns

         CS = ZERO;
         SN = ONE;
         TEMP = D;
         D = A;
         A = TEMP;
         B = -C;
         C = ZERO;

      } else if ( ( A-D ) == ZERO && SIGN( ONE, B ) != SIGN( ONE, C ) ) {
         CS = ONE;
         SN = ZERO;

      } else {

         TEMP = A - D;
         P = HALF*TEMP;
         BCMAX = max( ABS( B ), ABS( C ) );
         BCMIS = min( ABS( B ), ABS( C ) )*SIGN( ONE, B )*SIGN( ONE, C );
         SCALE = max( ABS( P ), BCMAX );
         Z = ( P / SCALE )*P + ( BCMAX / SCALE )*BCMIS;

         // If Z is of the order of the machine accuracy, postpone the
         // decision on the nature of eigenvalues

         if ( Z >= MULTPL*EPS ) {

            // Real eigenvalues. Compute A and D.

            Z = P + SIGN( sqrt( SCALE )*sqrt( Z ), P );
            A = D + Z;
            D = D - ( BCMAX / Z )*BCMIS;

            // Compute B and the rotation matrix

            TAU = DLAPY2( C, Z );
            CS = Z / TAU;
            SN = C / TAU;
            B = B - C;
            C = ZERO;

         } else {

            // Complex eigenvalues, or real (almost) equal eigenvalues.
            // Make diagonal elements equal.

            COUNT = 0;
            SIGMA = B + C;
            } // 10
            COUNT = COUNT + 1;
            SCALE = max( ABS(TEMP), ABS(SIGMA) );
            if ( SCALE >= SAFMX2 ) {
               SIGMA = SIGMA * SAFMN2;
               TEMP = TEMP * SAFMN2;
               if (COUNT <= 20) GOTO 10;
            }
            if ( SCALE <= SAFMN2 ) {
               SIGMA = SIGMA * SAFMX2;
               TEMP = TEMP * SAFMX2;
               if (COUNT <= 20) GOTO 10;
            }
            P = HALF*TEMP;
            TAU = DLAPY2( SIGMA, TEMP );
            CS = sqrt( HALF*( ONE+ABS( SIGMA ) / TAU ) );
            SN = -( P / ( TAU*CS ) )*SIGN( ONE, SIGMA );

            // Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
                    // [ CC  DD ]   [ C  D ] [ SN  CS ]

            AA = A*CS + B*SN;
            BB = -A*SN + B*CS;
            CC = C*CS + D*SN;
            DD = -C*SN + D*CS;

            // Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
                    // [ C  D ]   [-SN  CS ] [ CC  DD ]

            A = AA*CS + CC*SN;
            B = BB*CS + DD*SN;
            C = -AA*SN + CC*CS;
            D = -BB*SN + DD*CS;

            TEMP = HALF*( A+D );
            A = TEMP;
            D = TEMP;

            if ( C != ZERO ) {
               if ( B != ZERO ) {
                  if ( SIGN( ONE, B ) == SIGN( ONE, C ) ) {

                     // Real eigenvalues: reduce to upper triangular form

                     SAB = sqrt( ABS( B ) );
                     SAC = sqrt( ABS( C ) );
                     P = SIGN( SAB*SAC, C );
                     TAU = ONE / sqrt( ABS( B+C ) );
                     A = TEMP + P;
                     D = TEMP - P;
                     B = B - C;
                     C = ZERO;
                     CS1 = SAB*TAU;
                     SN1 = SAC*TAU;
                     TEMP = CS*CS1 - SN*SN1;
                     SN = CS*SN1 + SN*CS1;
                     CS = TEMP;
                  }
               } else {
                  B = -C;
                  C = ZERO;
                  TEMP = CS;
                  CS = -SN;
                  SN = TEMP;
               }
            }
         }

      }

      // Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).

      RT1R = A;
      RT2R = D;
      if ( C == ZERO ) {
         RT1I = ZERO;
         RT2I = ZERO;
      } else {
         RT1I = sqrt( ABS( B ) )*sqrt( ABS( C ) );
         RT2I = -RT1I;
      }
      return;

      // End of DLANV2

      }
