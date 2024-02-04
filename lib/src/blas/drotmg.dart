      void drotmg(DD1,DD2,DX1,DY1,DPARAM) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double           DD1,DD2,DX1,DY1;
      // ..
      // .. Array Arguments ..
      double           DPARAM(5);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double           DFLAG,DH11,DH12,DH21,DH22,DP1,DP2,DQ1,DQ2,DTEMP, DU,GAM,GAMSQ,ONE,RGAMSQ,TWO,ZERO;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DABS
      // ..
      // .. Data statements ..

      final (ZERO,ONE,TWO) = (0.0,1.0,2.0);
      final (GAM,GAMSQ,RGAMSQ) = (4096.0,16777216.0,5.9604645e-8);
      // ..

      if (DD1 < ZERO) {
         // GO ZERO-H-D-AND-DX1..
         DFLAG = -ONE;
         DH11 = ZERO;
         DH12 = ZERO;
         DH21 = ZERO;
         DH22 = ZERO;

         DD1 = ZERO;
         DD2 = ZERO;
         DX1 = ZERO;
      } else {
         // CASE-DD1-NONNEGATIVE
         DP2 = DD2*DY1;
         if (DP2 == ZERO) {
            DFLAG = -TWO;
            DPARAM[1] = DFLAG;
            return;
         }
         // REGULAR-CASE..
         DP1 = DD1*DX1;
         DQ2 = DP2*DY1;
         DQ1 = DP1*DX1;

         if ((DQ1).abs() > (DQ2).abs()) {
            DH21 = -DY1/DX1;
            DH12 = DP2/DP1;

            DU = ONE - DH12*DH21;

           if (DU > ZERO) {
             DFLAG = ZERO;
             DD1 = DD1/DU;
             DD2 = DD2/DU;
             DX1 = DX1*DU;
           } else {
             // This code path if here for safety. We do not expect this
             // condition to ever hold except in edge cases with rounding
             // errors. See DOI: 10.1145/355841.355847
             DFLAG = -ONE;
             DH11 = ZERO;
             DH12 = ZERO;
             DH21 = ZERO;
             DH22 = ZERO;

             DD1 = ZERO;
             DD2 = ZERO;
             DX1 = ZERO;
           }
         } else {

            if (DQ2 < ZERO) {
               // GO ZERO-H-D-AND-DX1..
               DFLAG = -ONE;
               DH11 = ZERO;
               DH12 = ZERO;
               DH21 = ZERO;
               DH22 = ZERO;

               DD1 = ZERO;
               DD2 = ZERO;
               DX1 = ZERO;
            } else {
               DFLAG = ONE;
               DH11 = DP1/DP2;
               DH22 = DX1/DY1;
               DU = ONE + DH11*DH22;
               DTEMP = DD2/DU;
               DD2 = DD1/DU;
               DD1 = DTEMP;
               DX1 = DY1*DU;
            }
         }

      // PROCEDURE..SCALE-CHECK
         if (DD1 != ZERO) {
            DO WHILE ((DD1 <= RGAMSQ) || (DD1 >= GAMSQ));
               if (DFLAG == ZERO) {
                  DH11 = ONE;
                  DH22 = ONE;
                  DFLAG = -ONE;
               } else {
                  DH21 = -ONE;
                  DH12 = ONE;
                  DFLAG = -ONE;
               }
               if (DD1 <= RGAMSQ) {
                  DD1 = DD1*GAM**2;
                  DX1 = DX1/GAM;
                  DH11 = DH11/GAM;
                  DH12 = DH12/GAM;
               } else {
                  DD1 = DD1/GAM**2;
                  DX1 = DX1*GAM;
                  DH11 = DH11*GAM;
                  DH12 = DH12*GAM;
               }
            }
         }

         if (DD2 != ZERO) {
            DO WHILE ( ((DD2).abs() <= RGAMSQ) || ((DD2).abs() >= GAMSQ) );
               if (DFLAG == ZERO) {
                  DH11 = ONE;
                  DH22 = ONE;
                  DFLAG = -ONE;
               } else {
                  DH21 = -ONE;
                  DH12 = ONE;
                  DFLAG = -ONE;
               }
               if ((DD2).abs() <= RGAMSQ) {
                  DD2 = DD2*GAM**2;
                  DH21 = DH21/GAM;
                  DH22 = DH22/GAM;
               } else {
                  DD2 = DD2/GAM**2;
                  DH21 = DH21*GAM;
                  DH22 = DH22*GAM;
               }
            }
         }

      }

      if (DFLAG < ZERO) {
         DPARAM[2] = DH11;
         DPARAM[3] = DH21;
         DPARAM[4] = DH12;
         DPARAM[5] = DH22;
      } else if (DFLAG == ZERO) {
         DPARAM[3] = DH21;
         DPARAM[4] = DH12;
      } else {
         DPARAM[2] = DH11;
         DPARAM[5] = DH22;
      }

      DPARAM[1] = DFLAG;
      return;
      }
