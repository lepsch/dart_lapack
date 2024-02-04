import 'package:lapack/src/box.dart';
import 'package:lapack/src/matrix.dart';

void drotmg(final Box<double> DD1,final Box<double> DD2,final Box<double> DX1,final double DY1,final Array<double> DPARAM,) {
// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double           DFLAG,DH11=0,DH12=0,DH21=0,DH22=0,DP1,DP2,DQ1,DQ2,DTEMP, DU;

      final (ZERO,ONE,TWO) = (0.0,1.0,2.0);
      final (GAM,GAMSQ,RGAMSQ) = (4096.0,16777216.0,5.9604645e-8);

      if (DD1.value < ZERO) {
         // GO ZERO-H-D-AND-DX1.value..
         DFLAG = -ONE;
         DH11 = ZERO;
         DH12 = ZERO;
         DH21 = ZERO;
         DH22 = ZERO;

         DD1.value = ZERO;
         DD2.value = ZERO;
         DX1.value = ZERO;
      } else {
         // CASE-DD1.value-NONNEGATIVE
         DP2 = DD2.value*DY1;
         if (DP2 == ZERO) {
            DFLAG = -TWO;
            DPARAM[1] = DFLAG;
            return;
         }
         // REGULAR-CASE..
         DP1 = DD1.value*DX1.value;
         DQ2 = DP2*DY1;
         DQ1 = DP1*DX1.value;

         if ((DQ1).abs() > (DQ2).abs()) {
            DH21 = -DY1/DX1.value;
            DH12 = DP2/DP1;

            DU = ONE - DH12*DH21;

           if (DU > ZERO) {
             DFLAG = ZERO;
             DD1.value = DD1.value/DU;
             DD2.value = DD2.value/DU;
             DX1.value = DX1.value*DU;
           } else {
             // This code path if here for safety. We do not expect this
             // condition to ever hold except in edge cases with rounding
             // errors. See DOI: 10.1145/355841.355847
             DFLAG = -ONE;
             DH11 = ZERO;
             DH12 = ZERO;
             DH21 = ZERO;
             DH22 = ZERO;

             DD1.value = ZERO;
             DD2.value = ZERO;
             DX1.value = ZERO;
           }
         } else {

            if (DQ2 < ZERO) {
               // GO ZERO-H-D-AND-DX1.value..
               DFLAG = -ONE;
               DH11 = ZERO;
               DH12 = ZERO;
               DH21 = ZERO;
               DH22 = ZERO;

               DD1.value = ZERO;
               DD2.value = ZERO;
               DX1.value = ZERO;
            } else {
               DFLAG = ONE;
               DH11 = DP1/DP2;
               DH22 = DX1.value/DY1;
               DU = ONE + DH11*DH22;
               DTEMP = DD2.value/DU;
               DD2.value = DD1.value/DU;
               DD1.value = DTEMP;
               DX1.value = DY1*DU;
            }
         }

      // PROCEDURE..SCALE-CHECK
         if (DD1.value != ZERO) {
            while ((DD1.value <= RGAMSQ) || (DD1.value >= GAMSQ)) {
               if (DFLAG == ZERO) {
                  DH11 = ONE;
                  DH22 = ONE;
                  DFLAG = -ONE;
               } else {
                  DH21 = -ONE;
                  DH12 = ONE;
                  DFLAG = -ONE;
               }
               if (DD1.value <= RGAMSQ) {
                  DD1.value = DD1.value*GAM**2;
                  DX1.value = DX1.value/GAM;
                  DH11 = DH11/GAM;
                  DH12 = DH12/GAM;
               } else {
                  DD1.value = DD1.value/GAM**2;
                  DX1.value = DX1.value*GAM;
                  DH11 = DH11*GAM;
                  DH12 = DH12*GAM;
               }
            }
         }

         if (DD2.value != ZERO) {
            while ( ((DD2.value).abs() <= RGAMSQ) || ((DD2.value).abs() >= GAMSQ) ) {
               if (DFLAG == ZERO) {
                  DH11 = ONE;
                  DH22 = ONE;
                  DFLAG = -ONE;
               } else {
                  DH21 = -ONE;
                  DH12 = ONE;
                  DFLAG = -ONE;
               }
               if ((DD2.value).abs() <= RGAMSQ) {
                  DD2.value = DD2.value*GAM**2;
                  DH21 = DH21/GAM;
                  DH22 = DH22/GAM;
               } else {
                  DD2.value = DD2.value/GAM**2;
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
