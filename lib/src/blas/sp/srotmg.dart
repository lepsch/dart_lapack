      void srotmg(SD1,SD2,SX1,SY1,SPARAM) {

// -- Reference BLAS level1 routine --
// -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double SD1,SD2,SX1,SY1;
      double SPARAM(5);
      // ..

// =====================================================================

      // .. Local Scalars ..
      double GAM,GAMSQ,ONE,RGAMSQ,SFLAG,SH11,SH12,SH21,SH22,SP1,SP2,SQ1, SQ2,STEMP,SU,TWO,ZERO;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Data statements ..

      final (ZERO,ONE,TWO) = (0.0,1.0,2.0);
      final (GAM,GAMSQ,RGAMSQ) = (4096.0,1.67772e7,5.96046e-8);
      // ..

      if (SD1 < ZERO) {
         // GO ZERO-H-D-AND-SX1..
         SFLAG = -ONE;
         SH11 = ZERO;
         SH12 = ZERO;
         SH21 = ZERO;
         SH22 = ZERO;

         SD1 = ZERO;
         SD2 = ZERO;
         SX1 = ZERO;
      } else {
         // CASE-SD1-NONNEGATIVE
         SP2 = SD2*SY1;
         if (SP2 == ZERO) {
            SFLAG = -TWO;
            SPARAM[1] = SFLAG;
            return;
         }
         // REGULAR-CASE..
         SP1 = SD1*SX1;
         SQ2 = SP2*SY1;
         SQ1 = SP1*SX1;

         if ((SQ1).abs() > (SQ2).abs()) {
            SH21 = -SY1/SX1;
            SH12 = SP2/SP1;

            SU = ONE - SH12*SH21;

           if (SU > ZERO) {
             SFLAG = ZERO;
             SD1 = SD1/SU;
             SD2 = SD2/SU;
             SX1 = SX1*SU;
           } else {
             // This code path if here for safety. We do not expect this
             // condition to ever hold except in edge cases with rounding
             // errors. See DOI: 10.1145/355841.355847
             SFLAG = -ONE;
             SH11 = ZERO;
             SH12 = ZERO;
             SH21 = ZERO;
             SH22 = ZERO;

             SD1 = ZERO;
             SD2 = ZERO;
             SX1 = ZERO;
           }
         } else {

            if (SQ2 < ZERO) {
               // GO ZERO-H-D-AND-SX1..
               SFLAG = -ONE;
               SH11 = ZERO;
               SH12 = ZERO;
               SH21 = ZERO;
               SH22 = ZERO;

               SD1 = ZERO;
               SD2 = ZERO;
               SX1 = ZERO;
            } else {
               SFLAG = ONE;
               SH11 = SP1/SP2;
               SH22 = SX1/SY1;
               SU = ONE + SH11*SH22;
               STEMP = SD2/SU;
               SD2 = SD1/SU;
               SD1 = STEMP;
               SX1 = SY1*SU;
            }
         }

      // PROCEDURE..SCALE-CHECK
         if (SD1 != ZERO) {
            while ((SD1 <= RGAMSQ) || (SD1 >= GAMSQ)) {
               if (SFLAG == ZERO) {
                  SH11 = ONE;
                  SH22 = ONE;
                  SFLAG = -ONE;
               } else {
                  SH21 = -ONE;
                  SH12 = ONE;
                  SFLAG = -ONE;
               }
               if (SD1 <= RGAMSQ) {
                  SD1 = SD1*GAM**2;
                  SX1 = SX1/GAM;
                  SH11 = SH11/GAM;
                  SH12 = SH12/GAM;
               } else {
                  SD1 = SD1/GAM**2;
                  SX1 = SX1*GAM;
                  SH11 = SH11*GAM;
                  SH12 = SH12*GAM;
               }
            }
         }

         if (SD2 != ZERO) {
            while (((SD2).abs() <= RGAMSQ) || ((SD2).abs() >= GAMSQ)) {
               if (SFLAG == ZERO) {
                  SH11 = ONE;
                  SH22 = ONE;
                  SFLAG = -ONE;
               } else {
                  SH21 = -ONE;
                  SH12 = ONE;
                  SFLAG = -ONE;
               }
               if ((SD2).abs() <= RGAMSQ) {
                  SD2 = SD2*GAM**2;
                  SH21 = SH21/GAM;
                  SH22 = SH22/GAM;
               } else {
                  SD2 = SD2/GAM**2;
                  SH21 = SH21*GAM;
                  SH22 = SH22*GAM;
               }
            }
         }

      }

      if (SFLAG < ZERO) {
         SPARAM[2] = SH11;
         SPARAM[3] = SH21;
         SPARAM[4] = SH12;
         SPARAM[5] = SH22;
      } else if (SFLAG == ZERO) {
         SPARAM[3] = SH21;
         SPARAM[4] = SH12;
      } else {
         SPARAM[2] = SH11;
         SPARAM[5] = SH22;
      }

      SPARAM[1] = SFLAG;
      }
