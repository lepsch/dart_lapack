      SUBROUTINE SROTMG(SD1,SD2,SX1,SY1,SPARAM)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL SD1,SD2,SX1,SY1
*     ..
*     .. Array Arguments ..
      REAL SPARAM(5)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL GAM,GAMSQ,ONE,RGAMSQ,SFLAG,SH11,SH12,SH21,SH22,SP1,SP2,SQ1, SQ2,STEMP,SU,TWO,ZERO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS
*     ..
*     .. Data statements ..
*
      DATA ZERO,ONE,TWO/0.E0,1.E0,2.E0/
      DATA GAM,GAMSQ,RGAMSQ/4096.E0,1.67772E7,5.96046E-8/
*     ..

      IF (SD1.LT.ZERO) THEN
*        GO ZERO-H-D-AND-SX1..
         SFLAG = -ONE
         SH11 = ZERO
         SH12 = ZERO
         SH21 = ZERO
         SH22 = ZERO
*
         SD1 = ZERO
         SD2 = ZERO
         SX1 = ZERO
      ELSE
*        CASE-SD1-NONNEGATIVE
         SP2 = SD2*SY1
         IF (SP2.EQ.ZERO) THEN
            SFLAG = -TWO
            SPARAM(1) = SFLAG
            RETURN
         END IF
*        REGULAR-CASE..
         SP1 = SD1*SX1
         SQ2 = SP2*SY1
         SQ1 = SP1*SX1
*
         IF (ABS(SQ1).GT.ABS(SQ2)) THEN
            SH21 = -SY1/SX1
            SH12 = SP2/SP1
*
            SU = ONE - SH12*SH21
*
           IF (SU.GT.ZERO) THEN
             SFLAG = ZERO
             SD1 = SD1/SU
             SD2 = SD2/SU
             SX1 = SX1*SU
           ELSE
*            This code path if here for safety. We do not expect this
*            condition to ever hold except in edge cases with rounding
*            errors. See DOI: 10.1145/355841.355847
             SFLAG = -ONE
             SH11 = ZERO
             SH12 = ZERO
             SH21 = ZERO
             SH22 = ZERO
*
             SD1 = ZERO
             SD2 = ZERO
             SX1 = ZERO
           END IF
         ELSE

            IF (SQ2.LT.ZERO) THEN
*              GO ZERO-H-D-AND-SX1..
               SFLAG = -ONE
               SH11 = ZERO
               SH12 = ZERO
               SH21 = ZERO
               SH22 = ZERO
*
               SD1 = ZERO
               SD2 = ZERO
               SX1 = ZERO
            ELSE
               SFLAG = ONE
               SH11 = SP1/SP2
               SH22 = SX1/SY1
               SU = ONE + SH11*SH22
               STEMP = SD2/SU
               SD2 = SD1/SU
               SD1 = STEMP
               SX1 = SY1*SU
            END IF
         END IF

*     PROCEDURE..SCALE-CHECK
         IF (SD1.NE.ZERO) THEN
            DO WHILE ((SD1.LE.RGAMSQ) .OR. (SD1.GE.GAMSQ))
               IF (SFLAG.EQ.ZERO) THEN
                  SH11 = ONE
                  SH22 = ONE
                  SFLAG = -ONE
               ELSE
                  SH21 = -ONE
                  SH12 = ONE
                  SFLAG = -ONE
               END IF
               IF (SD1.LE.RGAMSQ) THEN
                  SD1 = SD1*GAM**2
                  SX1 = SX1/GAM
                  SH11 = SH11/GAM
                  SH12 = SH12/GAM
               ELSE
                  SD1 = SD1/GAM**2
                  SX1 = SX1*GAM
                  SH11 = SH11*GAM
                  SH12 = SH12*GAM
               END IF
            ENDDO
         END IF

         IF (SD2.NE.ZERO) THEN
            DO WHILE ( (ABS(SD2).LE.RGAMSQ) .OR. (ABS(SD2).GE.GAMSQ) )
               IF (SFLAG.EQ.ZERO) THEN
                  SH11 = ONE
                  SH22 = ONE
                  SFLAG = -ONE
               ELSE
                  SH21 = -ONE
                  SH12 = ONE
                  SFLAG = -ONE
               END IF
               IF (ABS(SD2).LE.RGAMSQ) THEN
                  SD2 = SD2*GAM**2
                  SH21 = SH21/GAM
                  SH22 = SH22/GAM
               ELSE
                  SD2 = SD2/GAM**2
                  SH21 = SH21*GAM
                  SH22 = SH22*GAM
               END IF
            END DO
         END IF

      END IF

      IF (SFLAG.LT.ZERO) THEN
         SPARAM(2) = SH11
         SPARAM(3) = SH21
         SPARAM(4) = SH12
         SPARAM(5) = SH22
      ELSE IF (SFLAG.EQ.ZERO) THEN
         SPARAM(3) = SH21
         SPARAM(4) = SH12
      ELSE
         SPARAM(2) = SH11
         SPARAM(5) = SH22
      END IF

      SPARAM(1) = SFLAG
      RETURN
*
*     End of SROTMG
*
      END
