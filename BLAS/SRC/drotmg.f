      SUBROUTINE DROTMG(DD1,DD2,DX1,DY1,DPARAM)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      double           DD1,DD2,DX1,DY1;
*     ..
*     .. Array Arguments ..
      double           DPARAM(5);
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      double           DFLAG,DH11,DH12,DH21,DH22,DP1,DP2,DQ1,DQ2,DTEMP, DU,GAM,GAMSQ,ONE,RGAMSQ,TWO,ZERO;
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DABS
*     ..
*     .. Data statements ..
*
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/
      DATA GAM,GAMSQ,RGAMSQ/4096.D0,16777216.D0,5.9604645D-8/
*     ..

      IF (DD1.LT.ZERO) THEN
*        GO ZERO-H-D-AND-DX1..
         DFLAG = -ONE
         DH11 = ZERO
         DH12 = ZERO
         DH21 = ZERO
         DH22 = ZERO
*
         DD1 = ZERO
         DD2 = ZERO
         DX1 = ZERO
      ELSE
*        CASE-DD1-NONNEGATIVE
         DP2 = DD2*DY1
         IF (DP2.EQ.ZERO) THEN
            DFLAG = -TWO
            DPARAM(1) = DFLAG
            RETURN
         END IF
*        REGULAR-CASE..
         DP1 = DD1*DX1
         DQ2 = DP2*DY1
         DQ1 = DP1*DX1
*
         IF (DABS(DQ1).GT.DABS(DQ2)) THEN
            DH21 = -DY1/DX1
            DH12 = DP2/DP1
*
            DU = ONE - DH12*DH21
*
           IF (DU.GT.ZERO) THEN
             DFLAG = ZERO
             DD1 = DD1/DU
             DD2 = DD2/DU
             DX1 = DX1*DU
           ELSE
*            This code path if here for safety. We do not expect this
*            condition to ever hold except in edge cases with rounding
*            errors. See DOI: 10.1145/355841.355847
             DFLAG = -ONE
             DH11 = ZERO
             DH12 = ZERO
             DH21 = ZERO
             DH22 = ZERO
*
             DD1 = ZERO
             DD2 = ZERO
             DX1 = ZERO
           END IF
         ELSE

            IF (DQ2.LT.ZERO) THEN
*              GO ZERO-H-D-AND-DX1..
               DFLAG = -ONE
               DH11 = ZERO
               DH12 = ZERO
               DH21 = ZERO
               DH22 = ZERO
*
               DD1 = ZERO
               DD2 = ZERO
               DX1 = ZERO
            ELSE
               DFLAG = ONE
               DH11 = DP1/DP2
               DH22 = DX1/DY1
               DU = ONE + DH11*DH22
               DTEMP = DD2/DU
               DD2 = DD1/DU
               DD1 = DTEMP
               DX1 = DY1*DU
            END IF
         END IF

*     PROCEDURE..SCALE-CHECK
         IF (DD1.NE.ZERO) THEN
            DO WHILE ((DD1.LE.RGAMSQ) .OR. (DD1.GE.GAMSQ))
               IF (DFLAG.EQ.ZERO) THEN
                  DH11 = ONE
                  DH22 = ONE
                  DFLAG = -ONE
               ELSE
                  DH21 = -ONE
                  DH12 = ONE
                  DFLAG = -ONE
               END IF
               IF (DD1.LE.RGAMSQ) THEN
                  DD1 = DD1*GAM**2
                  DX1 = DX1/GAM
                  DH11 = DH11/GAM
                  DH12 = DH12/GAM
               ELSE
                  DD1 = DD1/GAM**2
                  DX1 = DX1*GAM
                  DH11 = DH11*GAM
                  DH12 = DH12*GAM
               END IF
            ENDDO
         END IF

         IF (DD2.NE.ZERO) THEN
            DO WHILE ( (DABS(DD2).LE.RGAMSQ) .OR. (DABS(DD2).GE.GAMSQ) )
               IF (DFLAG.EQ.ZERO) THEN
                  DH11 = ONE
                  DH22 = ONE
                  DFLAG = -ONE
               ELSE
                  DH21 = -ONE
                  DH12 = ONE
                  DFLAG = -ONE
               END IF
               IF (DABS(DD2).LE.RGAMSQ) THEN
                  DD2 = DD2*GAM**2
                  DH21 = DH21/GAM
                  DH22 = DH22/GAM
               ELSE
                  DD2 = DD2/GAM**2
                  DH21 = DH21*GAM
                  DH22 = DH22*GAM
               END IF
            END DO
         END IF

      END IF

      IF (DFLAG.LT.ZERO) THEN
         DPARAM(2) = DH11
         DPARAM(3) = DH21
         DPARAM(4) = DH12
         DPARAM(5) = DH22
      ELSE IF (DFLAG.EQ.ZERO) THEN
         DPARAM(3) = DH21
         DPARAM(4) = DH12
      ELSE
         DPARAM(2) = DH11
         DPARAM(5) = DH22
      END IF

      DPARAM(1) = DFLAG
      RETURN
*
*     End of DROTMG
*
      END
