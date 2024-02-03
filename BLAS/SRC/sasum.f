      REAL FUNCTION SASUM(N,SX,INCX)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int     INCX,N;
      // ..
      // .. Array Arguments ..
      REAL SX(*)
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      REAL STEMP
      int     I,M,MP1,NINCX;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS,MOD
      // ..
      SASUM = 0.0e0
      STEMP = 0.0e0
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
         // code for increment equal to 1
*
*
         // clean-up loop
*
         M = MOD(N,6)
         IF (M.NE.0) THEN
            DO I = 1,M
               STEMP = STEMP + ABS(SX(I))
            END DO
            IF (N.LT.6) THEN
               SASUM = STEMP
               RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,6
            STEMP = STEMP + ABS(SX(I)) + ABS(SX(I+1)) + ABS(SX(I+2)) + ABS(SX(I+3)) + ABS(SX(I+4)) + ABS(SX(I+5))
         END DO
      ELSE
*
         // code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            STEMP = STEMP + ABS(SX(I))
         END DO
      END IF
      SASUM = STEMP
      RETURN
*
      // End of SASUM
*
      END
