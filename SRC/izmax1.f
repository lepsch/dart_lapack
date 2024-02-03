      int     FUNCTION IZMAX1( N, ZX, INCX )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INCX, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         ZX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      double             DMAX;
      int                I, IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IZMAX1 = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IZMAX1 = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         DMAX = ABS(ZX(1))
         DO I = 2,N
            IF (ABS(ZX(I)).GT.DMAX) THEN
               IZMAX1 = I
               DMAX = ABS(ZX(I))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         IX = 1
         DMAX = ABS(ZX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (ABS(ZX(IX)).GT.DMAX) THEN
               IZMAX1 = I
               DMAX = ABS(ZX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
*
*     End of IZMAX1
*
      END
