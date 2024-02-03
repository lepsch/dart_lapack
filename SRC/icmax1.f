      int     FUNCTION ICMAX1( N, CX, INCX );
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INCX, N;
*     ..
*     .. Array Arguments ..
      COMPLEX            CX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL               SMAX
      int                I, IX;
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS
*     ..
*     .. Executable Statements ..
*
      ICMAX1 = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      ICMAX1 = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         SMAX = ABS(CX(1))
         DO I = 2,N
            IF (ABS(CX(I)).GT.SMAX) THEN
               ICMAX1 = I
               SMAX = ABS(CX(I))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         IX = 1
         SMAX = ABS(CX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (ABS(CX(IX)).GT.SMAX) THEN
               ICMAX1 = I
               SMAX = ABS(CX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
*
*     End of ICMAX1
*
      END
