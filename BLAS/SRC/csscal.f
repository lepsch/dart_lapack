      SUBROUTINE CSSCAL(N,SA,CX,INCX)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL SA
      int     INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX CX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      int     I,NINCX
*     ..
*     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E+0)
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC AIMAG,CMPLX,REAL
*     ..
      IF (N.LE.0 .OR. INCX.LE.0 .OR. SA.EQ.ONE) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         DO I = 1,N
            CX(I) = CMPLX(SA*REAL(CX(I)),SA*AIMAG(CX(I)))
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            CX(I) = CMPLX(SA*REAL(CX(I)),SA*AIMAG(CX(I)))
         END DO
      END IF
      RETURN
*
*     End of CSSCAL
*
      END
