      SUBROUTINE ZLARSCL2 ( M, N, D, X, LDX )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                M, N, LDX
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( LDX, * )
      DOUBLE PRECISION   D( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      int                I, J
*     ..
*     .. Executable Statements ..
*
      DO J = 1, N
         DO I = 1, M
            X( I, J ) = X( I, J ) / D( I )
         END DO
      END DO

      RETURN
      END
