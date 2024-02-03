      SUBROUTINE DLAPLL( N, X, INCX, Y, INCY, SSMIN )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INCX, INCY, N
      DOUBLE PRECISION   SSMIN
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * ), Y( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   A11, A12, A22, C, SSMAX, TAU
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DLARFG, DLAS2
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.1 ) THEN
         SSMIN = ZERO
         RETURN
      END IF
*
*     Compute the QR factorization of the N-by-2 matrix ( X Y )
*
      CALL DLARFG( N, X( 1 ), X( 1+INCX ), INCX, TAU )
      A11 = X( 1 )
      X( 1 ) = ONE
*
      C = -TAU*DDOT( N, X, INCX, Y, INCY )
      CALL DAXPY( N, C, X, INCX, Y, INCY )
*
      CALL DLARFG( N-1, Y( 1+INCY ), Y( 1+2*INCY ), INCY, TAU )
*
      A12 = Y( 1 )
      A22 = Y( 1+INCY )
*
*     Compute the SVD of 2-by-2 Upper triangular matrix.
*
      CALL DLAS2( A11, A12, A22, SSMIN, SSMAX )
*
      RETURN
*
*     End of DLAPLL
*
      END
