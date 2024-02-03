      SUBROUTINE DLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int                INCV, LDC, N
      DOUBLE PRECISION   TAU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0, HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DSYMV, DSYR2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DDOT
      EXTERNAL           DDOT
*     ..
*     .. Executable Statements ..
*
      IF( TAU.EQ.ZERO ) RETURN
*
*     Form  w:= C * v
*
      CALL DSYMV( UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 )
*
      ALPHA = -HALF*TAU*DDOT( N, WORK, 1, V, INCV )
      CALL DAXPY( N, ALPHA, V, INCV, WORK, 1 )
*
*     C := C - v * w' - w * v'
*
      CALL DSYR2( UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC )
*
      RETURN
*
*     End of DLARFY
*
      END
