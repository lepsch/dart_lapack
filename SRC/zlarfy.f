      SUBROUTINE ZLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INCV, LDC, N;
      COMPLEX*16         TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO, HALF
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ), HALF = ( 0.5D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZHEMV, ZHER2
*     ..
*     .. External Functions ..
      COMPLEX*16         ZDOTC
      EXTERNAL           ZDOTC
*     ..
*     .. Executable Statements ..
*
      IF( TAU.EQ.ZERO ) RETURN
*
*     Form  w:= C * v
*
      CALL ZHEMV( UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 )
*
      ALPHA = -HALF*TAU*ZDOTC( N, WORK, 1, V, INCV )
      CALL ZAXPY( N, ALPHA, V, INCV, WORK, 1 )
*
*     C := C - v * w' - w * v'
*
      CALL ZHER2( UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC )
*
      RETURN
*
*     End of ZLARFY
*
      END
