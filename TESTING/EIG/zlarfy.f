      SUBROUTINE ZLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INCV, LDC, N;
      COMPLEX*16         TAU;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE, ZERO, HALF;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      COMPLEX*16         ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZHEMV, ZHER2
      // ..
      // .. External Functions ..
      COMPLEX*16         ZDOTC;
      // EXTERNAL ZDOTC
      // ..
      // .. Executable Statements ..

      if (TAU == ZERO) RETURN;

      // Form  w:= C * v

      zhemv(UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 );

      ALPHA = -HALF*TAU*ZDOTC( N, WORK, 1, V, INCV );
      zaxpy(N, ALPHA, V, INCV, WORK, 1 );

      // C := C - v * w' - w * v'

      zher2(UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC );

      return;

      // End of ZLARFY

      }
