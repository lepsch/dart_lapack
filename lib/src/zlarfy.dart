      void zlarfy(UPLO, N, V, INCV, TAU, C, LDC, WORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INCV, LDC, N;
      Complex         TAU;
      Complex         C( LDC, * ), V( * ), WORK( * );
      // ..

      Complex         ONE, ZERO, HALF;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      Complex         ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZHEMV, ZHER2
      // ..
      // .. External Functions ..
      //- Complex         ZDOTC;
      // EXTERNAL ZDOTC

      if (TAU == ZERO) return;

      // Form  w:= C * v

      zhemv(UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 );

      ALPHA = -HALF*TAU*ZDOTC( N, WORK, 1, V, INCV );
      zaxpy(N, ALPHA, V, INCV, WORK, 1 );

      // C := C - v * w' - w * v'

      zher2(UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC );

      return;
      }
