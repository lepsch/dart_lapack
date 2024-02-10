      void dlarfy(UPLO, N, V, INCV, TAU, final Matrix<double> C, final int LDC, final Array<double> WORK) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INCV, LDC, N;
      double             TAU;
      double             C( LDC, * ), V( * ), WORK( * );
      // ..

      double             ONE, ZERO, HALF;
      const              ONE = 1.0, ZERO = 0.0, HALF = 0.5 ;
      double             ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DSYMV, DSYR2
      // ..
      // .. External Functions ..
      //- double             DDOT;
      // EXTERNAL DDOT

      if (TAU == ZERO) return;

      // Form  w:= C * v

      dsymv(UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 );

      ALPHA = -HALF*TAU*ddot( N, WORK, 1, V, INCV );
      daxpy(N, ALPHA, V, INCV, WORK, 1 );

      // C := C - v * w' - w * v'

      dsyr2(UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC );

      }
