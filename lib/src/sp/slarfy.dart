      void slarfy(UPLO, N, V, INCV, TAU, C, LDC, WORK ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INCV, LDC, N;
      double               TAU;
      double               C( LDC, * ), V( * ), WORK( * );
      // ..

      double               ONE, ZERO, HALF;
      const              ONE = 1.0, ZERO = 0.0, HALF = 0.5 ;
      double               ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SSYMV, SSYR2
      // ..
      // .. External Functions ..
      //- REAL               SDOT;
      // EXTERNAL SDOT

      if (TAU == ZERO) return;

      // Form  w:= C * v

      ssymv(UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 );

      ALPHA = -HALF*TAU*SDOT( N, WORK, 1, V, INCV );
      saxpy(N, ALPHA, V, INCV, WORK, 1 );

      // C := C - v * w' - w * v'

      ssyr2(UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC );

      }
