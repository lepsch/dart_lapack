      SUBROUTINE DLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INCV, LDC, N;
      double             TAU;
      // ..
      // .. Array Arguments ..
      double             C( LDC, * ), V( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO, HALF;
      const              ONE = 1.0D+0, ZERO = 0.0D+0, HALF = 0.5D+0 ;
      // ..
      // .. Local Scalars ..
      double             ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DSYMV, DSYR2
      // ..
      // .. External Functions ..
      double             DDOT;
      // EXTERNAL DDOT
      // ..
      // .. Executable Statements ..

      if (TAU.EQ.ZERO) RETURN;

      // Form  w:= C * v

      dsymv(UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 );

      ALPHA = -HALF*TAU*DDOT( N, WORK, 1, V, INCV )
      daxpy(N, ALPHA, V, INCV, WORK, 1 );

      // C := C - v * w' - w * v'

      dsyr2(UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC );

      RETURN

      // End of DLARFY

      }
