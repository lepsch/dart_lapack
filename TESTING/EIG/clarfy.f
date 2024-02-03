      SUBROUTINE CLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INCV, LDC, N;
      COMPLEX            TAU
      // ..
      // .. Array Arguments ..
      COMPLEX            C( LDC, * ), V( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, ZERO, HALF
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      COMPLEX            ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CHEMV, CHER2
      // ..
      // .. External Functions ..
      COMPLEX            CDOTC
      // EXTERNAL CDOTC
      // ..
      // .. Executable Statements ..

      if (TAU == ZERO) RETURN;

      // Form  w:= C * v

      chemv(UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 );

      ALPHA = -HALF*TAU*CDOTC( N, WORK, 1, V, INCV )
      caxpy(N, ALPHA, V, INCV, WORK, 1 );

      // C := C - v * w' - w * v'

      cher2(UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC );

      RETURN

      // End of CLARFY

      }
