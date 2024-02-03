      SUBROUTINE SLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INCV, LDC, N;
      REAL               TAU
      // ..
      // .. Array Arguments ..
      REAL               C( LDC, * ), V( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO, HALF
      const              ONE = 1.0E+0, ZERO = 0.0E+0, HALF = 0.5E+0 ;
      // ..
      // .. Local Scalars ..
      REAL               ALPHA
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SSYMV, SSYR2
      // ..
      // .. External Functions ..
      REAL               SDOT
      // EXTERNAL SDOT
      // ..
      // .. Executable Statements ..

      IF( TAU.EQ.ZERO ) RETURN

      // Form  w:= C * v

      CALL SSYMV( UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 )

      ALPHA = -HALF*TAU*SDOT( N, WORK, 1, V, INCV )
      CALL SAXPY( N, ALPHA, V, INCV, WORK, 1 )

      // C := C - v * w' - w * v'

      CALL SSYR2( UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC )

      RETURN

      // End of SLARFY

      }
