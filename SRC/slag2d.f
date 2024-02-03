      SUBROUTINE SLAG2D( M, N, SA, LDSA, A, LDA, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDSA, M, N;
      // ..
      // .. Array Arguments ..
      REAL               SA( LDSA, * )
      double             A( LDA, * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. Executable Statements ..

      INFO = 0
      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            A( I, J ) = SA( I, J )
         } // 10
      } // 20
      RETURN

      // End of SLAG2D

      }
