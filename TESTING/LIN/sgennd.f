      bool    FUNCTION SGENND (M, N, A, LDA);

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     M, N, LDA;
      // ..
      // .. Array Arguments ..
      REAL A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      int     I, K;
      // ..
      // .. Intrinsics ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..
      K = MIN( M, N )
      for (I = 1; I <= K; I++) {
         if ( A( I, I ) < ZERO ) {
            SGENND = false;
            RETURN
         }
      }
      SGENND = true;
      RETURN
      }
