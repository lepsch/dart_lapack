      bool    FUNCTION CGENND (M, N, A, LDA);

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     M, N, LDA;
      // ..
      // .. Array Arguments ..
      COMPLEX A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      int     I, K;
      COMPLEX AII
      // ..
      // .. Intrinsics ..
      // INTRINSIC MIN, REAL, AIMAG
      // ..
      // .. Executable Statements ..
      K = MIN( M, N )
      for (I = 1; I <= K; I++) {
         AII = A( I, I )
         if ( REAL( AII ) < ZERO || AIMAG( AII ) != ZERO ) {
            CGENND = false;
            RETURN
         }
      }
      CGENND = true;
      RETURN
      }
