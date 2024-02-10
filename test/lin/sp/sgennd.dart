      bool sgennd(final int M, final int N, final int A, final int LDA) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     M, N, LDA;
      double A( LDA, * );
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      int     I, K;
      // ..
      // .. Intrinsics ..
      // INTRINSIC MIN
      K = min( M, N );
      for (I = 1; I <= K; I++) {
         if ( A( I, I ) < ZERO ) {
            SGENND = false;
            return;
         }
      }
      SGENND = true;
      }
