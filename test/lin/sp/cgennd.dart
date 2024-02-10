      bool cgennd(M, N, A, final int LDA) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int     M, N, LDA;
      Complex A( LDA, * );
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      int     I, K;
      Complex AII;
      // ..
      // .. Intrinsics ..
      // INTRINSIC MIN, REAL, AIMAG
      K = min( M, N );
      for (I = 1; I <= K; I++) {
         AII = A( I, I );
         if ( double( AII ) < ZERO || AIMAG( AII ) != ZERO ) {
            CGENND = false;
            return;
         }
      }
      CGENND = true;
      }
