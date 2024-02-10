      bool zgennd(M, N, A, final int LDA) {

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
      // INTRINSIC MIN, DBLE, DIMAG
      K = min( M, N );
      for (I = 1; I <= K; I++) {
         AII = A( I, I );
         if ( AII.toDouble() < ZERO || DIMAG( AII ) != ZERO ) {
            ZGENND = false;
            return;
         }
      }
      ZGENND = true;
      }
