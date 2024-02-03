      bool zgennd(M, N, A, LDA) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     M, N, LDA;
      // ..
      // .. Array Arguments ..
      Complex A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int     I, K;
      Complex AII;
      // ..
      // .. Intrinsics ..
      // INTRINSIC MIN, DBLE, DIMAG
      // ..
      // .. Executable Statements ..
      K = min( M, N );
      for (I = 1; I <= K; I++) {
         AII = A( I, I );
         if ( DBLE( AII ) < ZERO || DIMAG( AII ) != ZERO ) {
            ZGENND = false;
            return;
         }
      }
      ZGENND = true;
      return;
      }
