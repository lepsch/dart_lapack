      int ilazlc(M, N, A, final int LDA) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                M, N, LDA;
      Complex         A( LDA, * );
      // ..

      Complex       ZERO;
      const     ZERO = (0.0, 0.0) ;
      int     I;

      // Quick test for the common case where one corner is non-zero.
      if ( N == 0 ) {
         ILAZLC = N;
      } else if ( A(1, N) != ZERO || A(M, N) != ZERO ) {
         ILAZLC = N;
      } else {
      // Now scan each column from the end, returning with the first non-zero.
         for (ILAZLC = N; ILAZLC >= 1; ILAZLC--) {
            for (I = 1; I <= M; I++) {
               if( A(I, ILAZLC) != ZERO ) return;
            }
         }
      }
      }
