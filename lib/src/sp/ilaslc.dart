      int ilaslc(M, N, A, LDA ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                M, N, LDA;
      // ..
      // .. Array Arguments ..
      double               A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const     ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int     I;
      // ..
      // .. Executable Statements ..

      // Quick test for the common case where one corner is non-zero.
      if ( N == 0 ) {
         ILASLC = N;
      } else if ( A(1, N) != ZERO || A(M, N) != ZERO ) {
         ILASLC = N;
      } else {
      // Now scan each column from the end, returning with the first non-zero.
         for (ILASLC = N; ILASLC >= 1; ILASLC--) {
            for (I = 1; I <= M; I++) {
               if( A(I, ILASLC) != ZERO ) return;
            }
         }
      }
      return;
      }