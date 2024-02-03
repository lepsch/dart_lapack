      int     FUNCTION ILAZLC( M, N, A, LDA );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                M, N, LDA;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX*16       ZERO;
      const     ZERO = (0.0, 0.0) ;
      // ..
      // .. Local Scalars ..
      int     I;
      // ..
      // .. Executable Statements ..

      // Quick test for the common case where one corner is non-zero.
      if ( N == 0 ) {
         ILAZLC = N;
      } else if ( A(1, N) != ZERO || A(M, N) != ZERO ) {
         ILAZLC = N;
      } else {
      // Now scan each column from the end, returning with the first non-zero.
         DO ILAZLC = N, 1, -1;
            for (I = 1; I <= M; I++) {
               IF( A(I, ILAZLC) != ZERO ) RETURN;
            }
         }
      }
      return;
      }
