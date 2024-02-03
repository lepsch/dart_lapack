      int     FUNCTION ILACLC( M, N, A, LDA );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                M, N, LDA;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX          ZERO
      const     ZERO = (0.0E+0, 0.0E+0) ;
      // ..
      // .. Local Scalars ..
      int     I;
      // ..
      // .. Executable Statements ..

      // Quick test for the common case where one corner is non-zero.
      if ( N == 0 ) {
         ILACLC = N
      } else if ( A(1, N) != ZERO || A(M, N) != ZERO ) {
         ILACLC = N
      } else {
      // Now scan each column from the end, returning with the first non-zero.
         DO ILACLC = N, 1, -1
            for (I = 1; I <= M; I++) {
               IF( A(I, ILACLC) != ZERO ) RETURN
            }
         }
      }
      RETURN
      }
