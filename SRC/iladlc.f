      int     FUNCTION ILADLC( M, N, A, LDA );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                M, N, LDA;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double           ZERO;
      const     ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int     I;
      // ..
      // .. Executable Statements ..

      // Quick test for the common case where one corner is non-zero.
      if ( N.EQ.0 ) {
         ILADLC = N
      } else if ( A(1, N).NE.ZERO .OR. A(M, N).NE.ZERO ) {
         ILADLC = N
      } else {
      // Now scan each column from the end, returning with the first non-zero.
         DO ILADLC = N, 1, -1
            for (I = 1; I <= M; I++) {
               IF( A(I, ILADLC).NE.ZERO ) RETURN
            }
         }
      }
      RETURN
      }
