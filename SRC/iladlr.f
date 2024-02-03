      int     FUNCTION ILADLR( M, N, A, LDA );

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
      const     ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int     I, J;
      // ..
      // .. Executable Statements ..

      // Quick test for the common case where one corner is non-zero.
      if ( M == 0 ) {
         ILADLR = M
      } else if ( A(M, 1) != ZERO || A(M, N) != ZERO ) {
         ILADLR = M
      } else {
      // Scan up each column tracking the last zero row seen.
         ILADLR = 0
         for (J = 1; J <= N; J++) {
            I=M
            DO WHILE((A(MAX(I,1),J) == ZERO) && (I >= 1))
               I=I-1
            }
            ILADLR = MAX( ILADLR, I )
         }
      }
      RETURN
      }
