      int     FUNCTION ILACLR( M, N, A, LDA );

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
      int     I, J;
      // ..
      // .. Executable Statements ..

      // Quick test for the common case where one corner is non-zero.
      if ( M == 0 ) {
         ILACLR = M
      } else if ( A(M, 1) != ZERO .OR. A(M, N) != ZERO ) {
         ILACLR = M
      } else {
      // Scan up each column tracking the last zero row seen.
         ILACLR = 0
         for (J = 1; J <= N; J++) {
            I=M
            DO WHILE((A(MAX(I,1),J) == ZERO).AND.(I.GE.1))
               I=I-1
            }
            ILACLR = MAX( ILACLR, I )
         }
      }
      RETURN
      }
