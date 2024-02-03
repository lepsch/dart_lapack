      int     FUNCTION ILAZLR( M, N, A, LDA );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                M, N, LDA;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16       ZERO
      const     ZERO = (0.0D+0, 0.0D+0) ;
      // ..
      // .. Local Scalars ..
      int     I, J;
      // ..
      // .. Executable Statements ..

      // Quick test for the common case where one corner is non-zero.
      if ( M.EQ.0 ) {
         ILAZLR = M
      } else if ( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) {
         ILAZLR = M
      } else {
      // Scan up each column tracking the last zero row seen.
         ILAZLR = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ILAZLR = MAX( ILAZLR, I )
         END DO
      }
      RETURN
      }
