      int     FUNCTION ILASLR( M, N, A, LDA );
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                M, N, LDA;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL             ZERO
      PARAMETER ( ZERO = 0.0E+0 )
      // ..
      // .. Local Scalars ..
      int     I, J;
      // ..
      // .. Executable Statements ..
*
      // Quick test for the common case where one corner is non-zero.
      IF( M.EQ.0 ) THEN
         ILASLR = M
      ELSEIF( A(M, 1).NE.ZERO .OR. A(M, N).NE.ZERO ) THEN
         ILASLR = M
      ELSE
      // Scan up each column tracking the last zero row seen.
         ILASLR = 0
         DO J = 1, N
            I=M
            DO WHILE((A(MAX(I,1),J).EQ.ZERO).AND.(I.GE.1))
               I=I-1
            ENDDO
            ILASLR = MAX( ILASLR, I )
         END DO
      END IF
      RETURN
      END
