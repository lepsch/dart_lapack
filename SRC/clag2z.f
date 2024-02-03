      SUBROUTINE CLAG2Z( M, N, SA, LDSA, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDA, LDSA, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            SA( LDSA, * )
      COMPLEX*16         A( LDA, * )
      // ..
*
*  =====================================================================
*
      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. Executable Statements ..
*
      INFO = 0
      DO 20 J = 1, N
         DO 10 I = 1, M
            A( I, J ) = SA( I, J )
   10    CONTINUE
   20 CONTINUE
      RETURN
*
      // End of CLAG2Z
*
      END
