      SUBROUTINE STRTI2( UPLO, DIAG, N, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, UPPER;
      int                J;
      REAL               AJJ
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, STRMV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'STRTI2', -INFO )
         RETURN
      }

      if ( UPPER ) {

         // Compute inverse of upper triangular matrix.

         DO 10 J = 1, N
            if ( NOUNIT ) {
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            } else {
               AJJ = -ONE
            }

            // Compute elements 1:j-1 of j-th column.

            CALL STRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA, A( 1, J ), 1 )
            CALL SSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      } else {

         // Compute inverse of lower triangular matrix.

         DO 20 J = N, 1, -1
            if ( NOUNIT ) {
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            } else {
               AJJ = -ONE
            }
            if ( J.LT.N ) {

               // Compute elements j+1:n of j-th column.

               CALL STRMV( 'Lower', 'No transpose', DIAG, N-J, A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL SSCAL( N-J, AJJ, A( J+1, J ), 1 )
            }
   20    CONTINUE
      }

      RETURN

      // End of STRTI2

      }
