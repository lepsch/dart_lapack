      SUBROUTINE SPOTF2( UPLO, N, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J;
      REAL               AJJ
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      REAL               SDOT
      // EXTERNAL LSAME, SDOT, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('SPOTF2', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( UPPER ) {

         // Compute the Cholesky factorization A = U**T *U.

         DO 10 J = 1, N

            // Compute U(J,J) and test for non-positive-definiteness.

            AJJ = A( J, J ) - SDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
            if ( AJJ.LE.ZERO.OR.SISNAN( AJJ ) ) {
               A( J, J ) = AJJ
               GO TO 30
            }
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ

            // Compute elements J+1:N of row J.

            if ( J.LT.N ) {
               sgemv('Transpose', J-1, N-J, -ONE, A( 1, J+1 ), LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA );
               sscal(N-J, ONE / AJJ, A( J, J+1 ), LDA );
            }
   10    CONTINUE
      } else {

         // Compute the Cholesky factorization A = L*L**T.

         DO 20 J = 1, N

            // Compute L(J,J) and test for non-positive-definiteness.

            AJJ = A( J, J ) - SDOT( J-1, A( J, 1 ), LDA, A( J, 1 ), LDA )
            if ( AJJ.LE.ZERO.OR.SISNAN( AJJ ) ) {
               A( J, J ) = AJJ
               GO TO 30
            }
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ

            // Compute elements J+1:N of column J.

            if ( J.LT.N ) {
               sgemv('No transpose', N-J, J-1, -ONE, A( J+1, 1 ), LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 );
               sscal(N-J, ONE / AJJ, A( J+1, J ), 1 );
            }
   20    CONTINUE
      }
      GO TO 40

   30 CONTINUE
      INFO = J

   40 CONTINUE
      RETURN

      // End of SPOTF2

      }
