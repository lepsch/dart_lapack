      void ztrtri(UPLO, DIAG, N, A, LDA, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, UPPER;
      int                J, JB, NB, NN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZTRMM, ZTRSM, ZTRTI2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      NOUNIT = LSAME( DIAG, 'N' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( !NOUNIT && !LSAME( DIAG, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('ZTRTRI', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Check for singularity if non-unit.

      if ( NOUNIT ) {
         for (INFO = 1; INFO <= N; INFO++) { // 10
            if( A( INFO, INFO ) == ZERO ) return;
         } // 10
         INFO = 0;
      }

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'ZTRTRI', UPLO // DIAG, N, -1, -1, -1 );
      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code

         ztrti2(UPLO, DIAG, N, A, LDA, INFO );
      } else {

         // Use blocked code

         if ( UPPER ) {

            // Compute inverse of upper triangular matrix

            DO 20 J = 1, N, NB;
               JB = min( NB, N-J+1 );

               // Compute rows 1:j-1 of current block column

               ztrmm('Left', 'Upper', 'No transpose', DIAG, J-1, JB, ONE, A, LDA, A( 1, J ), LDA );
               ztrsm('Right', 'Upper', 'No transpose', DIAG, J-1, JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA );

               // Compute inverse of current diagonal block

               ztrti2('Upper', DIAG, JB, A( J, J ), LDA, INFO );
            } // 20
         } else {

            // Compute inverse of lower triangular matrix

            NN = ( ( N-1 ) / NB )*NB + 1;
            DO 30 J = NN, 1, -NB;
               JB = min( NB, N-J+1 );
               if ( J+JB <= N ) {

                  // Compute rows j+jb:n of current block column

                  ztrmm('Left', 'Lower', 'No transpose', DIAG, N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA, A( J+JB, J ), LDA );
                  ztrsm('Right', 'Lower', 'No transpose', DIAG, N-J-JB+1, JB, -ONE, A( J, J ), LDA, A( J+JB, J ), LDA );
               }

               // Compute inverse of current diagonal block

               ztrti2('Lower', DIAG, JB, A( J, J ), LDA, INFO );
            } // 30
         }
      }

      return;

      // End of ZTRTRI

      }
