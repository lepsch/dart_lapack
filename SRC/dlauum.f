      SUBROUTINE DLAUUM( UPLO, N, A, LDA, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, IB, NB;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLAUU2, DSYRK, DTRMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
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
         xerbla('DLAUUM', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0) RETURN;

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'DLAUUM', UPLO, N, -1, -1, -1 )

      if ( NB.LE.1 .OR. NB.GE.N ) {

         // Use unblocked code

         dlauu2(UPLO, N, A, LDA, INFO );
      } else {

         // Use blocked code

         if ( UPPER ) {

            // Compute the product U * U**T.

            DO 10 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
               dtrmm('Right', 'Upper', 'Transpose', 'Non-unit', I-1, IB, ONE, A( I, I ), LDA, A( 1, I ), LDA );
               dlauu2('Upper', IB, A( I, I ), LDA, INFO );
               if ( I+IB.LE.N ) {
                  dgemm('No transpose', 'Transpose', I-1, IB, N-I-IB+1, ONE, A( 1, I+IB ), LDA, A( I, I+IB ), LDA, ONE, A( 1, I ), LDA );
                  dsyrk('Upper', 'No transpose', IB, N-I-IB+1, ONE, A( I, I+IB ), LDA, ONE, A( I, I ), LDA );
               }
            } // 10
         } else {

            // Compute the product L**T * L.

            DO 20 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
               dtrmm('Left', 'Lower', 'Transpose', 'Non-unit', IB, I-1, ONE, A( I, I ), LDA, A( I, 1 ), LDA );
               dlauu2('Lower', IB, A( I, I ), LDA, INFO );
               if ( I+IB.LE.N ) {
                  dgemm('Transpose', 'No transpose', IB, I-1, N-I-IB+1, ONE, A( I+IB, I ), LDA, A( I+IB, 1 ), LDA, ONE, A( I, 1 ), LDA );
                  dsyrk('Lower', 'Transpose', IB, N-I-IB+1, ONE, A( I+IB, I ), LDA, ONE, A( I, I ), LDA );
               }
            } // 20
         }
      }

      RETURN

      // End of DLAUUM

      }
