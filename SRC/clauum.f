      SUBROUTINE CLAUUM( UPLO, N, A, LDA, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      COMPLEX            CONE
      const              CONE = ( 1.0E+0, 0.0E+0 ) ;
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
      // EXTERNAL CGEMM, CHERK, CLAUU2, CTRMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO != 0 ) {
         xerbla('CLAUUM', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'CLAUUM', UPLO, N, -1, -1, -1 )

      if ( NB.LE.1 || NB >= N ) {

         // Use unblocked code

         clauu2(UPLO, N, A, LDA, INFO );
      } else {

         // Use blocked code

         if ( UPPER ) {

            // Compute the product U * U**H.

            DO 10 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
               ctrmm('Right', 'Upper', 'Conjugate transpose', 'Non-unit', I-1, IB, CONE, A( I, I ), LDA, A( 1, I ), LDA );
               clauu2('Upper', IB, A( I, I ), LDA, INFO );
               if ( I+IB.LE.N ) {
                  cgemm('No transpose', 'Conjugate transpose', I-1, IB, N-I-IB+1, CONE, A( 1, I+IB ), LDA, A( I, I+IB ), LDA, CONE, A( 1, I ), LDA );
                  cherk('Upper', 'No transpose', IB, N-I-IB+1, ONE, A( I, I+IB ), LDA, ONE, A( I, I ), LDA );
               }
            } // 10
         } else {

            // Compute the product L**H * L.

            DO 20 I = 1, N, NB
               IB = MIN( NB, N-I+1 )
               ctrmm('Left', 'Lower', 'Conjugate transpose', 'Non-unit', IB, I-1, CONE, A( I, I ), LDA, A( I, 1 ), LDA );
               clauu2('Lower', IB, A( I, I ), LDA, INFO );
               if ( I+IB.LE.N ) {
                  cgemm('Conjugate transpose', 'No transpose', IB, I-1, N-I-IB+1, CONE, A( I+IB, I ), LDA, A( I+IB, 1 ), LDA, CONE, A( I, 1 ), LDA );
                  cherk('Lower', 'Conjugate transpose', IB, N-I-IB+1, ONE, A( I+IB, I ), LDA, ONE, A( I, I ), LDA );
               }
            } // 20
         }
      }

      RETURN

      // End of CLAUUM

      }
