      SUBROUTINE SLAUU2( UPLO, N, A, LDA, INFO )

*  -- LAPACK auxiliary routine --
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
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I;
      REAL               AII
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SDOT
      // EXTERNAL LSAME, SDOT
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
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
         xerbla('SLAUU2', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( UPPER ) {

         // Compute the product U * U**T.

         DO 10 I = 1, N
            AII = A( I, I )
            if ( I.LT.N ) {
               A( I, I ) = SDOT( N-I+1, A( I, I ), LDA, A( I, I ), LDA )
               sgemv('No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, AII, A( 1, I ), 1 );
            } else {
               sscal(I, AII, A( 1, I ), 1 );
            }
   10    CONTINUE

      } else {

         // Compute the product L**T * L.

         DO 20 I = 1, N
            AII = A( I, I )
            if ( I.LT.N ) {
               A( I, I ) = SDOT( N-I+1, A( I, I ), 1, A( I, I ), 1 )
               sgemv('Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, AII, A( I, 1 ), LDA );
            } else {
               sscal(I, AII, A( I, 1 ), LDA );
            }
   20    CONTINUE
      }

      RETURN

      // End of SLAUU2

      }
