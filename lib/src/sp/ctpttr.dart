      void ctpttr(UPLO, N, AP, A, LDA, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N, LDA;
      // ..
      // .. Array Arguments ..
      Complex            A( LDA, * ), AP( * );
      // ..

// =====================================================================

      // .. Parameters ..
      // ..
      // .. Local Scalars ..
      bool               LOWER;
      int                I, J, K;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      LOWER = lsame( UPLO, 'L' );
      if ( !LOWER && !lsame( UPLO, 'U' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('CTPTTR', -INFO );
         return;
      }

      if ( LOWER ) {
         K = 0;
         for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
               K = K + 1;
               A[I][J] = AP( K );
            }
         }
      } else {
         K = 0;
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
               K = K + 1;
               A[I][J] = AP( K );
            }
         }
      }


      return;
      }
