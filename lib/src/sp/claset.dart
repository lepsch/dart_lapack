      void claset(final int UPLO, final int M, final int N, final int ALPHA, final int BETA, final int A, final int LDA,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, M, N;
      Complex            ALPHA, BETA;
      Complex            A( LDA, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN

      if ( lsame( UPLO, 'U' ) ) {

         // Set the diagonal to BETA and the strictly upper triangular
         // part of the array to ALPHA.

         for (J = 2; J <= N; J++) { // 20
            for (I = 1; I <= min( J-1, M ); I++) { // 10
               A[I][J] = ALPHA;
            } // 10
         } // 20
         for (I = 1; I <= min( N, M ); I++) { // 30
            A[I][I] = BETA;
         } // 30

      } else if ( lsame( UPLO, 'L' ) ) {

         // Set the diagonal to BETA and the strictly lower triangular
         // part of the array to ALPHA.

         for (J = 1; J <= min( M, N ); J++) { // 50
            for (I = J + 1; I <= M; I++) { // 40
               A[I][J] = ALPHA;
            } // 40
         } // 50
         for (I = 1; I <= min( N, M ); I++) { // 60
            A[I][I] = BETA;
         } // 60

      } else {

         // Set the array to BETA on the diagonal and ALPHA on the
         // offdiagonal.

         for (J = 1; J <= N; J++) { // 80
            for (I = 1; I <= M; I++) { // 70
               A[I][J] = ALPHA;
            } // 70
         } // 80
         for (I = 1; I <= min( M, N ); I++) { // 90
            A[I][I] = BETA;
         } // 90
      }

      }
