      void spbtf2(final int UPLO, final int N, final int KD, final Matrix<double> AB, final int LDAB, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, KD, LDAB, N;
      double               AB( LDAB, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               UPPER;
      int                J, KLD, KN;
      double               AJJ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSYR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KD < 0 ) {
         INFO = -3;
      } else if ( LDAB < KD+1 ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('SPBTF2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      KLD = max( 1, LDAB-1 );

      if ( UPPER ) {

         // Compute the Cholesky factorization A = U**T*U.

         for (J = 1; J <= N; J++) { // 10

            // Compute U(J,J) and test for non-positive-definiteness.

            AJJ = AB( KD+1, J );
            if (AJJ <= ZERO) GO TO 30;
            AJJ = sqrt( AJJ );
            AB[KD+1][J] = AJJ;

            // Compute elements J+1:J+KN of row J and update the
            // trailing submatrix within the band.

            KN = min( KD, N-J );
            if ( KN > 0 ) {
               sscal(KN, ONE / AJJ, AB( KD, J+1 ), KLD );
               ssyr('Upper', KN, -ONE, AB( KD, J+1 ), KLD, AB( KD+1, J+1 ), KLD );
            }
         } // 10
      } else {

         // Compute the Cholesky factorization A = L*L**T.

         for (J = 1; J <= N; J++) { // 20

            // Compute L(J,J) and test for non-positive-definiteness.

            AJJ = AB( 1, J );
            if (AJJ <= ZERO) GO TO 30;
            AJJ = sqrt( AJJ );
            AB[1][J] = AJJ;

            // Compute elements J+1:J+KN of column J and update the
            // trailing submatrix within the band.

            KN = min( KD, N-J );
            if ( KN > 0 ) {
               sscal(KN, ONE / AJJ, AB( 2, J ), 1 );
               ssyr('Lower', KN, -ONE, AB( 2, J ), 1, AB( 1, J+1 ), KLD );
            }
         } // 20
      }
      return;

      } // 30
      INFO = J;
      }
