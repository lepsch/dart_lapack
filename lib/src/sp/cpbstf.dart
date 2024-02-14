      void cpbstf(final int UPLO, final int N, final int KD, final Matrix<double> AB_, final int LDAB, final Box<int> INFO,) {
  final AB = AB_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, KD, LDAB, N;
      Complex            AB( LDAB, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               UPPER;
      int                J, KLD, KM, M;
      double               AJJ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHER, CLACGV, CSSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL, SQRT

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
         xerbla('CPBSTF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      KLD = max( 1, LDAB-1 );

      // Set the splitting point m.

      M = ( N+KD ) / 2;

      if ( UPPER ) {

         // Factorize A(m+1:n,m+1:n) as L**H*L, and update A(1:m,1:m).

         for (J = N; J >= M + 1; J--) { // 10

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = double( AB( KD+1, J ) );
            if ( AJJ <= ZERO ) {
               AB[KD+1][J] = AJJ;
               GO TO 50;
            }
            AJJ = sqrt( AJJ );
            AB[KD+1][J] = AJJ;
            KM = min( J-1, KD );

            // Compute elements j-km:j-1 of the j-th column and update the
            // the leading submatrix within the band.

            csscal(KM, ONE / AJJ, AB( KD+1-KM, J ), 1 );
            cher('Upper', KM, -ONE, AB( KD+1-KM, J ), 1, AB( KD+1, J-KM ), KLD );
         } // 10

         // Factorize the updated submatrix A(1:m,1:m) as U**H*U.

         for (J = 1; J <= M; J++) { // 20

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = double( AB( KD+1, J ) );
            if ( AJJ <= ZERO ) {
               AB[KD+1][J] = AJJ;
               GO TO 50;
            }
            AJJ = sqrt( AJJ );
            AB[KD+1][J] = AJJ;
            KM = min( KD, M-J );

            // Compute elements j+1:j+km of the j-th row and update the
            // trailing submatrix within the band.

            if ( KM > 0 ) {
               csscal(KM, ONE / AJJ, AB( KD, J+1 ), KLD );
               clacgv(KM, AB( KD, J+1 ), KLD );
               cher('Upper', KM, -ONE, AB( KD, J+1 ), KLD, AB( KD+1, J+1 ), KLD );
               clacgv(KM, AB( KD, J+1 ), KLD );
            }
         } // 20
      } else {

         // Factorize A(m+1:n,m+1:n) as L**H*L, and update A(1:m,1:m).

         for (J = N; J >= M + 1; J--) { // 30

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = double( AB( 1, J ) );
            if ( AJJ <= ZERO ) {
               AB[1][J] = AJJ;
               GO TO 50;
            }
            AJJ = sqrt( AJJ );
            AB[1][J] = AJJ;
            KM = min( J-1, KD );

            // Compute elements j-km:j-1 of the j-th row and update the
            // trailing submatrix within the band.

            csscal(KM, ONE / AJJ, AB( KM+1, J-KM ), KLD );
            clacgv(KM, AB( KM+1, J-KM ), KLD );
            cher('Lower', KM, -ONE, AB( KM+1, J-KM ), KLD, AB( 1, J-KM ), KLD );
            clacgv(KM, AB( KM+1, J-KM ), KLD );
         } // 30

         // Factorize the updated submatrix A(1:m,1:m) as U**H*U.

         for (J = 1; J <= M; J++) { // 40

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = double( AB( 1, J ) );
            if ( AJJ <= ZERO ) {
               AB[1][J] = AJJ;
               GO TO 50;
            }
            AJJ = sqrt( AJJ );
            AB[1][J] = AJJ;
            KM = min( KD, M-J );

            // Compute elements j+1:j+km of the j-th column and update the
            // trailing submatrix within the band.

            if ( KM > 0 ) {
               csscal(KM, ONE / AJJ, AB( 2, J ), 1 );
               cher('Lower', KM, -ONE, AB( 2, J ), 1, AB( 1, J+1 ), KLD );
            }
         } // 40
      }
      return;

      } // 50
      INFO = J;
      }
