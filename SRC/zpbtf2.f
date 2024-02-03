      SUBROUTINE ZPBTF2( UPLO, N, KD, AB, LDAB, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         AB( LDAB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, KLD, KN;
      double             AJJ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZHER, ZLACGV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KD < 0 ) {
         INFO = -3;
      } else if ( LDAB < KD+1 ) {
         INFO = -5;
      }
      if ( INFO != 0 ) {
         xerbla('ZPBTF2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      KLD = max( 1, LDAB-1 );

      if ( UPPER ) {

         // Compute the Cholesky factorization A = U**H * U.

         for (J = 1; J <= N; J++) { // 10

            // Compute U(J,J) and test for non-positive-definiteness.

            AJJ = DBLE( AB( KD+1, J ) );
            if ( AJJ <= ZERO ) {
               AB( KD+1, J ) = AJJ;
               GO TO 30;
            }
            AJJ = sqrt( AJJ );
            AB( KD+1, J ) = AJJ;

            // Compute elements J+1:J+KN of row J and update the
            // trailing submatrix within the band.

            KN = min( KD, N-J );
            if ( KN > 0 ) {
               zdscal(KN, ONE / AJJ, AB( KD, J+1 ), KLD );
               zlacgv(KN, AB( KD, J+1 ), KLD );
               zher('Upper', KN, -ONE, AB( KD, J+1 ), KLD, AB( KD+1, J+1 ), KLD );
               zlacgv(KN, AB( KD, J+1 ), KLD );
            }
         } // 10
      } else {

         // Compute the Cholesky factorization A = L*L**H.

         for (J = 1; J <= N; J++) { // 20

            // Compute L(J,J) and test for non-positive-definiteness.

            AJJ = DBLE( AB( 1, J ) );
            if ( AJJ <= ZERO ) {
               AB( 1, J ) = AJJ;
               GO TO 30;
            }
            AJJ = sqrt( AJJ );
            AB( 1, J ) = AJJ;

            // Compute elements J+1:J+KN of column J and update the
            // trailing submatrix within the band.

            KN = min( KD, N-J );
            if ( KN > 0 ) {
               zdscal(KN, ONE / AJJ, AB( 2, J ), 1 );
               zher('Lower', KN, -ONE, AB( 2, J ), 1, AB( 1, J+1 ), KLD );
            }
         } // 20
      }
      return;

      } // 30
      INFO = J;
      return;

      // End of ZPBTF2

      }
