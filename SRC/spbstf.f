      SUBROUTINE SPBSTF( UPLO, N, KD, AB, LDAB, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, N;
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, KLD, KM, M;
      REAL               AJJ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSYR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
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
         xerbla('SPBSTF', -INFO );
         RETURN;
      }

      // Quick return if possible

      if (N == 0) RETURN;

      KLD = MAX( 1, LDAB-1 );

      // Set the splitting point m.

      M = ( N+KD ) / 2;

      if ( UPPER ) {

         // Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m).

         DO 10 J = N, M + 1, -1;

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = AB( KD+1, J );
            if (AJJ <= ZERO) GO TO 50;
            AJJ = SQRT( AJJ );
            AB( KD+1, J ) = AJJ;
            KM = MIN( J-1, KD );

            // Compute elements j-km:j-1 of the j-th column and update the
            // the leading submatrix within the band.

            sscal(KM, ONE / AJJ, AB( KD+1-KM, J ), 1 );
            ssyr('Upper', KM, -ONE, AB( KD+1-KM, J ), 1, AB( KD+1, J-KM ), KLD );
         } // 10

         // Factorize the updated submatrix A(1:m,1:m) as U**T*U.

         for (J = 1; J <= M; J++) { // 20

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = AB( KD+1, J );
            if (AJJ <= ZERO) GO TO 50;
            AJJ = SQRT( AJJ );
            AB( KD+1, J ) = AJJ;
            KM = MIN( KD, M-J );

            // Compute elements j+1:j+km of the j-th row and update the
            // trailing submatrix within the band.

            if ( KM > 0 ) {
               sscal(KM, ONE / AJJ, AB( KD, J+1 ), KLD );
               ssyr('Upper', KM, -ONE, AB( KD, J+1 ), KLD, AB( KD+1, J+1 ), KLD );
            }
         } // 20
      } else {

         // Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m).

         DO 30 J = N, M + 1, -1;

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = AB( 1, J );
            if (AJJ <= ZERO) GO TO 50;
            AJJ = SQRT( AJJ );
            AB( 1, J ) = AJJ;
            KM = MIN( J-1, KD );

            // Compute elements j-km:j-1 of the j-th row and update the
            // trailing submatrix within the band.

            sscal(KM, ONE / AJJ, AB( KM+1, J-KM ), KLD );
            ssyr('Lower', KM, -ONE, AB( KM+1, J-KM ), KLD, AB( 1, J-KM ), KLD );
         } // 30

         // Factorize the updated submatrix A(1:m,1:m) as U**T*U.

         for (J = 1; J <= M; J++) { // 40

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = AB( 1, J );
            if (AJJ <= ZERO) GO TO 50;
            AJJ = SQRT( AJJ );
            AB( 1, J ) = AJJ;
            KM = MIN( KD, M-J );

            // Compute elements j+1:j+km of the j-th column and update the
            // trailing submatrix within the band.

            if ( KM > 0 ) {
               sscal(KM, ONE / AJJ, AB( 2, J ), 1 );
               ssyr('Lower', KM, -ONE, AB( 2, J ), 1, AB( 1, J+1 ), KLD );
            }
         } // 40
      }
      RETURN;

      } // 50
      INFO = J;
      RETURN;

      // End of SPBSTF

      }
