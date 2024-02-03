      SUBROUTINE SPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               AP( * ), B( LDB, * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               UPPER;
      int                I;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL STPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('SPPTRS', -INFO );
         RETURN;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) RETURN;

      if ( UPPER ) {

         // Solve A*X = B where A = U**T * U.

         for (I = 1; I <= NRHS; I++) { // 10

            // Solve U**T *X = B, overwriting B with X.

            stpsv('Upper', 'Transpose', 'Non-unit', N, AP, B( 1, I ), 1 );

            // Solve U*X = B, overwriting B with X.

            stpsv('Upper', 'No transpose', 'Non-unit', N, AP, B( 1, I ), 1 );
         } // 10
      } else {

         // Solve A*X = B where A = L * L**T.

         for (I = 1; I <= NRHS; I++) { // 20

            // Solve L*Y = B, overwriting B with X.

            stpsv('Lower', 'No transpose', 'Non-unit', N, AP, B( 1, I ), 1 );

            // Solve L**T *X = Y, overwriting B with X.

            stpsv('Lower', 'Transpose', 'Non-unit', N, AP, B( 1, I ), 1 );
         } // 20
      }

      RETURN;

      // End of SPPTRS

      }
