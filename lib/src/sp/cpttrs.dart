      void cpttrs(UPLO, N, NRHS, D, E, B, LDB, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               D( * );
      Complex            B( LDB, * ), E( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER;
      int                IUPLO, J, JB, NB;
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CPTTS2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0;
      UPPER = ( UPLO == 'U' || UPLO == 'u' );
      if ( !UPPER && !( UPLO == 'L' || UPLO == 'l' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }
      if ( INFO != 0 ) {
         xerbla('CPTTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      // Determine the number of right-hand sides to solve at a time.

      if ( NRHS == 1 ) {
         NB = 1;
      } else {
         NB = max( 1, ILAENV( 1, 'CPTTRS', UPLO, N, NRHS, -1, -1 ) );
      }

      // Decode UPLO

      if ( UPPER ) {
         IUPLO = 1;
      } else {
         IUPLO = 0;
      }

      if ( NB >= NRHS ) {
         cptts2(IUPLO, N, NRHS, D, E, B, LDB );
      } else {
         for (J = 1; NB < 0 ? J >= NRHS : J <= NRHS; J += NB) { // 10
            JB = min( NRHS-J+1, NB );
            cptts2(IUPLO, N, JB, D, E, B( 1, J ), LDB );
         } // 10
      }

      return;
      }
