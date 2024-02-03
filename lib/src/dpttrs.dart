      void dpttrs(N, NRHS, D, E, B, LDB, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             B( LDB, * ), D( * ), E( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                J, JB, NB;
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DPTTS2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
      } else if ( NRHS < 0 ) {
         INFO = -2;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('DPTTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      // Determine the number of right-hand sides to solve at a time.

      if ( NRHS == 1 ) {
         NB = 1;
      } else {
         NB = max( 1, ILAENV( 1, 'DPTTRS', ' ', N, NRHS, -1, -1 ) );
      }

      if ( NB >= NRHS ) {
         dptts2(N, NRHS, D, E, B, LDB );
      } else {
         for (J = 1; NB < 0 ? J >= NRHS : J <= NRHS; J += NB) { // 10
            JB = min( NRHS-J+1, NB );
            dptts2(N, JB, D, E, B( 1, J ), LDB );
         } // 10
      }

      return;
      }
