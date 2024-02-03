      SUBROUTINE SPTTRS( N, NRHS, D, E, B, LDB, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               B( LDB, * ), D( * ), E( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                J, JB, NB;
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL SPTTS2, XERBLA
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
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('SPTTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      // Determine the number of right-hand sides to solve at a time.

      if ( NRHS == 1 ) {
         NB = 1;
      } else {
         NB = MAX( 1, ILAENV( 1, 'SPTTRS', ' ', N, NRHS, -1, -1 ) );
      }

      if ( NB >= NRHS ) {
         sptts2(N, NRHS, D, E, B, LDB );
      } else {
         DO 10 J = 1, NRHS, NB;
            JB = MIN( NRHS-J+1, NB );
            sptts2(N, JB, D, E, B( 1, J ), LDB );
         } // 10
      }

      return;

      // End of SPTTRS

      }
