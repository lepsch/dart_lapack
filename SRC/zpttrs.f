      SUBROUTINE ZPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      double             D( * );
      COMPLEX*16         B( LDB, * ), E( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               UPPER;
      int                IUPLO, J, JB, NB;
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZPTTS2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0
      UPPER = ( UPLO == 'U' || UPLO == 'u' )
      if ( .NOT.UPPER && .NOT.( UPLO == 'L' || UPLO == 'l' ) ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( NRHS < 0 ) {
         INFO = -3
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO != 0 ) {
         xerbla('ZPTTRS', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) RETURN;

      // Determine the number of right-hand sides to solve at a time.

      if ( NRHS == 1 ) {
         NB = 1
      } else {
         NB = MAX( 1, ILAENV( 1, 'ZPTTRS', UPLO, N, NRHS, -1, -1 ) )
      }

      // Decode UPLO

      if ( UPPER ) {
         IUPLO = 1
      } else {
         IUPLO = 0
      }

      if ( NB.GE.NRHS ) {
         zptts2(IUPLO, N, NRHS, D, E, B, LDB );
      } else {
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            zptts2(IUPLO, N, JB, D, E, B( 1, J ), LDB );
         } // 10
      }

      RETURN

      // End of ZPTTRS

      }
