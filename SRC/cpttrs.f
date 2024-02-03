      SUBROUTINE CPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               D( * )
      COMPLEX            B( LDB, * ), E( * )
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
      // EXTERNAL CPTTS2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments.

      INFO = 0
      UPPER = ( UPLO == 'U' .OR. UPLO == 'u' )
      if ( .NOT.UPPER && .NOT.( UPLO == 'L' .OR. UPLO == 'l' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      }
      if ( INFO != 0 ) {
         xerbla('CPTTRS', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0 .OR. NRHS == 0) RETURN;

      // Determine the number of right-hand sides to solve at a time.

      if ( NRHS == 1 ) {
         NB = 1
      } else {
         NB = MAX( 1, ILAENV( 1, 'CPTTRS', UPLO, N, NRHS, -1, -1 ) )
      }

      // Decode UPLO

      if ( UPPER ) {
         IUPLO = 1
      } else {
         IUPLO = 0
      }

      if ( NB.GE.NRHS ) {
         cptts2(IUPLO, N, NRHS, D, E, B, LDB );
      } else {
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            cptts2(IUPLO, N, JB, D, E, B( 1, J ), LDB );
         } // 10
      }

      RETURN

      // End of CPTTRS

      }
