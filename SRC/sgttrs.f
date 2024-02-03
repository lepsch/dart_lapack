      SUBROUTINE SGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               NOTRAN;
      int                ITRANS, J, JB, NB;
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGTTS2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0;
      NOTRAN = ( TRANS == 'N' || TRANS == 'n' );
      if ( !NOTRAN && !( TRANS == 'T' || TRANS == 't' ) && !( TRANS == 'C' || TRANS == 'c' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDB < MAX( N, 1 ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('SGTTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) RETURN;

      // Decode TRANS

      if ( NOTRAN ) {
         ITRANS = 0;
      } else {
         ITRANS = 1;
      }

      // Determine the number of right-hand sides to solve at a time.

      if ( NRHS == 1 ) {
         NB = 1;
      } else {
         NB = MAX( 1, ILAENV( 1, 'SGTTRS', TRANS, N, NRHS, -1, -1 ) );
      }

      if ( NB >= NRHS ) {
         sgtts2(ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB );
      } else {
         DO 10 J = 1, NRHS, NB;
            JB = MIN( NRHS-J+1, NB );
            sgtts2(ITRANS, N, JB, DL, D, DU, DU2, IPIV, B( 1, J ), LDB );
         } // 10
      }

      // End of SGTTRS

      }
