      void cgttrs(TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                INFO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               NOTRAN;
      int                ITRANS, J, JB, NB;
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGTTS2, XERBLA
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
      } else if ( LDB < max( N, 1 ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('CGTTRS', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      // Decode TRANS

      if ( NOTRAN ) {
         ITRANS = 0;
      } else if ( TRANS == 'T' || TRANS == 't' ) {
         ITRANS = 1;
      } else {
         ITRANS = 2;
      }

      // Determine the number of right-hand sides to solve at a time.

      if ( NRHS == 1 ) {
         NB = 1;
      } else {
         NB = max( 1, ILAENV( 1, 'CGTTRS', TRANS, N, NRHS, -1, -1 ) );
      }

      if ( NB >= NRHS ) {
         cgtts2(ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB );
      } else {
         DO 10 J = 1, NRHS, NB;
            JB = min( NRHS-J+1, NB );
            cgtts2(ITRANS, N, JB, DL, D, DU, DU2, IPIV, B( 1, J ), LDB );
         } // 10
      }

      // End of CGTTRS

      }
