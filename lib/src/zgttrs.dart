      void zgttrs(final int TRANS, final int N, final int NRHS, final int DL, final int D, final int DU, final int DU2, final Array<int> IPIV_, final Matrix<double> B_, final int LDB, final Box<int> INFO,) {
  final IPIV = IPIV_.dim();
  final B = B_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS;
      int                INFO, LDB, N, NRHS;
      int                IPIV( * );
      Complex         B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               NOTRAN;
      int                ITRANS, J, JB, NB;
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGTTS2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

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
         xerbla('ZGTTRS', -INFO );
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
         NB = max( 1, ilaenv( 1, 'ZGTTRS', TRANS, N, NRHS, -1, -1 ) );
      }

      if ( NB >= NRHS ) {
         zgtts2(ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB );
      } else {
         for (J = 1; NB < 0 ? J >= NRHS : J <= NRHS; J += NB) { // 10
            JB = min( NRHS-J+1, NB );
            zgtts2(ITRANS, N, JB, DL, D, DU, DU2, IPIV, B( 1, J ), LDB );
         } // 10
      }
      }
