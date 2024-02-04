      void sgemqr(SIDE, TRANS, M, N, K, A, LDA, T, TSIZE, C, LDC, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                INFO, LDA, M, N, K, TSIZE, LWORK, LDC;
      // ..
      // .. Array Arguments ..
      double               A( LDA, * ), T( * ), C( LDC, * ), WORK( * );
      // ..

// =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LEFT, RIGHT, TRAN, NOTRAN, LQUERY;
      int                MB, NB, LW, NBLCKS, MN, MINMNK, LWMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      double               SROUNDUP_LWORK;
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMQRT, SLAMTSQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN, MOD
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      LQUERY  = ( LWORK == -1 );
      NOTRAN  = lsame( TRANS, 'N' );
      TRAN    = lsame( TRANS, 'T' );
      LEFT    = lsame( SIDE, 'L' );
      RIGHT   = lsame( SIDE, 'R' );

      MB = INT( T( 2 ) );
      NB = INT( T( 3 ) );
      if ( LEFT ) {
        LW = N * NB;
        MN = M;
      } else {
        LW = MB * NB;
        MN = N;
      }

      MINMNK = min( M, N, K );
      if ( MINMNK == 0 ) {
         LWMIN = 1;
      } else {
         LWMIN = max( 1, LW );
      }

      if ( ( MB > K ) && ( MN > K ) ) {
        if ( (MN - K % MB - K) == 0 ) {
          NBLCKS = ( MN - K ) / ( MB - K );
        } else {
          NBLCKS = ( MN - K ) / ( MB - K ) + 1;
        }
      } else {
        NBLCKS = 1;
      }

      INFO = 0;
      if ( !LEFT && !RIGHT ) {
        INFO = -1;
      } else if ( !TRAN && !NOTRAN ) {
        INFO = -2;
      } else if ( M < 0 ) {
        INFO = -3;
      } else if ( N < 0 ) {
        INFO = -4;
      } else if ( K < 0 || K > MN ) {
        INFO = -5;
      } else if ( LDA < max( 1, MN ) ) {
        INFO = -7;
      } else if ( TSIZE < 5 ) {
        INFO = -9;
      } else if ( LDC < max( 1, M ) ) {
        INFO = -11;
      } else if ( LWORK < LWMIN && !LQUERY ) {
        INFO = -13;
      }

      if ( INFO == 0 ) {
        WORK[1] = SROUNDUP_LWORK( LWMIN );
      }

      if ( INFO != 0 ) {
        xerbla('SGEMQR', -INFO );
        return;
      } else if ( LQUERY ) {
        return;
      }

      // Quick return if possible

      if ( MINMNK == 0 ) {
        return;
      }

      if( ( LEFT && M <= K ) || ( RIGHT && N <= K ) || ( MB <= K ) || ( MB >= max( M, N, K ) ) ) {
        CALL SGEMQRT( SIDE, TRANS, M, N, K, NB, A, LDA, T( 6 ), NB, C, LDC, WORK, INFO );
      } else {
        slamtsqr(SIDE, TRANS, M, N, K, MB, NB, A, LDA, T( 6 ), NB, C, LDC, WORK, LWORK, INFO );
      }

      WORK[1] = SROUNDUP_LWORK( LWMIN );

      return;
      }
