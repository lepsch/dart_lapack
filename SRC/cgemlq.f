      SUBROUTINE CGEMLQ( SIDE, TRANS, M, N, K, A, LDA, T, TSIZE, C, LDC, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                INFO, LDA, M, N, K, TSIZE, LWORK, LDC;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), T( * ), C( LDC, * ), WORK( * )
      // ..

* =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LEFT, RIGHT, TRAN, NOTRAN, LQUERY;
      int                MB, NB, LW, NBLCKS, MN, MINMNK, LWMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLAMSWLQ, CGEMLQT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN, MOD
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      LQUERY  = ( LWORK == -1 )
      NOTRAN  = LSAME( TRANS, 'N' )
      TRAN    = LSAME( TRANS, 'C' )
      LEFT    = LSAME( SIDE, 'L' )
      RIGHT   = LSAME( SIDE, 'R' )

      MB = INT( T( 2 ) )
      NB = INT( T( 3 ) )
      if ( LEFT ) {
        LW = N * MB
        MN = M
      } else {
        LW = M * MB
        MN = N
      }

      MINMNK = MIN( M, N, K )
      if ( MINMNK == 0 ) {
         LWMIN = 1
      } else {
         LWMIN = MAX( 1, LW )
      }

      if ( ( NB > K ) && ( MN > K ) ) {
        if ( MOD( MN - K, NB - K ) == 0 ) {
          NBLCKS = ( MN - K ) / ( NB - K )
        } else {
          NBLCKS = ( MN - K ) / ( NB - K ) + 1
        }
      } else {
        NBLCKS = 1
      }

      INFO = 0
      if ( .NOT.LEFT && .NOT.RIGHT ) {
        INFO = -1
      } else if ( .NOT.TRAN && .NOT.NOTRAN ) {
        INFO = -2
      } else if ( M < 0 ) {
        INFO = -3
      } else if ( N < 0 ) {
        INFO = -4
      } else if ( K < 0 || K > MN ) {
        INFO = -5
      } else if ( LDA < MAX( 1, K ) ) {
        INFO = -7
      } else if ( TSIZE < 5 ) {
        INFO = -9
      } else if ( LDC < MAX( 1, M ) ) {
        INFO = -11
      } else if ( ( LWORK < LWMIN ) && ( .NOT.LQUERY ) ) {
        INFO = -13
      }

      if ( INFO == 0 ) {
        WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      }

      if ( INFO != 0 ) {
        xerbla('CGEMLQ', -INFO );
        RETURN
      } else if ( LQUERY ) {
        RETURN
      }

      // Quick return if possible

      if ( MINMNK == 0 ) {
        RETURN
      }

      IF( ( LEFT && M.LE.K ) || ( RIGHT && N.LE.K ) || ( NB.LE.K ) || ( NB.GE.MAX( M, N, K ) ) ) THEN         CALL CGEMLQT( SIDE, TRANS, M, N, K, MB, A, LDA, T( 6 ), MB, C, LDC, WORK, INFO )
      } else {
        clamswlq(SIDE, TRANS, M, N, K, MB, NB, A, LDA, T( 6 ), MB, C, LDC, WORK, LWORK, INFO );
      }

      WORK( 1 ) = SROUNDUP_LWORK( LWMIN )

      RETURN

      // End of CGEMLQ

      }
