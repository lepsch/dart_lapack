      SUBROUTINE SGEMLQ( SIDE, TRANS, M, N, K, A, LDA, T, TSIZE, C, LDC, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                INFO, LDA, M, N, K, TSIZE, LWORK, LDC;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), T( * ), C( LDC, * ), WORK( * )
      // ..

* =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LEFT, RIGHT, TRAN, NOTRAN, LQUERY;
      int                MB, NB, LW, NBLCKS, MN, MINMNK, LWMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Functions ..
      REAL               SROUNDUP_LWORK
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAMSWLQ, SGEMLQT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN, MOD
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      LQUERY  = ( LWORK.EQ.-1 )
      NOTRAN  = LSAME( TRANS, 'N' )
      TRAN    = LSAME( TRANS, 'T' )
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
      if ( MINMNK.EQ.0 ) {
         LWMIN = 1
      } else {
         LWMIN = MAX( 1, LW )
      }

      if ( ( NB.GT.K ) .AND. ( MN.GT.K ) ) {
        if ( MOD( MN - K, NB - K ) .EQ. 0 ) {
          NBLCKS = ( MN - K ) / ( NB - K )
        } else {
          NBLCKS = ( MN - K ) / ( NB - K ) + 1
        }
      } else {
        NBLCKS = 1
      }

      INFO = 0
      if ( .NOT.LEFT .AND. .NOT.RIGHT ) {
        INFO = -1
      } else if ( .NOT.TRAN .AND. .NOT.NOTRAN ) {
        INFO = -2
      } else if ( M.LT.0 ) {
        INFO = -3
      } else if ( N.LT.0 ) {
        INFO = -4
      } else if ( K.LT.0 .OR. K.GT.MN ) {
        INFO = -5
      } else if ( LDA.LT.MAX( 1, K ) ) {
        INFO = -7
      } else if ( TSIZE.LT.5 ) {
        INFO = -9
      } else if ( LDC.LT.MAX( 1, M ) ) {
        INFO = -11
      } else if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
        INFO = -13
      }

      if ( INFO.EQ.0 ) {
        WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      }

      if ( INFO.NE.0 ) {
        CALL XERBLA( 'SGEMLQ', -INFO )
        RETURN
      } else if ( LQUERY ) {
        RETURN
      }

      // Quick return if possible

      if ( MINMNK.EQ.0 ) {
        RETURN
      }

      IF( ( LEFT .AND. M.LE.K ) .OR. ( RIGHT .AND. N.LE.K ) .OR. ( NB.LE.K ) .OR. ( NB.GE.MAX( M, N, K ) ) ) THEN         CALL SGEMLQT( SIDE, TRANS, M, N, K, MB, A, LDA, T( 6 ), MB, C, LDC, WORK, INFO )
      } else {
        CALL SLAMSWLQ( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T( 6 ), MB, C, LDC, WORK, LWORK, INFO )
      }

      WORK( 1 ) = SROUNDUP_LWORK( LWMIN )

      RETURN

      // End of SGEMLQ

      }
