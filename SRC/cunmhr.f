      SUBROUTINE CUNMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                IHI, ILO, INFO, LDA, LDC, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LEFT, LQUERY;
      int                I1, I2, IINFO, LWKOPT, MI, NB, NH, NI, NQ, NW;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CUNMQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      NH = IHI - ILO
      LEFT = LSAME( SIDE, 'L' )
      LQUERY = ( LWORK == -1 )

      // NQ is the order of Q and NW is the minimum dimension of WORK

      if ( LEFT ) {
         NQ = M
         NW = MAX( 1, N )
      } else {
         NQ = N
         NW = MAX( 1, M )
      }
      if ( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( M.LT.0 ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( ILO.LT.1 .OR. ILO.GT.MAX( 1, NQ ) ) {
         INFO = -5
      } else if ( IHI.LT.MIN( ILO, NQ ) .OR. IHI.GT.NQ ) {
         INFO = -6
      } else if ( LDA.LT.MAX( 1, NQ ) ) {
         INFO = -8
      } else if ( LDC.LT.MAX( 1, M ) ) {
         INFO = -11
      } else if ( LWORK.LT.NW .AND. .NOT.LQUERY ) {
         INFO = -13
      }

      if ( INFO == 0 ) {
         if ( LEFT ) {
            NB = ILAENV( 1, 'CUNMQR', SIDE // TRANS, NH, N, NH, -1 )
         } else {
            NB = ILAENV( 1, 'CUNMQR', SIDE // TRANS, M, NH, NH, -1 )
         }
         LWKOPT = NW*NB
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      }

      if ( INFO.NE.0 ) {
         xerbla('CUNMHR', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M == 0 .OR. N == 0 .OR. NH == 0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      if ( LEFT ) {
         MI = NH
         NI = N
         I1 = ILO + 1
         I2 = 1
      } else {
         MI = M
         NI = NH
         I1 = 1
         I2 = ILO + 1
      }

      cunmqr(SIDE, TRANS, MI, NI, NH, A( ILO+1, ILO ), LDA, TAU( ILO ), C( I1, I2 ), LDC, WORK, LWORK, IINFO );

      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      RETURN

      // End of CUNMHR

      }
