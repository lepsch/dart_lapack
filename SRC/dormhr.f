      SUBROUTINE DORMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                IHI, ILO, INFO, LDA, LDC, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LEFT, LQUERY;
      int                I1, I2, IINFO, LWKOPT, MI, NB, NH, NI, NQ, NW;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL DORMQR, XERBLA
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
      if ( .NOT.LEFT && .NOT.LSAME( SIDE, 'R' ) ) {
         INFO = -1
      } else if ( .NOT.LSAME( TRANS, 'N' ) && .NOT.LSAME( TRANS, 'T' ) ) {
         INFO = -2
      } else if ( M < 0 ) {
         INFO = -3
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( ILO < 1 || ILO > MAX( 1, NQ ) ) {
         INFO = -5
      } else if ( IHI < MIN( ILO, NQ ) || IHI > NQ ) {
         INFO = -6
      } else if ( LDA < MAX( 1, NQ ) ) {
         INFO = -8
      } else if ( LDC < MAX( 1, M ) ) {
         INFO = -11
      } else if ( LWORK < NW && .NOT.LQUERY ) {
         INFO = -13
      }

      if ( INFO == 0 ) {
         if ( LEFT ) {
            NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, NH, N, NH, -1 )
         } else {
            NB = ILAENV( 1, 'DORMQR', SIDE // TRANS, M, NH, NH, -1 )
         }
         LWKOPT = NW*NB
         WORK( 1 ) = LWKOPT
      }

      if ( INFO != 0 ) {
         xerbla('DORMHR', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( M == 0 || N == 0 || NH == 0 ) {
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

      dormqr(SIDE, TRANS, MI, NI, NH, A( ILO+1, ILO ), LDA, TAU( ILO ), C( I1, I2 ), LDC, WORK, LWORK, IINFO );

      WORK( 1 ) = LWKOPT
      RETURN

      // End of DORMHR

      }
