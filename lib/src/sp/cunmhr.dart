      void cunmhr(SIDE, TRANS, M, N, ILO, IHI, final Matrix<double> A, final int LDA, TAU, final Matrix<double> C, final int LDC, final Array<double> WORK, final int LWORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE, TRANS;
      int                IHI, ILO, INFO, LDA, LDC, LWORK, M, N;
      Complex            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LEFT, LQUERY;
      int                I1, I2, IINFO, LWKOPT, MI, NB, NH, NI, NQ, NW;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CUNMQR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input arguments

      INFO = 0;
      NH = IHI - ILO;
      LEFT = lsame( SIDE, 'L' );
      LQUERY = ( LWORK == -1 );

      // NQ is the order of Q and NW is the minimum dimension of WORK

      if ( LEFT ) {
         NQ = M;
         NW = max( 1, N );
      } else {
         NQ = N;
         NW = max( 1, M );
      }
      if ( !LEFT && !lsame( SIDE, 'R' ) ) {
         INFO = -1;
      } else if ( !lsame( TRANS, 'N' ) && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( M < 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( ILO < 1 || ILO > max( 1, NQ ) ) {
         INFO = -5;
      } else if ( IHI < min( ILO, NQ ) || IHI > NQ ) {
         INFO = -6;
      } else if ( LDA < max( 1, NQ ) ) {
         INFO = -8;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -11;
      } else if ( LWORK < NW && !LQUERY ) {
         INFO = -13;
      }

      if ( INFO == 0 ) {
         if ( LEFT ) {
            NB = ilaenv( 1, 'CUNMQR', SIDE + TRANS, NH, N, NH, -1 );
         } else {
            NB = ilaenv( 1, 'CUNMQR', SIDE + TRANS, M, NH, NH, -1 );
         }
         LWKOPT = NW*NB;
         WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }

      if ( INFO != 0 ) {
         xerbla('CUNMHR', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( M == 0 || N == 0 || NH == 0 ) {
         WORK[1] = 1;
         return;
      }

      if ( LEFT ) {
         MI = NH;
         NI = N;
         I1 = ILO + 1;
         I2 = 1;
      } else {
         MI = M;
         NI = NH;
         I1 = 1;
         I2 = ILO + 1;
      }

      cunmqr(SIDE, TRANS, MI, NI, NH, A( ILO+1, ILO ), LDA, TAU( ILO ), C( I1, I2 ), LDC, WORK, LWORK, IINFO );

      WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }
