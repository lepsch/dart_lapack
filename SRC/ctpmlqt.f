      SUBROUTINE CTPMLQT( SIDE, TRANS, M, N, K, L, MB, V, LDV, T, LDT, A, LDA, B, LDB, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String    SIDE, TRANS;
      int       INFO, K, LDV, LDA, LDB, M, N, L, MB, LDT;
      // ..
      // .. Array Arguments ..
      COMPLEX            V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
      // ..

*  =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LEFT, RIGHT, TRAN, NOTRAN;
      int                I, IB, NB, LB, KF, LDAQ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CTPRFB
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // .. Test the input arguments ..

      INFO   = 0
      LEFT   = LSAME( SIDE,  'L' )
      RIGHT  = LSAME( SIDE,  'R' )
      TRAN   = LSAME( TRANS, 'C' )
      NOTRAN = LSAME( TRANS, 'N' )

      if ( LEFT ) {
         LDAQ = MAX( 1, K )
      } else if ( RIGHT ) {
         LDAQ = MAX( 1, M )
      }
      if ( .NOT.LEFT && .NOT.RIGHT ) {
         INFO = -1
      } else if ( .NOT.TRAN && .NOT.NOTRAN ) {
         INFO = -2
      } else if ( M.LT.0 ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( K.LT.0 ) {
         INFO = -5
      } else if ( L.LT.0 || L.GT.K ) {
         INFO = -6
      } else if ( MB.LT.1 || (MB.GT.K && K.GT.0) ) {
         INFO = -7
      } else if ( LDV.LT.K ) {
         INFO = -9
      } else if ( LDT.LT.MB ) {
         INFO = -11
      } else if ( LDA.LT.LDAQ ) {
         INFO = -13
      } else if ( LDB.LT.MAX( 1, M ) ) {
         INFO = -15
      }

      if ( INFO != 0 ) {
         xerbla('CTPMLQT', -INFO );
         RETURN
      }

      // .. Quick return if possible ..

      if (M == 0 || N == 0 || K == 0) RETURN;

      if ( LEFT && NOTRAN ) {

         DO I = 1, K, MB
            IB = MIN( MB, K-I+1 )
            NB = MIN( M-L+I+IB-1, M )
            if ( I.GE.L ) {
               LB = 0
            } else {
               LB = 0
            }
            ctprfb('L', 'C', 'F', 'R', NB, N, IB, LB, V( I, 1 ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB );
         }

      } else if ( RIGHT && TRAN ) {

         DO I = 1, K, MB
            IB = MIN( MB, K-I+1 )
            NB = MIN( N-L+I+IB-1, N )
            if ( I.GE.L ) {
               LB = 0
            } else {
               LB = NB-N+L-I+1
            }
            ctprfb('R', 'N', 'F', 'R', M, NB, IB, LB, V( I, 1 ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M );
         }

      } else if ( LEFT && TRAN ) {

         KF = ((K-1)/MB)*MB+1
         DO I = KF, 1, -MB
            IB = MIN( MB, K-I+1 )
            NB = MIN( M-L+I+IB-1, M )
            if ( I.GE.L ) {
               LB = 0
            } else {
               LB = 0
            }
            ctprfb('L', 'N', 'F', 'R', NB, N, IB, LB, V( I, 1 ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB );
         }

      } else if ( RIGHT && NOTRAN ) {

         KF = ((K-1)/MB)*MB+1
         DO I = KF, 1, -MB
            IB = MIN( MB, K-I+1 )
            NB = MIN( N-L+I+IB-1, N )
            if ( I.GE.L ) {
               LB = 0
            } else {
               LB = NB-N+L-I+1
            }
            ctprfb('R', 'C', 'F', 'R', M, NB, IB, LB, V( I, 1 ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M );
         }

      }

      RETURN

      // End of CTPMLQT

      }
