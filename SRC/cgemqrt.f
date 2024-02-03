      SUBROUTINE CGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT, C, LDC, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String    SIDE, TRANS;
      int       INFO, K, LDV, LDC, M, N, NB, LDT;
      // ..
      // .. Array Arguments ..
      COMPLEX   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
      // ..

*  =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LEFT, RIGHT, TRAN, NOTRAN;
      int                I, IB, LDWORK, KF, Q;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CLARFB
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
         LDWORK = MAX( 1, N )
         Q = M
      } else if ( RIGHT ) {
         LDWORK = MAX( 1, M )
         Q = N
      }
      if ( .NOT.LEFT && .NOT.RIGHT ) {
         INFO = -1
      } else if ( .NOT.TRAN && .NOT.NOTRAN ) {
         INFO = -2
      } else if ( M.LT.0 ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( K.LT.0 || K.GT.Q ) {
         INFO = -5
      } else if ( NB.LT.1 || (NB.GT.K && K.GT.0)) {
         INFO = -6
      } else if ( LDV.LT.MAX( 1, Q ) ) {
         INFO = -8
      } else if ( LDT.LT.NB ) {
         INFO = -10
      } else if ( LDC.LT.MAX( 1, M ) ) {
         INFO = -12
      }

      if ( INFO != 0 ) {
         xerbla('CGEMQRT', -INFO );
         RETURN
      }

      // .. Quick return if possible ..

      if (M == 0 || N == 0 || K == 0) RETURN;

      if ( LEFT && TRAN ) {

         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            clarfb('L', 'C', 'F', 'C', M-I+1, N, IB, V( I, I ), LDV, T( 1, I ), LDT, C( I, 1 ), LDC, WORK, LDWORK );
         }

      } else if ( RIGHT && NOTRAN ) {

         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            clarfb('R', 'N', 'F', 'C', M, N-I+1, IB, V( I, I ), LDV, T( 1, I ), LDT, C( 1, I ), LDC, WORK, LDWORK );
         }

      } else if ( LEFT && NOTRAN ) {

         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            clarfb('L', 'N', 'F', 'C', M-I+1, N, IB, V( I, I ), LDV, T( 1, I ), LDT, C( I, 1 ), LDC, WORK, LDWORK );
         }

      } else if ( RIGHT && TRAN ) {

         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            clarfb('R', 'C', 'F', 'C', M, N-I+1, IB, V( I, I ), LDV, T( 1, I ), LDT, C( 1, I ), LDC, WORK, LDWORK );
         }

      }

      RETURN

      // End of CGEMQRT

      }
