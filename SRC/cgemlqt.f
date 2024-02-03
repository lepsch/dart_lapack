      SUBROUTINE CGEMLQT( SIDE, TRANS, M, N, K, MB, V, LDV, T, LDT, C, LDC, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String    SIDE, TRANS;
      int       INFO, K, LDV, LDC, M, N, MB, LDT;
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
      if ( !LEFT && !RIGHT ) {
         INFO = -1
      } else if ( !TRAN && !NOTRAN ) {
         INFO = -2
      } else if ( M < 0 ) {
         INFO = -3
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( K < 0 || K > Q ) {
         INFO = -5
      } else if ( MB < 1 || (MB > K && K > 0)) {
         INFO = -6
      } else if ( LDV < MAX( 1, K ) ) {
          INFO = -8
      } else if ( LDT < MB ) {
         INFO = -10
      } else if ( LDC < MAX( 1, M ) ) {
         INFO = -12
      }

      if ( INFO != 0 ) {
         xerbla('CGEMLQT', -INFO );
         RETURN
      }

      // .. Quick return if possible ..

      if (M == 0 || N == 0 || K == 0) RETURN;

      if ( LEFT && NOTRAN ) {

         DO I = 1, K, MB
            IB = MIN( MB, K-I+1 )
            clarfb('L', 'C', 'F', 'R', M-I+1, N, IB, V( I, I ), LDV, T( 1, I ), LDT, C( I, 1 ), LDC, WORK, LDWORK );
         }

      } else if ( RIGHT && TRAN ) {

         DO I = 1, K, MB
            IB = MIN( MB, K-I+1 )
            clarfb('R', 'N', 'F', 'R', M, N-I+1, IB, V( I, I ), LDV, T( 1, I ), LDT, C( 1, I ), LDC, WORK, LDWORK );
         }

      } else if ( LEFT && TRAN ) {

         KF = ((K-1)/MB)*MB+1
         DO I = KF, 1, -MB
            IB = MIN( MB, K-I+1 )
            clarfb('L', 'N', 'F', 'R', M-I+1, N, IB, V( I, I ), LDV, T( 1, I ), LDT, C( I, 1 ), LDC, WORK, LDWORK );
         }

      } else if ( RIGHT && NOTRAN ) {

         KF = ((K-1)/MB)*MB+1
         DO I = KF, 1, -MB
            IB = MIN( MB, K-I+1 )
            clarfb('R', 'C', 'F', 'R', M, N-I+1, IB, V( I, I ), LDV, T( 1, I ), LDT, C( 1, I ), LDC, WORK, LDWORK );
         }

      }

      RETURN

      // End of CGEMLQT

      }
