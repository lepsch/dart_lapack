      SUBROUTINE STPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT, A, LDA, B, LDB, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String    SIDE, TRANS;
      int       INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT;
      // ..
      // .. Array Arguments ..
      REAL   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
      // ..

*  =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LEFT, RIGHT, TRAN, NOTRAN;
      int                I, IB, MB, LB, KF, LDAQ, LDVQ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL STPRFB, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // .. Test the input arguments ..

      INFO   = 0
      LEFT   = LSAME( SIDE,  'L' )
      RIGHT  = LSAME( SIDE,  'R' )
      TRAN   = LSAME( TRANS, 'T' )
      NOTRAN = LSAME( TRANS, 'N' )

      if ( LEFT ) {
         LDVQ = MAX( 1, M )
         LDAQ = MAX( 1, K )
      } else if ( RIGHT ) {
         LDVQ = MAX( 1, N )
         LDAQ = MAX( 1, M )
      }
      if ( .NOT.LEFT && .NOT.RIGHT ) {
         INFO = -1
      } else if ( .NOT.TRAN && .NOT.NOTRAN ) {
         INFO = -2
      } else if ( M < 0 ) {
         INFO = -3
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( K < 0 ) {
         INFO = -5
      } else if ( L < 0 || L > K ) {
         INFO = -6
      } else if ( NB < 1 || (NB > K && K > 0) ) {
         INFO = -7
      } else if ( LDV < LDVQ ) {
         INFO = -9
      } else if ( LDT < NB ) {
         INFO = -11
      } else if ( LDA < LDAQ ) {
         INFO = -13
      } else if ( LDB < MAX( 1, M ) ) {
         INFO = -15
      }

      if ( INFO != 0 ) {
         xerbla('STPMQRT', -INFO );
         RETURN
      }

      // .. Quick return if possible ..

      if (M == 0 || N == 0 || K == 0) RETURN;

      if ( LEFT && TRAN ) {

         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( M-L+I+IB-1, M )
            if ( I >= L ) {
               LB = 0
            } else {
               LB = MB-M+L-I+1
            }
            stprfb('L', 'T', 'F', 'C', MB, N, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB );
         }

      } else if ( RIGHT && NOTRAN ) {

         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( N-L+I+IB-1, N )
            if ( I >= L ) {
               LB = 0
            } else {
               LB = MB-N+L-I+1
            }
            stprfb('R', 'N', 'F', 'C', M, MB, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M );
         }

      } else if ( LEFT && NOTRAN ) {

         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( M-L+I+IB-1, M )
            if ( I >= L ) {
               LB = 0
            } else {
               LB = MB-M+L-I+1
            }
            stprfb('L', 'N', 'F', 'C', MB, N, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB );
         }

      } else if ( RIGHT && TRAN ) {

         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( N-L+I+IB-1, N )
            if ( I >= L ) {
               LB = 0
            } else {
               LB = MB-N+L-I+1
            }
            stprfb('R', 'T', 'F', 'C', M, MB, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M );
         }

      }

      RETURN

      // End of STPMQRT

      }
