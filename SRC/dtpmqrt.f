      SUBROUTINE DTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT, A, LDA, B, LDB, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String    SIDE, TRANS;
      int       INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT;
      // ..
      // .. Array Arguments ..
      double             V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * );
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
      // EXTERNAL DTPRFB, XERBLA
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
      if ( .NOT.LEFT .AND. .NOT.RIGHT ) {
         INFO = -1
      } else if ( .NOT.TRAN .AND. .NOT.NOTRAN ) {
         INFO = -2
      } else if ( M.LT.0 ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( K.LT.0 ) {
         INFO = -5
      } else if ( L.LT.0 .OR. L.GT.K ) {
         INFO = -6
      } else if ( NB.LT.1 .OR. (NB.GT.K .AND. K.GT.0) ) {
         INFO = -7
      } else if ( LDV.LT.LDVQ ) {
         INFO = -9
      } else if ( LDT.LT.NB ) {
         INFO = -11
      } else if ( LDA.LT.LDAQ ) {
         INFO = -13
      } else if ( LDB.LT.MAX( 1, M ) ) {
         INFO = -15
      }

      if ( INFO.NE.0 ) {
         xerbla('DTPMQRT', -INFO );
         RETURN
      }

      // .. Quick return if possible ..

      if (M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0) RETURN;

      if ( LEFT .AND. TRAN ) {

         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( M-L+I+IB-1, M )
            if ( I.GE.L ) {
               LB = 0
            } else {
               LB = MB-M+L-I+1
            }
            dtprfb('L', 'T', 'F', 'C', MB, N, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB );
         }

      } else if ( RIGHT .AND. NOTRAN ) {

         DO I = 1, K, NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( N-L+I+IB-1, N )
            if ( I.GE.L ) {
               LB = 0
            } else {
               LB = MB-N+L-I+1
            }
            dtprfb('R', 'N', 'F', 'C', M, MB, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M );
         }

      } else if ( LEFT .AND. NOTRAN ) {

         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( M-L+I+IB-1, M )
            if ( I.GE.L ) {
               LB = 0
            } else {
               LB = MB-M+L-I+1
            }
            dtprfb('L', 'N', 'F', 'C', MB, N, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB );
         }

      } else if ( RIGHT .AND. TRAN ) {

         KF = ((K-1)/NB)*NB+1
         DO I = KF, 1, -NB
            IB = MIN( NB, K-I+1 )
            MB = MIN( N-L+I+IB-1, N )
            if ( I.GE.L ) {
               LB = 0
            } else {
               LB = MB-N+L-I+1
            }
            dtprfb('R', 'T', 'F', 'C', M, MB, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M );
         }

      }

      RETURN

      // End of DTPMQRT

      }
