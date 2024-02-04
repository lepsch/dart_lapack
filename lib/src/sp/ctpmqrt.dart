      void ctpmqrt(SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT, A, LDA, B, LDB, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String    SIDE, TRANS;
      int       INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT;
      // ..
      // .. Array Arguments ..
      Complex   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * );
      // ..

// =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LEFT, RIGHT, TRAN, NOTRAN;
      int                I, IB, MB, LB, KF, LDAQ, LDVQ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTPRFB, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // .. Test the input arguments ..

      INFO   = 0;
      LEFT   = lsame( SIDE,  'L' );
      RIGHT  = lsame( SIDE,  'R' );
      TRAN   = lsame( TRANS, 'C' );
      NOTRAN = lsame( TRANS, 'N' );

      if ( LEFT ) {
         LDVQ = max( 1, M );
         LDAQ = max( 1, K );
      } else if ( RIGHT ) {
         LDVQ = max( 1, N );
         LDAQ = max( 1, M );
      }
      if ( !LEFT && !RIGHT ) {
         INFO = -1;
      } else if ( !TRAN && !NOTRAN ) {
         INFO = -2;
      } else if ( M < 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( K < 0 ) {
         INFO = -5;
      } else if ( L < 0 || L > K ) {
         INFO = -6;
      } else if ( NB < 1 || (NB > K && K > 0) ) {
         INFO = -7;
      } else if ( LDV < LDVQ ) {
         INFO = -9;
      } else if ( LDT < NB ) {
         INFO = -11;
      } else if ( LDA < LDAQ ) {
         INFO = -13;
      } else if ( LDB < max( 1, M ) ) {
         INFO = -15;
      }

      if ( INFO != 0 ) {
         xerbla('CTPMQRT', -INFO );
         return;
      }

      // .. Quick return if possible ..

      if (M == 0 || N == 0 || K == 0) return;

      if ( LEFT && TRAN ) {

         for (I = 1; NB < 0 ? I >= K : I <= K; I += NB) {
            IB = min( NB, K-I+1 );
            MB = min( M-L+I+IB-1, M );
            if ( I >= L ) {
               LB = 0;
            } else {
               LB = MB-M+L-I+1;
            }
            ctprfb('L', 'C', 'F', 'C', MB, N, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB );
         }

      } else if ( RIGHT && NOTRAN ) {

         for (I = 1; NB < 0 ? I >= K : I <= K; I += NB) {
            IB = min( NB, K-I+1 );
            MB = min( N-L+I+IB-1, N );
            if ( I >= L ) {
               LB = 0;
            } else {
               LB = MB-N+L-I+1;
            }
            ctprfb('R', 'N', 'F', 'C', M, MB, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M );
         }

      } else if ( LEFT && NOTRAN ) {

         KF = ((K-1)/NB)*NB+1;
         for (I = KF; -NB < 0 ? I >= 1 : I <= 1; I += -NB) {
            IB = min( NB, K-I+1 );
            MB = min( M-L+I+IB-1, M );
            if ( I >= L ) {
               LB = 0;
            } else {
               LB = MB-M+L-I+1;
            }
            ctprfb('L', 'N', 'F', 'C', MB, N, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB );
         }

      } else if ( RIGHT && TRAN ) {

         KF = ((K-1)/NB)*NB+1;
         for (I = KF; -NB < 0 ? I >= 1 : I <= 1; I += -NB) {
            IB = min( NB, K-I+1 );
            MB = min( N-L+I+IB-1, N );
            if ( I >= L ) {
               LB = 0;
            } else {
               LB = MB-N+L-I+1;
            }
            ctprfb('R', 'C', 'F', 'C', M, MB, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M );
         }

      }

      return;
      }
