      void stpmqrt(final int SIDE, final int TRANS, final int M, final int N, final int K, final int L, final int NB, final Matrix<double> V_, final int LDV, final Matrix<double> T_, final int LDT, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final Array<double> _WORK_, final Box<int> INFO,) {
  final V = V_.dim();
  final T = T_.dim();
  final A = A_.dim();
  final B = B_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String    SIDE, TRANS;
      int       INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT;
      double   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * );
      // ..

// =====================================================================

      bool               LEFT, RIGHT, TRAN, NOTRAN;
      int                I, IB, MB, LB, KF, LDAQ, LDVQ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL STPRFB, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // .. Test the input arguments ..

      INFO   = 0;
      LEFT   = lsame( SIDE,  'L' );
      RIGHT  = lsame( SIDE,  'R' );
      TRAN   = lsame( TRANS, 'T' );
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
         xerbla('STPMQRT', -INFO );
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
            stprfb('L', 'T', 'F', 'C', MB, N, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB );
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
            stprfb('R', 'N', 'F', 'C', M, MB, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M );
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
            stprfb('L', 'N', 'F', 'C', MB, N, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( I, 1 ), LDA, B, LDB, WORK, IB );
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
            stprfb('R', 'T', 'F', 'C', M, MB, IB, LB, V( 1, I ), LDV, T( 1, I ), LDT, A( 1, I ), LDA, B, LDB, WORK, M );
         }

      }

      }
