      void zunmr2(final int SIDE, final int TRANS, final int M, final int N, final int K, final Matrix<double> A_, final int LDA, final int TAU, final Matrix<double> C_, final int LDC, final Array<double> _WORK_, final Box<int> INFO,) {
  final A = A_.dim();
  final C = C_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE, TRANS;
      int                INFO, K, LDA, LDC, M, N;
      Complex         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

      Complex         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      bool               LEFT, NOTRAN;
      int                I, I1, I2, I3, MI, NI, NQ;
      Complex         AII, TAUI;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACGV, ZLARF
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX

      // Test the input arguments

      INFO = 0;
      LEFT = lsame( SIDE, 'L' );
      NOTRAN = lsame( TRANS, 'N' );

      // NQ is the order of Q

      if ( LEFT ) {
         NQ = M;
      } else {
         NQ = N;
      }
      if ( !LEFT && !lsame( SIDE, 'R' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !lsame( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( M < 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( K < 0 || K > NQ ) {
         INFO = -5;
      } else if ( LDA < max( 1, K ) ) {
         INFO = -7;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('ZUNMR2', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0 || K == 0) return;

      if ( ( LEFT && !NOTRAN || !LEFT && NOTRAN ) ) {
         I1 = 1;
         I2 = K;
         I3 = 1;
      } else {
         I1 = K;
         I2 = 1;
         I3 = -1;
      }

      if ( LEFT ) {
         NI = N;
      } else {
         MI = M;
      }

      for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) { // 10
         if ( LEFT ) {

            // H(i) or H(i)**H is applied to C(1:m-k+i,1:n)

            MI = M - K + I;
         } else {

            // H(i) or H(i)**H is applied to C(1:m,1:n-k+i)

            NI = N - K + I;
         }

         // Apply H(i) or H(i)**H

         if ( NOTRAN ) {
            TAUI = DCONJG( TAU( I ) );
         } else {
            TAUI = TAU( I );
         }
         zlacgv(NQ-K+I-1, A( I, 1 ), LDA );
         AII = A( I, NQ-K+I );
         A[I][NQ-K+I] = ONE;
         zlarf(SIDE, MI, NI, A( I, 1 ), LDA, TAUI, C, LDC, WORK );
         A[I][NQ-K+I] = AII;
         zlacgv(NQ-K+I-1, A( I, 1 ), LDA );
      } // 10
      }
