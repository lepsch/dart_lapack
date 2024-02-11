      void cunm2r(final int SIDE, final int TRANS, final int M, final int N, final int K, final Matrix<double> A, final int LDA, final int TAU, final Matrix<double> C, final int LDC, final Array<double> _WORK, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE, TRANS;
      int                INFO, K, LDA, LDC, M, N;
      Complex            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

      Complex            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      bool               LEFT, NOTRAN;
      int                I, I1, I2, I3, IC, JC, MI, NI, NQ;
      Complex            AII, TAUI;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX

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
      } else if ( LDA < max( 1, NQ ) ) {
         INFO = -7;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('CUNM2R', -INFO );
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
         JC = 1;
      } else {
         MI = M;
         IC = 1;
      }

      for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) { // 10
         if ( LEFT ) {

            // H(i) or H(i)**H is applied to C(i:m,1:n)

            MI = M - I + 1;
            IC = I;
         } else {

            // H(i) or H(i)**H is applied to C(1:m,i:n)

            NI = N - I + 1;
            JC = I;
         }

         // Apply H(i) or H(i)**H

         if ( NOTRAN ) {
            TAUI = TAU( I );
         } else {
            TAUI = CONJG( TAU( I ) );
         }
         AII = A( I, I );
         A[I][I] = ONE;
         clarf(SIDE, MI, NI, A( I, I ), 1, TAUI, C( IC, JC ), LDC, WORK );
         A[I][I] = AII;
      } // 10
      }
