      void sorm2l(SIDE, TRANS, M, N, K, final Matrix<double> A, final int LDA, TAU, final Matrix<double> C, final int LDC, WORK, Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE, TRANS;
      int                INFO, K, LDA, LDC, M, N;
      double               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      bool               LEFT, NOTRAN;
      int                I, I1, I2, I3, MI, NI, NQ;
      double               AII;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

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
      } else if ( !NOTRAN && !lsame( TRANS, 'T' ) ) {
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
         xerbla('SORM2L', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0 || K == 0) return;

      if ( ( LEFT && NOTRAN ) || ( !LEFT && !NOTRAN ) ) {
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

            // H(i) is applied to C(1:m-k+i,1:n)

            MI = M - K + I;
         } else {

            // H(i) is applied to C(1:m,1:n-k+i)

            NI = N - K + I;
         }

         // Apply H(i)

         AII = A( NQ-K+I, I );
         A[NQ-K+I][I] = ONE;
         slarf(SIDE, MI, NI, A( 1, I ), 1, TAU( I ), C, LDC, WORK );
         A[NQ-K+I][I] = AII;
      } // 10
      }
