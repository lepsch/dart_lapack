      void dormr2(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                INFO, K, LDA, LDC, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LEFT, NOTRAN;
      int                I, I1, I2, I3, MI, NI, NQ;
      double             AII;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

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
      } else if ( LDA < max( 1, K ) ) {
         INFO = -7;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('DORMR2', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0 || K == 0) return;

      if ( ( LEFT && !NOTRAN ) || ( !LEFT && NOTRAN ) ) {
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

         AII = A( I, NQ-K+I );
         A[I, NQ-K+I] = ONE;
         dlarf(SIDE, MI, NI, A( I, 1 ), LDA, TAU( I ), C, LDC, WORK );
         A[I, NQ-K+I] = AII;
      } // 10
      return;
      }