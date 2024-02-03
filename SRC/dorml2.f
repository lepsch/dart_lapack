      void dorml2(SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO ) {

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
      int                I, I1, I2, I3, IC, JC, MI, NI, NQ;
      double             AII;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      // EXTERNAL LSAME
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
      LEFT = LSAME( SIDE, 'L' );
      NOTRAN = LSAME( TRANS, 'N' );

      // NQ is the order of Q

      if ( LEFT ) {
         NQ = M;
      } else {
         NQ = N;
      }
      if ( !LEFT && !LSAME( SIDE, 'R' ) ) {
         INFO = -1;
      } else if ( !NOTRAN && !LSAME( TRANS, 'T' ) ) {
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
         xerbla('DORML2', -INFO );
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
         JC = 1;
      } else {
         MI = M;
         IC = 1;
      }

      for (I = I1; I3 < 0 ? I >= I2 : I <= I2; I += I3) { // 10
         if ( LEFT ) {

            // H(i) is applied to C(i:m,1:n)

            MI = M - I + 1;
            IC = I;
         } else {

            // H(i) is applied to C(1:m,i:n)

            NI = N - I + 1;
            JC = I;
         }

         // Apply H(i)

         AII = A( I, I );
         A( I, I ) = ONE;
         dlarf(SIDE, MI, NI, A( I, I ), LDA, TAU( I ), C( IC, JC ), LDC, WORK );
         A( I, I ) = AII;
      } // 10
      return;
      }
