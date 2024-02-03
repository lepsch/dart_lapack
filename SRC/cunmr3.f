      void cunmr3(SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                INFO, K, L, LDA, LDC, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LEFT, NOTRAN;
      int                I, I1, I2, I3, IC, JA, JC, MI, NI, NQ;
      COMPLEX            TAUI;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARZ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX
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
      } else if ( !NOTRAN && !LSAME( TRANS, 'C' ) ) {
         INFO = -2;
      } else if ( M < 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( K < 0 || K > NQ ) {
         INFO = -5;
      } else if ( L < 0 || ( LEFT && ( L > M ) ) || ( !LEFT && ( L > N ) ) ) {
         INFO = -6;
      } else if ( LDA < max( 1, K ) ) {
         INFO = -8;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('CUNMR3', -INFO );
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
         JA = M - L + 1;
         JC = 1;
      } else {
         MI = M;
         JA = N - L + 1;
         IC = 1;
      }

      DO 10 I = I1, I2, I3;
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
         clarz(SIDE, MI, NI, L, A( I, JA ), LDA, TAUI, C( IC, JC ), LDC, WORK );

      } // 10

      return;
      }
