      SUBROUTINE ZUNM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                INFO, K, LDA, LDC, M, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LEFT, NOTRAN;
      int                I, I1, I2, I3, MI, NI, NQ;
      COMPLEX*16         AII, TAUI;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLARF
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX
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
      } else if ( LDA < MAX( 1, NQ ) ) {
         INFO = -7;
      } else if ( LDC < MAX( 1, M ) ) {
         INFO = -10;
      }
      if ( INFO != 0 ) {
         xerbla('ZUNM2L', -INFO );
         return;
      }

      // Quick return if possible

      if (M == 0 || N == 0 || K == 0) return;

      if ( ( LEFT && NOTRAN || !LEFT && !NOTRAN ) ) {
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

      DO 10 I = I1, I2, I3;
         if ( LEFT ) {

            // H(i) or H(i)**H is applied to C(1:m-k+i,1:n)

            MI = M - K + I;
         } else {

            // H(i) or H(i)**H is applied to C(1:m,1:n-k+i)

            NI = N - K + I;
         }

         // Apply H(i) or H(i)**H

         if ( NOTRAN ) {
            TAUI = TAU( I );
         } else {
            TAUI = DCONJG( TAU( I ) );
         }
         AII = A( NQ-K+I, I );
         A( NQ-K+I, I ) = ONE;
         zlarf(SIDE, MI, NI, A( 1, I ), 1, TAUI, C, LDC, WORK );
         A( NQ-K+I, I ) = AII;
      } // 10
      return;

      // End of ZUNM2L

      }
