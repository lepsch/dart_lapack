      void clatrz(M, N, L, A, LDA, TAU, WORK ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                L, LDA, M, N;
      // ..
      // .. Array Arguments ..
      Complex            A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      Complex            ZERO;
      const              ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I;
      Complex            ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLACGV, CLARFG, CLARZ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( M == 0 ) {
         return;
      } else if ( M == N ) {
         for (I = 1; I <= N; I++) { // 10
            TAU[I] = ZERO;
         } // 10
         return;
      }

      for (I = M; I >= 1; I--) { // 20

         // Generate elementary reflector H(i) to annihilate
         // [ A(i,i) A(i,n-l+1:n) ]

         clacgv(L, A( I, N-L+1 ), LDA );
         ALPHA = CONJG( A( I, I ) );
         clarfg(L+1, ALPHA, A( I, N-L+1 ), LDA, TAU( I ) );
         TAU[I] = CONJG( TAU( I ) );

         // Apply H(i) to A(1:i-1,i:n) from the right

         clarz('Right', I-1, N-I+1, L, A( I, N-L+1 ), LDA, CONJG( TAU( I ) ), A( 1, I ), LDA, WORK );
         A[I][I] = CONJG( ALPHA );

      } // 20

      return;
      }
