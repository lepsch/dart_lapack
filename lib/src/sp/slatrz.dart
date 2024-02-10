      void slatrz(M, N, L, final Matrix<double> A, final int LDA, TAU, WORK ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                L, LDA, M, N;
      double               A( LDA, * ), TAU( * ), WORK( * );
      // ..

      double               ZERO;
      const              ZERO = 0.0 ;
      int                I;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARFG, SLARZ

      // Test the input arguments

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

         slarfg(L+1, A( I, I ), A( I, N-L+1 ), LDA, TAU( I ) );

         // Apply H(i) to A(1:i-1,i:n) from the right

         slarz('Right', I-1, N-I+1, L, A( I, N-L+1 ), LDA, TAU( I ), A( 1, I ), LDA, WORK );

      } // 20

      }
