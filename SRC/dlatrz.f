      SUBROUTINE DLATRZ( M, N, L, A, LDA, TAU, WORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                L, LDA, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFG, DLARZ
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      // Quick return if possible

      if ( M == 0 ) {
         RETURN
      } else if ( M == N ) {
         for (I = 1; I <= N; I++) { // 10
            TAU( I ) = ZERO
         } // 10
         RETURN
      }

      DO 20 I = M, 1, -1

         // Generate elementary reflector H(i) to annihilate
         // [ A(i,i) A(i,n-l+1:n) ]

         dlarfg(L+1, A( I, I ), A( I, N-L+1 ), LDA, TAU( I ) );

         // Apply H(i) to A(1:i-1,i:n) from the right

         dlarz('Right', I-1, N-I+1, L, A( I, N-L+1 ), LDA, TAU( I ), A( 1, I ), LDA, WORK );

      } // 20

      RETURN

      // End of DLATRZ

      }
