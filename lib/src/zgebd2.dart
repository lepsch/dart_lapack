      void zgebd2(final int M, final int N, final Matrix<double> A, final int LDA, final int D, final int E, final int TAUQ, final int TAUP, final Array<double> _WORK, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, M, N;
      double             D( * ), E( * );
      Complex         A( LDA, * ), TAUP( * ), TAUQ( * ), WORK( * );
      // ..

      Complex         ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      int                I;
      Complex         ALPHA;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLACGV, ZLARF, ZLARFG
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX, MIN

      // Test the input parameters

      INFO = 0;
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -4;
      }
      if ( INFO < 0 ) {
         xerbla('ZGEBD2', -INFO );
         return;
      }

      if ( M >= N ) {

         // Reduce to upper bidiagonal form

         for (I = 1; I <= N; I++) { // 10

            // Generate elementary reflector H(i) to annihilate A(i+1:m,i)

            ALPHA = A( I, I );
            zlarfg(M-I+1, ALPHA, A( min( I+1, M ), I ), 1, TAUQ( I ) );
            D[I] = ALPHA.toDouble();
            A[I][I] = ONE;

            // Apply H(i)**H to A(i:m,i+1:n) from the left

            if (I < N) zlarf( 'Left', M-I+1, N-I, A( I, I ), 1, DCONJG( TAUQ( I ) ), A( I, I+1 ), LDA, WORK );
            A[I][I] = D( I );

            if ( I < N ) {

               // Generate elementary reflector G(i) to annihilate
               // A(i,i+2:n)

               zlacgv(N-I, A( I, I+1 ), LDA );
               ALPHA = A( I, I+1 );
               zlarfg(N-I, ALPHA, A( I, min( I+2, N ) ), LDA, TAUP( I ) );
               E[I] = ALPHA.toDouble();
               A[I][I+1] = ONE;

               // Apply G(i) to A(i+1:m,i+1:n) from the right

               zlarf('Right', M-I, N-I, A( I, I+1 ), LDA, TAUP( I ), A( I+1, I+1 ), LDA, WORK );
               zlacgv(N-I, A( I, I+1 ), LDA );
               A[I][I+1] = E( I );
            } else {
               TAUP[I] = ZERO;
            }
         } // 10
      } else {

         // Reduce to lower bidiagonal form

         for (I = 1; I <= M; I++) { // 20

            // Generate elementary reflector G(i) to annihilate A(i,i+1:n)

            zlacgv(N-I+1, A( I, I ), LDA );
            ALPHA = A( I, I );
            zlarfg(N-I+1, ALPHA, A( I, min( I+1, N ) ), LDA, TAUP( I ) );
            D[I] = ALPHA.toDouble();
            A[I][I] = ONE;

            // Apply G(i) to A(i+1:m,i:n) from the right

            if (I < M) zlarf( 'Right', M-I, N-I+1, A( I, I ), LDA, TAUP( I ), A( I+1, I ), LDA, WORK );
            zlacgv(N-I+1, A( I, I ), LDA );
            A[I][I] = D( I );

            if ( I < M ) {

               // Generate elementary reflector H(i) to annihilate
               // A(i+2:m,i)

               ALPHA = A( I+1, I );
               zlarfg(M-I, ALPHA, A( min( I+2, M ), I ), 1, TAUQ( I ) );
               E[I] = ALPHA.toDouble();
               A[I+1][I] = ONE;

               // Apply H(i)**H to A(i+1:m,i+1:n) from the left

               zlarf('Left', M-I, N-I, A( I+1, I ), 1, DCONJG( TAUQ( I ) ), A( I+1, I+1 ), LDA, WORK );
               A[I+1][I] = E( I );
            } else {
               TAUQ[I] = ZERO;
            }
         } // 20
      }
      }
