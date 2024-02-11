      void cgetc2(final int N, final Matrix<double> A, final int LDA, final Array<int> IPIV, final int JPIV, final Box<int> INFO,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, N;
      int                IPIV( * ), JPIV( * );
      Complex            A( LDA, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, IP, IPV, J, JP, JPV;
      double               BIGNUM, EPS, SMIN, SMLNUM, XMAX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGERU, CSWAP
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX

      INFO = 0;

      // Quick return if possible

      if (N == 0) return;

      // Set constants to control overflow

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Handle the case N=1 by itself

      if ( N == 1 ) {
         IPIV[1] = 1;
         JPIV[1] = 1;
         if ( ( A( 1, 1 ) ).abs() < SMLNUM ) {
            INFO = 1;
            A[1][1] = CMPLX( SMLNUM, ZERO );
         }
         return;
      }

      // Factorize A using complete pivoting.
      // Set pivots less than SMIN to SMIN

      for (I = 1; I <= N - 1; I++) { // 40

         // Find max element in matrix A

         XMAX = ZERO;
         for (IP = I; IP <= N; IP++) { // 20
            for (JP = I; JP <= N; JP++) { // 10
               if ( ( A( IP, JP ) ).abs() >= XMAX ) {
                  XMAX = ( A( IP, JP ) ).abs();
                  IPV = IP;
                  JPV = JP;
               }
            } // 10
         } // 20
         if (I == 1) SMIN = max( EPS*XMAX, SMLNUM );

         // Swap rows

         if (IPV != I) cswap( N, A( IPV, 1 ), LDA, A( I, 1 ), LDA );
         IPIV[I] = IPV;

         // Swap columns

         if (JPV != I) cswap( N, A( 1, JPV ), 1, A( 1, I ), 1 );
         JPIV[I] = JPV;

         // Check for singularity

         if ( ( A( I, I ) ).abs() < SMIN ) {
            INFO = I;
            A[I][I] = CMPLX( SMIN, ZERO );
         }
         for (J = I + 1; J <= N; J++) { // 30
            A[J][I] = A( J, I ) / A( I, I );
         } // 30
         cgeru(N-I, N-I, -CMPLX( ONE ), A( I+1, I ), 1, A( I, I+1 ), LDA, A( I+1, I+1 ), LDA );
      } // 40

      if ( ( A( N, N ) ).abs() < SMIN ) {
         INFO = N;
         A[N][N] = CMPLX( SMIN, ZERO );
      }

      // Set last pivots to N

      IPIV[N] = N;
      JPIV[N] = N;

      }
