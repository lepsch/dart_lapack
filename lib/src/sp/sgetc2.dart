      void sgetc2(N, A, LDA, IPIV, JPIV, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), JPIV( * );
      double               A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, IP, IPV, J, JP, JPV;
      double               BIGNUM, EPS, SMIN, SMLNUM, XMAX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGER, SSWAP
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

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
            A[1][1] = SMLNUM;
         }
         return;
      }

      // Factorize A using complete pivoting.
      // Set pivots less than SMIN to SMIN.

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

         if (IPV != I) sswap( N, A( IPV, 1 ), LDA, A( I, 1 ), LDA );
         IPIV[I] = IPV;

         // Swap columns

         if (JPV != I) sswap( N, A( 1, JPV ), 1, A( 1, I ), 1 );
         JPIV[I] = JPV;

         // Check for singularity

         if ( ( A( I, I ) ).abs() < SMIN ) {
            INFO = I;
            A[I][I] = SMIN;
         }
         for (J = I + 1; J <= N; J++) { // 30
            A[J][I] = A( J, I ) / A( I, I );
         } // 30
         sger(N-I, N-I, -ONE, A( I+1, I ), 1, A( I, I+1 ), LDA, A( I+1, I+1 ), LDA );
      } // 40

      if ( ( A( N, N ) ).abs() < SMIN ) {
         INFO = N;
         A[N][N] = SMIN;
      }

      // Set last pivots to N

      IPIV[N] = N;
      JPIV[N] = N;

      return;
      }
