      SUBROUTINE CGETC2( N, A, LDA, IPIV, JPIV, INFO );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), JPIV( * );
      COMPLEX            A( LDA, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, IP, IPV, J, JP, JPV;
      REAL               BIGNUM, EPS, SMIN, SMLNUM, XMAX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGERU, CSWAP
      // ..
      // .. External Functions ..
      REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, MAX
      // ..
      // .. Executable Statements ..

      INFO = 0;

      // Quick return if possible

      if (N == 0) RETURN;

      // Set constants to control overflow

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' ) / EPS;
      BIGNUM = ONE / SMLNUM;

      // Handle the case N=1 by itself

      if ( N == 1 ) {
         IPIV( 1 ) = 1;
         JPIV( 1 ) = 1;
         if ( ABS( A( 1, 1 ) ) < SMLNUM ) {
            INFO = 1;
            A( 1, 1 ) = CMPLX( SMLNUM, ZERO );
         }
         RETURN;
      }

      // Factorize A using complete pivoting.
      // Set pivots less than SMIN to SMIN

      for (I = 1; I <= N - 1; I++) { // 40

         // Find max element in matrix A

         XMAX = ZERO;
         for (IP = I; IP <= N; IP++) { // 20
            for (JP = I; JP <= N; JP++) { // 10
               if ( ABS( A( IP, JP ) ) >= XMAX ) {
                  XMAX = ABS( A( IP, JP ) );
                  IPV = IP;
                  JPV = JP;
               }
            } // 10
         } // 20
         if (I == 1) SMIN = MAX( EPS*XMAX, SMLNUM );

         // Swap rows

         if (IPV != I) CALL CSWAP( N, A( IPV, 1 ), LDA, A( I, 1 ), LDA );
         IPIV( I ) = IPV;

         // Swap columns

         if (JPV != I) CALL CSWAP( N, A( 1, JPV ), 1, A( 1, I ), 1 );
         JPIV( I ) = JPV;

         // Check for singularity

         if ( ABS( A( I, I ) ) < SMIN ) {
            INFO = I;
            A( I, I ) = CMPLX( SMIN, ZERO );
         }
         for (J = I + 1; J <= N; J++) { // 30
            A( J, I ) = A( J, I ) / A( I, I );
         } // 30
         cgeru(N-I, N-I, -CMPLX( ONE ), A( I+1, I ), 1, A( I, I+1 ), LDA, A( I+1, I+1 ), LDA );
      } // 40

      if ( ABS( A( N, N ) ) < SMIN ) {
         INFO = N;
         A( N, N ) = CMPLX( SMIN, ZERO );
      }

      // Set last pivots to N

      IPIV( N ) = N;
      JPIV( N ) = N;

      RETURN;

      // End of CGETC2

      }
