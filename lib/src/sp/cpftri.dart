      void cpftri(TRANSR, UPLO, N, A, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                INFO, N;
      // .. Array Arguments ..
      Complex            A( 0: * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE;
      Complex            CONE;
      const              ONE = 1.0, CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, NISODD, NORMALTRANSR;
      int                N1, N2, K;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CTFTRI, CLAUUM, CTRMM, CHERK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NORMALTRANSR = lsame( TRANSR, 'N' );
      LOWER = lsame( UPLO, 'L' );
      if ( !NORMALTRANSR && !lsame( TRANSR, 'C' ) ) {
         INFO = -1;
      } else if ( !LOWER && !lsame( UPLO, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      }
      if ( INFO != 0 ) {
         xerbla('CPFTRI', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Invert the triangular Cholesky factor U or L.

      ctftri(TRANSR, UPLO, 'N', N, A, INFO );
      if (INFO > 0) return;

      // If N is odd, set NISODD = true;
      // If N is even, set K = N/2 and NISODD = false;

      if ( (N % 2) == 0 ) {
         K = N / 2;
         NISODD = false;
      } else {
         NISODD = true;
      }

      // Set N1 and N2 depending on LOWER

      if ( LOWER ) {
         N2 = N / 2;
         N1 = N - N2;
      } else {
         N1 = N / 2;
         N2 = N - N1;
      }

      // Start execution of triangular matrix multiply: inv(U)*inv(U)^C or
      // inv(L)^C*inv(L). There are eight cases.

      if ( NISODD ) {

         // N is odd

         if ( NORMALTRANSR ) {

            // N is odd and TRANSR = 'N'

            if ( LOWER ) {

               // SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:N1-1) )
               // T1 -> a(0,0), T2 -> a(0,1), S -> a(N1,0)
               // T1 -> a(0), T2 -> a(n), S -> a(N1)

               clauum('L', N1, A( 0 ), N, INFO );
               cherk('L', 'C', N1, N2, ONE, A( N1 ), N, ONE, A( 0 ), N );
               ctrmm('L', 'U', 'N', 'N', N2, N1, CONE, A( N ), N, A( N1 ), N );
               clauum('U', N2, A( N ), N, INFO );

            } else {

               // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:N2-1)
               // T1 -> a(N1+1,0), T2 -> a(N1,0), S -> a(0,0)
               // T1 -> a(N2), T2 -> a(N1), S -> a(0)

               clauum('L', N1, A( N2 ), N, INFO );
               cherk('L', 'N', N1, N2, ONE, A( 0 ), N, ONE, A( N2 ), N );
               ctrmm('R', 'U', 'C', 'N', N1, N2, CONE, A( N1 ), N, A( 0 ), N );
               clauum('U', N2, A( N1 ), N, INFO );

            }

         } else {

            // N is odd and TRANSR = 'C'

            if ( LOWER ) {

               // SRPA for LOWER, TRANSPOSE, and N is odd
               // T1 -> a(0), T2 -> a(1), S -> a(0+N1*N1)

               clauum('U', N1, A( 0 ), N1, INFO );
               cherk('U', 'N', N1, N2, ONE, A( N1*N1 ), N1, ONE, A( 0 ), N1 );
               ctrmm('R', 'L', 'N', 'N', N1, N2, CONE, A( 1 ), N1, A( N1*N1 ), N1 );
               clauum('L', N2, A( 1 ), N1, INFO );

            } else {

               // SRPA for UPPER, TRANSPOSE, and N is odd
               // T1 -> a(0+N2*N2), T2 -> a(0+N1*N2), S -> a(0)

               clauum('U', N1, A( N2*N2 ), N2, INFO );
               cherk('U', 'C', N1, N2, ONE, A( 0 ), N2, ONE, A( N2*N2 ), N2 );
               ctrmm('L', 'L', 'C', 'N', N2, N1, CONE, A( N1*N2 ), N2, A( 0 ), N2 );
               clauum('L', N2, A( N1*N2 ), N2, INFO );

            }

         }

      } else {

         // N is even

         if ( NORMALTRANSR ) {

            // N is even and TRANSR = 'N'

            if ( LOWER ) {

               // SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
               // T1 -> a(1), T2 -> a(0), S -> a(k+1)

               clauum('L', K, A( 1 ), N+1, INFO );
               cherk('L', 'C', K, K, ONE, A( K+1 ), N+1, ONE, A( 1 ), N+1 );
               ctrmm('L', 'U', 'N', 'N', K, K, CONE, A( 0 ), N+1, A( K+1 ), N+1 );
               clauum('U', K, A( 0 ), N+1, INFO );

            } else {

               // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
               // T1 -> a(k+1), T2 -> a(k), S -> a(0)

               clauum('L', K, A( K+1 ), N+1, INFO );
               cherk('L', 'N', K, K, ONE, A( 0 ), N+1, ONE, A( K+1 ), N+1 );
               ctrmm('R', 'U', 'C', 'N', K, K, CONE, A( K ), N+1, A( 0 ), N+1 );
               clauum('U', K, A( K ), N+1, INFO );

            }

         } else {

            // N is even and TRANSR = 'C'

            if ( LOWER ) {

               // SRPA for LOWER, TRANSPOSE, and N is even (see paper)
               // T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1),
               // T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k

               clauum('U', K, A( K ), K, INFO );
               cherk('U', 'N', K, K, ONE, A( K*( K+1 ) ), K, ONE, A( K ), K );
               ctrmm('R', 'L', 'N', 'N', K, K, CONE, A( 0 ), K, A( K*( K+1 ) ), K );
               clauum('L', K, A( 0 ), K, INFO );

            } else {

               // SRPA for UPPER, TRANSPOSE, and N is even (see paper)
               // T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0),
               // T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k

               clauum('U', K, A( K*( K+1 ) ), K, INFO );
               cherk('U', 'C', K, K, ONE, A( 0 ), K, ONE, A( K*( K+1 ) ), K );
               ctrmm('L', 'L', 'C', 'N', K, K, CONE, A( K*K ), K, A( 0 ), K );
               clauum('L', K, A( K*K ), K, INFO );

            }

         }

      }

      return;
      }
