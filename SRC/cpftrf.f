      SUBROUTINE CPFTRF( TRANSR, UPLO, N, A, INFO );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                N, INFO;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( 0: * );

*  =====================================================================

      // .. Parameters ..
      REAL               ONE;
      COMPLEX            CONE;
      const              ONE = 1.0, CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, NISODD, NORMALTRANSR;
      int                N1, N2, K;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CHERK, CPOTRF, CTRSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NORMALTRANSR = LSAME( TRANSR, 'N' );
      LOWER = LSAME( UPLO, 'L' );
      if ( !NORMALTRANSR && !LSAME( TRANSR, 'C' ) ) {
         INFO = -1;
      } else if ( !LOWER && !LSAME( UPLO, 'U' ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      }
      if ( INFO != 0 ) {
         xerbla('CPFTRF', -INFO );
         RETURN;
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // If N is odd, set NISODD = true;
      // If N is even, set K = N/2 and NISODD = false;

      if ( MOD( N, 2 ) == 0 ) {
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

      // start execution: there are eight cases

      if ( NISODD ) {

         // N is odd

         if ( NORMALTRANSR ) {

            // N is odd and TRANSR = 'N'

            if ( LOWER ) {

              // SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
              // T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
              // T1 -> a(0), T2 -> a(n), S -> a(n1)

               cpotrf('L', N1, A( 0 ), N, INFO );
               if (INFO > 0) RETURN;
               ctrsm('R', 'L', 'C', 'N', N2, N1, CONE, A( 0 ), N, A( N1 ), N );
               cherk('U', 'N', N2, N1, -ONE, A( N1 ), N, ONE, A( N ), N );
               cpotrf('U', N2, A( N ), N, INFO );
               if (INFO > 0) INFO = INFO + N1;

            } else {

              // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
              // T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
              // T1 -> a(n2), T2 -> a(n1), S -> a(0)

               cpotrf('L', N1, A( N2 ), N, INFO );
               if (INFO > 0) RETURN;
               ctrsm('L', 'L', 'N', 'N', N1, N2, CONE, A( N2 ), N, A( 0 ), N );
               cherk('U', 'C', N2, N1, -ONE, A( 0 ), N, ONE, A( N1 ), N );
               cpotrf('U', N2, A( N1 ), N, INFO );
               if (INFO > 0) INFO = INFO + N1;

            }

         } else {

            // N is odd and TRANSR = 'C'

            if ( LOWER ) {

               // SRPA for LOWER, TRANSPOSE and N is odd
               // T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
               // T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1

               cpotrf('U', N1, A( 0 ), N1, INFO );
               if (INFO > 0) RETURN;
               ctrsm('L', 'U', 'C', 'N', N1, N2, CONE, A( 0 ), N1, A( N1*N1 ), N1 );
               cherk('L', 'C', N2, N1, -ONE, A( N1*N1 ), N1, ONE, A( 1 ), N1 );
               cpotrf('L', N2, A( 1 ), N1, INFO );
               if (INFO > 0) INFO = INFO + N1;

            } else {

               // SRPA for UPPER, TRANSPOSE and N is odd
               // T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
               // T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2

               cpotrf('U', N1, A( N2*N2 ), N2, INFO );
               if (INFO > 0) RETURN;
               ctrsm('R', 'U', 'N', 'N', N2, N1, CONE, A( N2*N2 ), N2, A( 0 ), N2 );
               cherk('L', 'N', N2, N1, -ONE, A( 0 ), N2, ONE, A( N1*N2 ), N2 );
               cpotrf('L', N2, A( N1*N2 ), N2, INFO );
               if (INFO > 0) INFO = INFO + N1;

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

               cpotrf('L', K, A( 1 ), N+1, INFO );
               if (INFO > 0) RETURN;
               ctrsm('R', 'L', 'C', 'N', K, K, CONE, A( 1 ), N+1, A( K+1 ), N+1 );
               cherk('U', 'N', K, K, -ONE, A( K+1 ), N+1, ONE, A( 0 ), N+1 );
               cpotrf('U', K, A( 0 ), N+1, INFO );
               if (INFO > 0) INFO = INFO + K;

            } else {

               // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
               // T1 -> a(k+1), T2 -> a(k), S -> a(0)

               cpotrf('L', K, A( K+1 ), N+1, INFO );
               if (INFO > 0) RETURN;
               ctrsm('L', 'L', 'N', 'N', K, K, CONE, A( K+1 ), N+1, A( 0 ), N+1 );
               cherk('U', 'C', K, K, -ONE, A( 0 ), N+1, ONE, A( K ), N+1 );
               cpotrf('U', K, A( K ), N+1, INFO );
               if (INFO > 0) INFO = INFO + K;

            }

         } else {

            // N is even and TRANSR = 'C'

            if ( LOWER ) {

               // SRPA for LOWER, TRANSPOSE and N is even (see paper)
               // T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
               // T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k

               cpotrf('U', K, A( 0+K ), K, INFO );
               if (INFO > 0) RETURN;
               ctrsm('L', 'U', 'C', 'N', K, K, CONE, A( K ), N1, A( K*( K+1 ) ), K );
               cherk('L', 'C', K, K, -ONE, A( K*( K+1 ) ), K, ONE, A( 0 ), K );
               cpotrf('L', K, A( 0 ), K, INFO );
               if (INFO > 0) INFO = INFO + K;

            } else {

               // SRPA for UPPER, TRANSPOSE and N is even (see paper)
               // T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
               // T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k

               cpotrf('U', K, A( K*( K+1 ) ), K, INFO );
               if (INFO > 0) RETURN;
               ctrsm('R', 'U', 'N', 'N', K, K, CONE, A( K*( K+1 ) ), K, A( 0 ), K );
               cherk('L', 'N', K, K, -ONE, A( 0 ), K, ONE, A( K*K ), K );
               cpotrf('L', K, A( K*K ), K, INFO );
               if (INFO > 0) INFO = INFO + K;

            }

         }

      }

      RETURN;

      // End of CPFTRF

      }
