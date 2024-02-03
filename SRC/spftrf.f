      SUBROUTINE SPFTRF( TRANSR, UPLO, N, A, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                N, INFO;
      // ..
      // .. Array Arguments ..
      REAL               A( 0: * )

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
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
      // EXTERNAL XERBLA, SSYRK, SPOTRF, STRSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      if ( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'T' ) ) {
         INFO = -1
      } else if ( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      }
      if ( INFO.NE.0 ) {
         xerbla('SPFTRF', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // If N is odd, set NISODD = true;
      // If N is even, set K = N/2 and NISODD = false;

      if ( MOD( N, 2 ) == 0 ) {
         K = N / 2
         NISODD = false;
      } else {
         NISODD = true;
      }

      // Set N1 and N2 depending on LOWER

      if ( LOWER ) {
         N2 = N / 2
         N1 = N - N2
      } else {
         N1 = N / 2
         N2 = N - N1
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

               spotrf('L', N1, A( 0 ), N, INFO );
               if (INFO.GT.0) RETURN;
               strsm('R', 'L', 'T', 'N', N2, N1, ONE, A( 0 ), N, A( N1 ), N );
               ssyrk('U', 'N', N2, N1, -ONE, A( N1 ), N, ONE, A( N ), N );
               spotrf('U', N2, A( N ), N, INFO );
               if (INFO.GT.0) INFO = INFO + N1;

            } else {

              // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
              // T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
              // T1 -> a(n2), T2 -> a(n1), S -> a(0)

               spotrf('L', N1, A( N2 ), N, INFO );
               if (INFO.GT.0) RETURN;
               strsm('L', 'L', 'N', 'N', N1, N2, ONE, A( N2 ), N, A( 0 ), N );
               ssyrk('U', 'T', N2, N1, -ONE, A( 0 ), N, ONE, A( N1 ), N );
               spotrf('U', N2, A( N1 ), N, INFO );
               if (INFO.GT.0) INFO = INFO + N1;

            }

         } else {

            // N is odd and TRANSR = 'T'

            if ( LOWER ) {

               // SRPA for LOWER, TRANSPOSE and N is odd
               // T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
               // T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1

               spotrf('U', N1, A( 0 ), N1, INFO );
               if (INFO.GT.0) RETURN;
               strsm('L', 'U', 'T', 'N', N1, N2, ONE, A( 0 ), N1, A( N1*N1 ), N1 );
               ssyrk('L', 'T', N2, N1, -ONE, A( N1*N1 ), N1, ONE, A( 1 ), N1 );
               spotrf('L', N2, A( 1 ), N1, INFO );
               if (INFO.GT.0) INFO = INFO + N1;

            } else {

               // SRPA for UPPER, TRANSPOSE and N is odd
               // T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
               // T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2

               spotrf('U', N1, A( N2*N2 ), N2, INFO );
               if (INFO.GT.0) RETURN;
               strsm('R', 'U', 'N', 'N', N2, N1, ONE, A( N2*N2 ), N2, A( 0 ), N2 );
               ssyrk('L', 'N', N2, N1, -ONE, A( 0 ), N2, ONE, A( N1*N2 ), N2 );
               spotrf('L', N2, A( N1*N2 ), N2, INFO );
               if (INFO.GT.0) INFO = INFO + N1;

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

               spotrf('L', K, A( 1 ), N+1, INFO );
               if (INFO.GT.0) RETURN;
               strsm('R', 'L', 'T', 'N', K, K, ONE, A( 1 ), N+1, A( K+1 ), N+1 );
               ssyrk('U', 'N', K, K, -ONE, A( K+1 ), N+1, ONE, A( 0 ), N+1 );
               spotrf('U', K, A( 0 ), N+1, INFO );
               if (INFO.GT.0) INFO = INFO + K;

            } else {

               // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
               // T1 -> a(k+1), T2 -> a(k), S -> a(0)

               spotrf('L', K, A( K+1 ), N+1, INFO );
               if (INFO.GT.0) RETURN;
               strsm('L', 'L', 'N', 'N', K, K, ONE, A( K+1 ), N+1, A( 0 ), N+1 );
               ssyrk('U', 'T', K, K, -ONE, A( 0 ), N+1, ONE, A( K ), N+1 );
               spotrf('U', K, A( K ), N+1, INFO );
               if (INFO.GT.0) INFO = INFO + K;

            }

         } else {

            // N is even and TRANSR = 'T'

            if ( LOWER ) {

               // SRPA for LOWER, TRANSPOSE and N is even (see paper)
               // T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
               // T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k

               spotrf('U', K, A( 0+K ), K, INFO );
               if (INFO.GT.0) RETURN;
               strsm('L', 'U', 'T', 'N', K, K, ONE, A( K ), N1, A( K*( K+1 ) ), K );
               ssyrk('L', 'T', K, K, -ONE, A( K*( K+1 ) ), K, ONE, A( 0 ), K );
               spotrf('L', K, A( 0 ), K, INFO );
               if (INFO.GT.0) INFO = INFO + K;

            } else {

               // SRPA for UPPER, TRANSPOSE and N is even (see paper)
               // T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
               // T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k

               spotrf('U', K, A( K*( K+1 ) ), K, INFO );
               if (INFO.GT.0) RETURN;
               strsm('R', 'U', 'N', 'N', K, K, ONE, A( K*( K+1 ) ), K, A( 0 ), K );
               ssyrk('L', 'N', K, K, -ONE, A( 0 ), K, ONE, A( K*K ), K );
               spotrf('L', K, A( K*K ), K, INFO );
               if (INFO.GT.0) INFO = INFO + K;

            }

         }

      }

      RETURN

      // End of SPFTRF

      }
