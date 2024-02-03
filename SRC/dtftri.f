      SUBROUTINE DTFTRI( TRANSR, UPLO, DIAG, N, A, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, UPLO, DIAG;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      double             A( 0: * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D+0 ;
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
      // EXTERNAL XERBLA, DTRMM, DTRTRI
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
      } else if ( .NOT.LSAME( DIAG, 'N' ) .AND. .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DTFTRI', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // If N is odd, set NISODD = .TRUE.
      // If N is even, set K = N/2 and NISODD = .FALSE.

      if ( MOD( N, 2 ).EQ.0 ) {
         K = N / 2
         NISODD = .FALSE.
      } else {
         NISODD = .TRUE.
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

               CALL DTRTRI( 'L', DIAG, N1, A( 0 ), N, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'R', 'L', 'N', DIAG, N2, N1, -ONE, A( 0 ), N, A( N1 ), N )
               CALL DTRTRI( 'U', DIAG, N2, A( N ), N, INFO )
               IF( INFO.GT.0 ) INFO = INFO + N1                IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'L', 'U', 'T', DIAG, N2, N1, ONE, A( N ), N, A( N1 ), N )

            } else {

              // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
              // T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
              // T1 -> a(n2), T2 -> a(n1), S -> a(0)

               CALL DTRTRI( 'L', DIAG, N1, A( N2 ), N, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'L', 'L', 'T', DIAG, N1, N2, -ONE, A( N2 ), N, A( 0 ), N )
               CALL DTRTRI( 'U', DIAG, N2, A( N1 ), N, INFO )
               IF( INFO.GT.0 ) INFO = INFO + N1                IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'R', 'U', 'N', DIAG, N1, N2, ONE, A( N1 ), N, A( 0 ), N )

            }

         } else {

            // N is odd and TRANSR = 'T'

            if ( LOWER ) {

               // SRPA for LOWER, TRANSPOSE and N is odd
               // T1 -> a(0), T2 -> a(1), S -> a(0+n1*n1)

               CALL DTRTRI( 'U', DIAG, N1, A( 0 ), N1, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'L', 'U', 'N', DIAG, N1, N2, -ONE, A( 0 ), N1, A( N1*N1 ), N1 )
               CALL DTRTRI( 'L', DIAG, N2, A( 1 ), N1, INFO )
               IF( INFO.GT.0 ) INFO = INFO + N1                IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'R', 'L', 'T', DIAG, N1, N2, ONE, A( 1 ), N1, A( N1*N1 ), N1 )

            } else {

               // SRPA for UPPER, TRANSPOSE and N is odd
               // T1 -> a(0+n2*n2), T2 -> a(0+n1*n2), S -> a(0)

               CALL DTRTRI( 'U', DIAG, N1, A( N2*N2 ), N2, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'R', 'U', 'T', DIAG, N2, N1, -ONE, A( N2*N2 ), N2, A( 0 ), N2 )
               CALL DTRTRI( 'L', DIAG, N2, A( N1*N2 ), N2, INFO )
               IF( INFO.GT.0 ) INFO = INFO + N1                IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'L', 'L', 'N', DIAG, N2, N1, ONE, A( N1*N2 ), N2, A( 0 ), N2 )
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

               CALL DTRTRI( 'L', DIAG, K, A( 1 ), N+1, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'R', 'L', 'N', DIAG, K, K, -ONE, A( 1 ), N+1, A( K+1 ), N+1 )
               CALL DTRTRI( 'U', DIAG, K, A( 0 ), N+1, INFO )
               IF( INFO.GT.0 ) INFO = INFO + K                IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'L', 'U', 'T', DIAG, K, K, ONE, A( 0 ), N+1, A( K+1 ), N+1 )

            } else {

               // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
               // T1 -> a(k+1), T2 -> a(k), S -> a(0)

               CALL DTRTRI( 'L', DIAG, K, A( K+1 ), N+1, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'L', 'L', 'T', DIAG, K, K, -ONE, A( K+1 ), N+1, A( 0 ), N+1 )
               CALL DTRTRI( 'U', DIAG, K, A( K ), N+1, INFO )
               IF( INFO.GT.0 ) INFO = INFO + K                IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'R', 'U', 'N', DIAG, K, K, ONE, A( K ), N+1, A( 0 ), N+1 )
            }
         } else {

            // N is even and TRANSR = 'T'

            if ( LOWER ) {

               // SRPA for LOWER, TRANSPOSE and N is even (see paper)
               // T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
               // T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k

               CALL DTRTRI( 'U', DIAG, K, A( K ), K, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'L', 'U', 'N', DIAG, K, K, -ONE, A( K ), K, A( K*( K+1 ) ), K )
               CALL DTRTRI( 'L', DIAG, K, A( 0 ), K, INFO )
               IF( INFO.GT.0 ) INFO = INFO + K                IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'R', 'L', 'T', DIAG, K, K, ONE, A( 0 ), K, A( K*( K+1 ) ), K )
            } else {

               // SRPA for UPPER, TRANSPOSE and N is even (see paper)
               // T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
               // T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k

               CALL DTRTRI( 'U', DIAG, K, A( K*( K+1 ) ), K, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'R', 'U', 'T', DIAG, K, K, -ONE, A( K*( K+1 ) ), K, A( 0 ), K )
               CALL DTRTRI( 'L', DIAG, K, A( K*K ), K, INFO )
               IF( INFO.GT.0 ) INFO = INFO + K                IF( INFO.GT.0 ) RETURN                CALL DTRMM( 'L', 'L', 'N', DIAG, K, K, ONE, A( K*K ), K, A( 0 ), K )
            }
         }
      }

      RETURN

      // End of DTFTRI

      }
