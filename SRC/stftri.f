      SUBROUTINE STFTRI( TRANSR, UPLO, DIAG, N, A, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             TRANSR, UPLO, DIAG;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      REAL               A( 0: * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
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
      // EXTERNAL XERBLA, STRMM, STRTRI
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'T' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LSAME( DIAG, 'N' ) .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STFTRI', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      // If N is odd, set NISODD = .TRUE.
      // If N is even, set K = N/2 and NISODD = .FALSE.
*
      IF( MOD( N, 2 ).EQ.0 ) THEN
         K = N / 2
         NISODD = .FALSE.
      ELSE
         NISODD = .TRUE.
      END IF
*
      // Set N1 and N2 depending on LOWER
*
      IF( LOWER ) THEN
         N2 = N / 2
         N1 = N - N2
      ELSE
         N1 = N / 2
         N2 = N - N1
      END IF
*
*
      // start execution: there are eight cases
*
      IF( NISODD ) THEN
*
         // N is odd
*
         IF( NORMALTRANSR ) THEN
*
            // N is odd and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
              // SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
              // T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
              // T1 -> a(0), T2 -> a(n), S -> a(n1)
*
               CALL STRTRI( 'L', DIAG, N1, A( 0 ), N, INFO )
               IF( INFO.GT.0 ) RETURN                CALL STRMM( 'R', 'L', 'N', DIAG, N2, N1, -ONE, A( 0 ), N, A( N1 ), N )
               CALL STRTRI( 'U', DIAG, N2, A( N ), N, INFO )
               IF( INFO.GT.0 ) INFO = INFO + N1                IF( INFO.GT.0 ) RETURN                CALL STRMM( 'L', 'U', 'T', DIAG, N2, N1, ONE, A( N ), N, A( N1 ), N )
*
            ELSE
*
              // SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
              // T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
              // T1 -> a(n2), T2 -> a(n1), S -> a(0)
*
               CALL STRTRI( 'L', DIAG, N1, A( N2 ), N, INFO )
               IF( INFO.GT.0 ) RETURN                CALL STRMM( 'L', 'L', 'T', DIAG, N1, N2, -ONE, A( N2 ), N, A( 0 ), N )
               CALL STRTRI( 'U', DIAG, N2, A( N1 ), N, INFO )
               IF( INFO.GT.0 ) INFO = INFO + N1                IF( INFO.GT.0 ) RETURN                CALL STRMM( 'R', 'U', 'N', DIAG, N1, N2, ONE, A( N1 ), N, A( 0 ), N )
*
            END IF
*
         ELSE
*
            // N is odd and TRANSR = 'T'
*
            IF( LOWER ) THEN
*
               // SRPA for LOWER, TRANSPOSE and N is odd
               // T1 -> a(0), T2 -> a(1), S -> a(0+n1*n1)
*
               CALL STRTRI( 'U', DIAG, N1, A( 0 ), N1, INFO )
               IF( INFO.GT.0 ) RETURN                CALL STRMM( 'L', 'U', 'N', DIAG, N1, N2, -ONE, A( 0 ), N1, A( N1*N1 ), N1 )
               CALL STRTRI( 'L', DIAG, N2, A( 1 ), N1, INFO )
               IF( INFO.GT.0 ) INFO = INFO + N1                IF( INFO.GT.0 ) RETURN                CALL STRMM( 'R', 'L', 'T', DIAG, N1, N2, ONE, A( 1 ), N1, A( N1*N1 ), N1 )
*
            ELSE
*
               // SRPA for UPPER, TRANSPOSE and N is odd
               // T1 -> a(0+n2*n2), T2 -> a(0+n1*n2), S -> a(0)
*
               CALL STRTRI( 'U', DIAG, N1, A( N2*N2 ), N2, INFO )
               IF( INFO.GT.0 ) RETURN                CALL STRMM( 'R', 'U', 'T', DIAG, N2, N1, -ONE, A( N2*N2 ), N2, A( 0 ), N2 )
               CALL STRTRI( 'L', DIAG, N2, A( N1*N2 ), N2, INFO )
               IF( INFO.GT.0 ) INFO = INFO + N1                IF( INFO.GT.0 ) RETURN                CALL STRMM( 'L', 'L', 'N', DIAG, N2, N1, ONE, A( N1*N2 ), N2, A( 0 ), N2 )
            END IF
*
         END IF
*
      ELSE
*
         // N is even
*
         IF( NORMALTRANSR ) THEN
*
            // N is even and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
               // SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
               // T1 -> a(1), T2 -> a(0), S -> a(k+1)
*
               CALL STRTRI( 'L', DIAG, K, A( 1 ), N+1, INFO )
               IF( INFO.GT.0 ) RETURN                CALL STRMM( 'R', 'L', 'N', DIAG, K, K, -ONE, A( 1 ), N+1, A( K+1 ), N+1 )
               CALL STRTRI( 'U', DIAG, K, A( 0 ), N+1, INFO )
               IF( INFO.GT.0 ) INFO = INFO + K                IF( INFO.GT.0 ) RETURN                CALL STRMM( 'L', 'U', 'T', DIAG, K, K, ONE, A( 0 ), N+1, A( K+1 ), N+1 )
*
            ELSE
*
               // SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
               // T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
               // T1 -> a(k+1), T2 -> a(k), S -> a(0)
*
               CALL STRTRI( 'L', DIAG, K, A( K+1 ), N+1, INFO )
               IF( INFO.GT.0 ) RETURN                CALL STRMM( 'L', 'L', 'T', DIAG, K, K, -ONE, A( K+1 ), N+1, A( 0 ), N+1 )
               CALL STRTRI( 'U', DIAG, K, A( K ), N+1, INFO )
               IF( INFO.GT.0 ) INFO = INFO + K                IF( INFO.GT.0 ) RETURN                CALL STRMM( 'R', 'U', 'N', DIAG, K, K, ONE, A( K ), N+1, A( 0 ), N+1 )
            END IF
         ELSE
*
            // N is even and TRANSR = 'T'
*
            IF( LOWER ) THEN
*
               // SRPA for LOWER, TRANSPOSE and N is even (see paper)
               // T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
               // T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k
*
               CALL STRTRI( 'U', DIAG, K, A( K ), K, INFO )
               IF( INFO.GT.0 ) RETURN                CALL STRMM( 'L', 'U', 'N', DIAG, K, K, -ONE, A( K ), K, A( K*( K+1 ) ), K )
               CALL STRTRI( 'L', DIAG, K, A( 0 ), K, INFO )
               IF( INFO.GT.0 ) INFO = INFO + K                IF( INFO.GT.0 ) RETURN                CALL STRMM( 'R', 'L', 'T', DIAG, K, K, ONE, A( 0 ), K, A( K*( K+1 ) ), K )
            ELSE
*
               // SRPA for UPPER, TRANSPOSE and N is even (see paper)
               // T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
               // T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k
*
               CALL STRTRI( 'U', DIAG, K, A( K*( K+1 ) ), K, INFO )
               IF( INFO.GT.0 ) RETURN                CALL STRMM( 'R', 'U', 'T', DIAG, K, K, -ONE, A( K*( K+1 ) ), K, A( 0 ), K )
               CALL STRTRI( 'L', DIAG, K, A( K*K ), K, INFO )
               IF( INFO.GT.0 ) INFO = INFO + K                IF( INFO.GT.0 ) RETURN                CALL STRMM( 'L', 'L', 'N', DIAG, K, K, ONE, A( K*K ), K, A( 0 ), K )
            END IF
         END IF
      END IF
*
      RETURN
*
      // End of STFTRI
*
      END
