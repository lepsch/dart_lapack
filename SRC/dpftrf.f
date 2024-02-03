      SUBROUTINE DPFTRF( TRANSR, UPLO, N, A, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TRANSR, UPLO;
      int                N, INFO
*     ..
*     .. Array Arguments ..
      double             A( 0: * );
*
*  =====================================================================
*
*     .. Parameters ..
      double             ONE;
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      bool               LOWER, NISODD, NORMALTRANSR;
      int                N1, N2, K
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, DSYRK, DPOTRF, DTRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'T' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPFTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     If N is odd, set NISODD = .TRUE.
*     If N is even, set K = N/2 and NISODD = .FALSE.
*
      IF( MOD( N, 2 ).EQ.0 ) THEN
         K = N / 2
         NISODD = .FALSE.
      ELSE
         NISODD = .TRUE.
      END IF
*
*     Set N1 and N2 depending on LOWER
*
      IF( LOWER ) THEN
         N2 = N / 2
         N1 = N - N2
      ELSE
         N1 = N / 2
         N2 = N - N1
      END IF
*
*     start execution: there are eight cases
*
      IF( NISODD ) THEN
*
*        N is odd
*
         IF( NORMALTRANSR ) THEN
*
*           N is odd and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
*             SRPA for LOWER, NORMAL and N is odd ( a(0:n-1,0:n1-1) )
*             T1 -> a(0,0), T2 -> a(0,1), S -> a(n1,0)
*             T1 -> a(0), T2 -> a(n), S -> a(n1)
*
               CALL DPOTRF( 'L', N1, A( 0 ), N, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRSM( 'R', 'L', 'T', 'N', N2, N1, ONE, A( 0 ), N, A( N1 ), N )                CALL DSYRK( 'U', 'N', N2, N1, -ONE, A( N1 ), N, ONE, A( N ), N )
               CALL DPOTRF( 'U', N2, A( N ), N, INFO )
               IF( INFO.GT.0 ) INFO = INFO + N1
*
            ELSE
*
*             SRPA for UPPER, NORMAL and N is odd ( a(0:n-1,0:n2-1)
*             T1 -> a(n1+1,0), T2 -> a(n1,0), S -> a(0,0)
*             T1 -> a(n2), T2 -> a(n1), S -> a(0)
*
               CALL DPOTRF( 'L', N1, A( N2 ), N, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRSM( 'L', 'L', 'N', 'N', N1, N2, ONE, A( N2 ), N, A( 0 ), N )                CALL DSYRK( 'U', 'T', N2, N1, -ONE, A( 0 ), N, ONE, A( N1 ), N )
               CALL DPOTRF( 'U', N2, A( N1 ), N, INFO )
               IF( INFO.GT.0 ) INFO = INFO + N1
*
            END IF
*
         ELSE
*
*           N is odd and TRANSR = 'T'
*
            IF( LOWER ) THEN
*
*              SRPA for LOWER, TRANSPOSE and N is odd
*              T1 -> A(0,0) , T2 -> A(1,0) , S -> A(0,n1)
*              T1 -> a(0+0) , T2 -> a(1+0) , S -> a(0+n1*n1); lda=n1
*
               CALL DPOTRF( 'U', N1, A( 0 ), N1, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRSM( 'L', 'U', 'T', 'N', N1, N2, ONE, A( 0 ), N1, A( N1*N1 ), N1 )                CALL DSYRK( 'L', 'T', N2, N1, -ONE, A( N1*N1 ), N1, ONE, A( 1 ), N1 )
               CALL DPOTRF( 'L', N2, A( 1 ), N1, INFO )
               IF( INFO.GT.0 ) INFO = INFO + N1
*
            ELSE
*
*              SRPA for UPPER, TRANSPOSE and N is odd
*              T1 -> A(0,n1+1), T2 -> A(0,n1), S -> A(0,0)
*              T1 -> a(n2*n2), T2 -> a(n1*n2), S -> a(0); lda = n2
*
               CALL DPOTRF( 'U', N1, A( N2*N2 ), N2, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRSM( 'R', 'U', 'N', 'N', N2, N1, ONE, A( N2*N2 ), N2, A( 0 ), N2 )                CALL DSYRK( 'L', 'N', N2, N1, -ONE, A( 0 ), N2, ONE, A( N1*N2 ), N2 )
               CALL DPOTRF( 'L', N2, A( N1*N2 ), N2, INFO )
               IF( INFO.GT.0 ) INFO = INFO + N1
*
            END IF
*
         END IF
*
      ELSE
*
*        N is even
*
         IF( NORMALTRANSR ) THEN
*
*           N is even and TRANSR = 'N'
*
            IF( LOWER ) THEN
*
*              SRPA for LOWER, NORMAL, and N is even ( a(0:n,0:k-1) )
*              T1 -> a(1,0), T2 -> a(0,0), S -> a(k+1,0)
*              T1 -> a(1), T2 -> a(0), S -> a(k+1)
*
               CALL DPOTRF( 'L', K, A( 1 ), N+1, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRSM( 'R', 'L', 'T', 'N', K, K, ONE, A( 1 ), N+1, A( K+1 ), N+1 )                CALL DSYRK( 'U', 'N', K, K, -ONE, A( K+1 ), N+1, ONE, A( 0 ), N+1 )
               CALL DPOTRF( 'U', K, A( 0 ), N+1, INFO )
               IF( INFO.GT.0 ) INFO = INFO + K
*
            ELSE
*
*              SRPA for UPPER, NORMAL, and N is even ( a(0:n,0:k-1) )
*              T1 -> a(k+1,0) ,  T2 -> a(k,0),   S -> a(0,0)
*              T1 -> a(k+1), T2 -> a(k), S -> a(0)
*
               CALL DPOTRF( 'L', K, A( K+1 ), N+1, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRSM( 'L', 'L', 'N', 'N', K, K, ONE, A( K+1 ), N+1, A( 0 ), N+1 )                CALL DSYRK( 'U', 'T', K, K, -ONE, A( 0 ), N+1, ONE, A( K ), N+1 )
               CALL DPOTRF( 'U', K, A( K ), N+1, INFO )
               IF( INFO.GT.0 ) INFO = INFO + K
*
            END IF
*
         ELSE
*
*           N is even and TRANSR = 'T'
*
            IF( LOWER ) THEN
*
*              SRPA for LOWER, TRANSPOSE and N is even (see paper)
*              T1 -> B(0,1), T2 -> B(0,0), S -> B(0,k+1)
*              T1 -> a(0+k), T2 -> a(0+0), S -> a(0+k*(k+1)); lda=k
*
               CALL DPOTRF( 'U', K, A( 0+K ), K, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRSM( 'L', 'U', 'T', 'N', K, K, ONE, A( K ), N1, A( K*( K+1 ) ), K )                CALL DSYRK( 'L', 'T', K, K, -ONE, A( K*( K+1 ) ), K, ONE, A( 0 ), K )
               CALL DPOTRF( 'L', K, A( 0 ), K, INFO )
               IF( INFO.GT.0 ) INFO = INFO + K
*
            ELSE
*
*              SRPA for UPPER, TRANSPOSE and N is even (see paper)
*              T1 -> B(0,k+1),     T2 -> B(0,k),   S -> B(0,0)
*              T1 -> a(0+k*(k+1)), T2 -> a(0+k*k), S -> a(0+0)); lda=k
*
               CALL DPOTRF( 'U', K, A( K*( K+1 ) ), K, INFO )
               IF( INFO.GT.0 ) RETURN                CALL DTRSM( 'R', 'U', 'N', 'N', K, K, ONE, A( K*( K+1 ) ), K, A( 0 ), K )                CALL DSYRK( 'L', 'N', K, K, -ONE, A( 0 ), K, ONE, A( K*K ), K )
               CALL DPOTRF( 'L', K, A( K*K ), K, INFO )
               IF( INFO.GT.0 ) INFO = INFO + K
*
            END IF
*
         END IF
*
      END IF
*
      RETURN
*
*     End of DPFTRF
*
      END
