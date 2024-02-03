      SUBROUTINE DSFRK( TRANSR, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      double             ALPHA, BETA;
      int                K, LDA, N;
      String             TRANS, TRANSR, UPLO;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), C( * );
      // ..
*
*  =====================================================================
*
      // ..
      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               LOWER, NORMALTRANSR, NISODD, NOTRANS;
      int                INFO, NROWA, J, NK, N1, N2;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DGEMM, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LOWER = LSAME( UPLO, 'L' )
      NOTRANS = LSAME( TRANS, 'N' )
*
      IF( NOTRANS ) THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
*
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'T' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOTRANS .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSFRK ', -INFO )
         RETURN
      END IF
*
      // Quick return if possible.
*
      // The quick return case: ((ALPHA.EQ.0).AND.(BETA.NE.ZERO)) is not
      // done (it is in DSYRK for example) and left in the general case.
*
      IF( ( N.EQ.0 ) .OR. ( ( ( ALPHA.EQ.ZERO ) .OR. ( K.EQ.0 ) ) .AND. ( BETA.EQ.ONE ) ) )RETURN
*
      IF( ( ALPHA.EQ.ZERO ) .AND. ( BETA.EQ.ZERO ) ) THEN
         DO J = 1, ( ( N*( N+1 ) ) / 2 )
            C( J ) = ZERO
         END DO
         RETURN
      END IF
*
      // C is N-by-N.
      // If N is odd, set NISODD = .TRUE., and N1 and N2.
      // If N is even, NISODD = .FALSE., and NK.
*
      IF( MOD( N, 2 ).EQ.0 ) THEN
         NISODD = .FALSE.
         NK = N / 2
      ELSE
         NISODD = .TRUE.
         IF( LOWER ) THEN
            N2 = N / 2
            N1 = N - N2
         ELSE
            N1 = N / 2
            N2 = N - N1
         END IF
      END IF
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
               // N is odd, TRANSR = 'N', and UPLO = 'L'
*
               IF( NOTRANS ) THEN
*
                  // N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'N'
*
                  CALL DSYRK( 'L', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N )                   CALL DSYRK( 'U', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( N+1 ), N )                   CALL DGEMM( 'N', 'T', N2, N1, K, ALPHA, A( N1+1, 1 ), LDA, A( 1, 1 ), LDA, BETA, C( N1+1 ), N )
*
               ELSE
*
                  // N is odd, TRANSR = 'N', UPLO = 'L', and TRANS = 'T'
*
                  CALL DSYRK( 'L', 'T', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N )                   CALL DSYRK( 'U', 'T', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( N+1 ), N )                   CALL DGEMM( 'T', 'N', N2, N1, K, ALPHA, A( 1, N1+1 ), LDA, A( 1, 1 ), LDA, BETA, C( N1+1 ), N )
*
               END IF
*
            ELSE
*
               // N is odd, TRANSR = 'N', and UPLO = 'U'
*
               IF( NOTRANS ) THEN
*
                  // N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'
*
                  CALL DSYRK( 'L', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2+1 ), N )                   CALL DSYRK( 'U', 'N', N2, K, ALPHA, A( N2, 1 ), LDA, BETA, C( N1+1 ), N )                   CALL DGEMM( 'N', 'T', N1, N2, K, ALPHA, A( 1, 1 ), LDA, A( N2, 1 ), LDA, BETA, C( 1 ), N )
*
               ELSE
*
                  // N is odd, TRANSR = 'N', UPLO = 'U', and TRANS = 'T'
*
                  CALL DSYRK( 'L', 'T', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2+1 ), N )                   CALL DSYRK( 'U', 'T', N2, K, ALPHA, A( 1, N2 ), LDA, BETA, C( N1+1 ), N )                   CALL DGEMM( 'T', 'N', N1, N2, K, ALPHA, A( 1, 1 ), LDA, A( 1, N2 ), LDA, BETA, C( 1 ), N )
*
               END IF
*
            END IF
*
         ELSE
*
            // N is odd, and TRANSR = 'T'
*
            IF( LOWER ) THEN
*
               // N is odd, TRANSR = 'T', and UPLO = 'L'
*
               IF( NOTRANS ) THEN
*
                  // N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'N'
*
                  CALL DSYRK( 'U', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N1 )                   CALL DSYRK( 'L', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( 2 ), N1 )                   CALL DGEMM( 'N', 'T', N1, N2, K, ALPHA, A( 1, 1 ), LDA, A( N1+1, 1 ), LDA, BETA, C( N1*N1+1 ), N1 )
*
               ELSE
*
                  // N is odd, TRANSR = 'T', UPLO = 'L', and TRANS = 'T'
*
                  CALL DSYRK( 'U', 'T', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 1 ), N1 )                   CALL DSYRK( 'L', 'T', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( 2 ), N1 )                   CALL DGEMM( 'T', 'N', N1, N2, K, ALPHA, A( 1, 1 ), LDA, A( 1, N1+1 ), LDA, BETA, C( N1*N1+1 ), N1 )
*
               END IF
*
            ELSE
*
               // N is odd, TRANSR = 'T', and UPLO = 'U'
*
               IF( NOTRANS ) THEN
*
                  // N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'N'
*
                  CALL DSYRK( 'U', 'N', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2*N2+1 ), N2 )                   CALL DSYRK( 'L', 'N', N2, K, ALPHA, A( N1+1, 1 ), LDA, BETA, C( N1*N2+1 ), N2 )                   CALL DGEMM( 'N', 'T', N2, N1, K, ALPHA, A( N1+1, 1 ), LDA, A( 1, 1 ), LDA, BETA, C( 1 ), N2 )
*
               ELSE
*
                  // N is odd, TRANSR = 'T', UPLO = 'U', and TRANS = 'T'
*
                  CALL DSYRK( 'U', 'T', N1, K, ALPHA, A( 1, 1 ), LDA, BETA, C( N2*N2+1 ), N2 )                   CALL DSYRK( 'L', 'T', N2, K, ALPHA, A( 1, N1+1 ), LDA, BETA, C( N1*N2+1 ), N2 )                   CALL DGEMM( 'T', 'N', N2, N1, K, ALPHA, A( 1, N1+1 ), LDA, A( 1, 1 ), LDA, BETA, C( 1 ), N2 )
*
               END IF
*
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
               // N is even, TRANSR = 'N', and UPLO = 'L'
*
               IF( NOTRANS ) THEN
*
                  // N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'N'
*
                  CALL DSYRK( 'L', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 2 ), N+1 )                   CALL DSYRK( 'U', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( 1 ), N+1 )                   CALL DGEMM( 'N', 'T', NK, NK, K, ALPHA, A( NK+1, 1 ), LDA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 )
*
               ELSE
*
                  // N is even, TRANSR = 'N', UPLO = 'L', and TRANS = 'T'
*
                  CALL DSYRK( 'L', 'T', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( 2 ), N+1 )                   CALL DSYRK( 'U', 'T', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( 1 ), N+1 )                   CALL DGEMM( 'T', 'N', NK, NK, K, ALPHA, A( 1, NK+1 ), LDA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 )
*
               END IF
*
            ELSE
*
               // N is even, TRANSR = 'N', and UPLO = 'U'
*
               IF( NOTRANS ) THEN
*
                  // N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'N'
*
                  CALL DSYRK( 'L', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 )                   CALL DSYRK( 'U', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( NK+1 ), N+1 )                   CALL DGEMM( 'N', 'T', NK, NK, K, ALPHA, A( 1, 1 ), LDA, A( NK+1, 1 ), LDA, BETA, C( 1 ), N+1 )
*
               ELSE
*
                  // N is even, TRANSR = 'N', UPLO = 'U', and TRANS = 'T'
*
                  CALL DSYRK( 'L', 'T', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+2 ), N+1 )                   CALL DSYRK( 'U', 'T', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( NK+1 ), N+1 )                   CALL DGEMM( 'T', 'N', NK, NK, K, ALPHA, A( 1, 1 ), LDA, A( 1, NK+1 ), LDA, BETA, C( 1 ), N+1 )
*
               END IF
*
            END IF
*
         ELSE
*
            // N is even, and TRANSR = 'T'
*
            IF( LOWER ) THEN
*
               // N is even, TRANSR = 'T', and UPLO = 'L'
*
               IF( NOTRANS ) THEN
*
                  // N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'N'
*
                  CALL DSYRK( 'U', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+1 ), NK )                   CALL DSYRK( 'L', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( 1 ), NK )                   CALL DGEMM( 'N', 'T', NK, NK, K, ALPHA, A( 1, 1 ), LDA, A( NK+1, 1 ), LDA, BETA, C( ( ( NK+1 )*NK )+1 ), NK )
*
               ELSE
*
                  // N is even, TRANSR = 'T', UPLO = 'L', and TRANS = 'T'
*
                  CALL DSYRK( 'U', 'T', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK+1 ), NK )                   CALL DSYRK( 'L', 'T', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( 1 ), NK )                   CALL DGEMM( 'T', 'N', NK, NK, K, ALPHA, A( 1, 1 ), LDA, A( 1, NK+1 ), LDA, BETA, C( ( ( NK+1 )*NK )+1 ), NK )
*
               END IF
*
            ELSE
*
               // N is even, TRANSR = 'T', and UPLO = 'U'
*
               IF( NOTRANS ) THEN
*
                  // N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'N'
*
                  CALL DSYRK( 'U', 'N', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK*( NK+1 )+1 ), NK )                   CALL DSYRK( 'L', 'N', NK, K, ALPHA, A( NK+1, 1 ), LDA, BETA, C( NK*NK+1 ), NK )                   CALL DGEMM( 'N', 'T', NK, NK, K, ALPHA, A( NK+1, 1 ), LDA, A( 1, 1 ), LDA, BETA, C( 1 ), NK )
*
               ELSE
*
                  // N is even, TRANSR = 'T', UPLO = 'U', and TRANS = 'T'
*
                  CALL DSYRK( 'U', 'T', NK, K, ALPHA, A( 1, 1 ), LDA, BETA, C( NK*( NK+1 )+1 ), NK )                   CALL DSYRK( 'L', 'T', NK, K, ALPHA, A( 1, NK+1 ), LDA, BETA, C( NK*NK+1 ), NK )                   CALL DGEMM( 'T', 'N', NK, NK, K, ALPHA, A( 1, NK+1 ), LDA, A( 1, 1 ), LDA, BETA, C( 1 ), NK )
*
               END IF
*
            END IF
*
         END IF
*
      END IF
*
      RETURN
*
      // End of DSFRK
*
      END
