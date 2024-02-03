      SUBROUTINE DTFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, B, LDB )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             TRANSR, DIAG, SIDE, TRANS, UPLO;
      int                LDB, M, N;
      double             ALPHA;
      // ..
      // .. Array Arguments ..
      double             A( 0: * ), B( 0: LDB-1, 0: * );
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
      bool               LOWER, LSIDE, MISODD, NISODD, NORMALTRANSR, NOTRANS;
      int                M1, M2, N1, N2, K, INFO, I, J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DGEMM, DTRSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MOD
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LSIDE = LSAME( SIDE, 'L' )
      LOWER = LSAME( UPLO, 'L' )
      NOTRANS = LSAME( TRANS, 'N' )
      IF( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'T' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSIDE .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.NOTRANS .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -4
      ELSE IF( .NOT.LSAME( DIAG, 'N' ) .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, M ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTFSM ', -INFO )
         RETURN
      END IF
*
      // Quick return when ( (N.EQ.0).OR.(M.EQ.0) )
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) RETURN
*
      // Quick return when ALPHA.EQ.(0D+0)
*
      IF( ALPHA.EQ.ZERO ) THEN
         DO 20 J = 0, N - 1
            DO 10 I = 0, M - 1
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
      IF( LSIDE ) THEN
*
         // SIDE = 'L'
*
         // A is M-by-M.
         // If M is odd, set NISODD = .TRUE., and M1 and M2.
         // If M is even, NISODD = .FALSE., and M.
*
         IF( MOD( M, 2 ).EQ.0 ) THEN
            MISODD = .FALSE.
            K = M / 2
         ELSE
            MISODD = .TRUE.
            IF( LOWER ) THEN
               M2 = M / 2
               M1 = M - M2
            ELSE
               M1 = M / 2
               M2 = M - M1
            END IF
         END IF
*
*
         IF( MISODD ) THEN
*
            // SIDE = 'L' and N is odd
*
            IF( NORMALTRANSR ) THEN
*
               // SIDE = 'L', N is odd, and TRANSR = 'N'
*
               IF( LOWER ) THEN
*
                  // SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
                     // TRANS = 'N'
*
                     IF( M.EQ.1 ) THEN
                        CALL DTRSM( 'L', 'L', 'N', DIAG, M1, N, ALPHA, A, M, B, LDB )
                     ELSE
                        CALL DTRSM( 'L', 'L', 'N', DIAG, M1, N, ALPHA, A( 0 ), M, B, LDB )                         CALL DGEMM( 'N', 'N', M2, N, M1, -ONE, A( M1 ), M, B, LDB, ALPHA, B( M1, 0 ), LDB )                         CALL DTRSM( 'L', 'U', 'T', DIAG, M2, N, ONE, A( M ), M, B( M1, 0 ), LDB )
                     END IF
*
                  ELSE
*
                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
                     // TRANS = 'T'
*
                     IF( M.EQ.1 ) THEN
                        CALL DTRSM( 'L', 'L', 'T', DIAG, M1, N, ALPHA, A( 0 ), M, B, LDB )
                     ELSE
                        CALL DTRSM( 'L', 'U', 'N', DIAG, M2, N, ALPHA, A( M ), M, B( M1, 0 ), LDB )                         CALL DGEMM( 'T', 'N', M1, N, M2, -ONE, A( M1 ), M, B( M1, 0 ), LDB, ALPHA, B, LDB )                         CALL DTRSM( 'L', 'L', 'T', DIAG, M1, N, ONE, A( 0 ), M, B, LDB )
                     END IF
*
                  END IF
*
               ELSE
*
                  // SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U'
*
                  IF( .NOT.NOTRANS ) THEN
*
                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'N'
*
                     CALL DTRSM( 'L', 'L', 'N', DIAG, M1, N, ALPHA, A( M2 ), M, B, LDB )                      CALL DGEMM( 'T', 'N', M2, N, M1, -ONE, A( 0 ), M, B, LDB, ALPHA, B( M1, 0 ), LDB )                      CALL DTRSM( 'L', 'U', 'T', DIAG, M2, N, ONE, A( M1 ), M, B( M1, 0 ), LDB )
*
                  ELSE
*
                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'T'
*
                     CALL DTRSM( 'L', 'U', 'N', DIAG, M2, N, ALPHA, A( M1 ), M, B( M1, 0 ), LDB )                      CALL DGEMM( 'N', 'N', M1, N, M2, -ONE, A( 0 ), M, B( M1, 0 ), LDB, ALPHA, B, LDB )                      CALL DTRSM( 'L', 'L', 'T', DIAG, M1, N, ONE, A( M2 ), M, B, LDB )
*
                  END IF
*
               END IF
*
            ELSE
*
               // SIDE = 'L', N is odd, and TRANSR = 'T'
*
               IF( LOWER ) THEN
*
                  // SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
                     // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and
                     // TRANS = 'N'
*
                     IF( M.EQ.1 ) THEN
                        CALL DTRSM( 'L', 'U', 'T', DIAG, M1, N, ALPHA, A( 0 ), M1, B, LDB )
                     ELSE
                        CALL DTRSM( 'L', 'U', 'T', DIAG, M1, N, ALPHA, A( 0 ), M1, B, LDB )                         CALL DGEMM( 'T', 'N', M2, N, M1, -ONE, A( M1*M1 ), M1, B, LDB, ALPHA, B( M1, 0 ), LDB )
                        CALL DTRSM( 'L', 'L', 'N', DIAG, M2, N, ONE, A( 1 ), M1, B( M1, 0 ), LDB )
                     END IF
*
                  ELSE
*
                     // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and
                     // TRANS = 'T'
*
                     IF( M.EQ.1 ) THEN
                        CALL DTRSM( 'L', 'U', 'N', DIAG, M1, N, ALPHA, A( 0 ), M1, B, LDB )
                     ELSE
                        CALL DTRSM( 'L', 'L', 'T', DIAG, M2, N, ALPHA, A( 1 ), M1, B( M1, 0 ), LDB )                         CALL DGEMM( 'N', 'N', M1, N, M2, -ONE, A( M1*M1 ), M1, B( M1, 0 ), LDB, ALPHA, B, LDB )
                        CALL DTRSM( 'L', 'U', 'N', DIAG, M1, N, ONE, A( 0 ), M1, B, LDB )
                     END IF
*
                  END IF
*
               ELSE
*
                  // SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'U'
*
                  IF( .NOT.NOTRANS ) THEN
*
                     // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and
                     // TRANS = 'N'
*
                     CALL DTRSM( 'L', 'U', 'T', DIAG, M1, N, ALPHA, A( M2*M2 ), M2, B, LDB )                      CALL DGEMM( 'N', 'N', M2, N, M1, -ONE, A( 0 ), M2, B, LDB, ALPHA, B( M1, 0 ), LDB )                      CALL DTRSM( 'L', 'L', 'N', DIAG, M2, N, ONE, A( M1*M2 ), M2, B( M1, 0 ), LDB )
*
                  ELSE
*
                     // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and
                     // TRANS = 'T'
*
                     CALL DTRSM( 'L', 'L', 'T', DIAG, M2, N, ALPHA, A( M1*M2 ), M2, B( M1, 0 ), LDB )                      CALL DGEMM( 'T', 'N', M1, N, M2, -ONE, A( 0 ), M2, B( M1, 0 ), LDB, ALPHA, B, LDB )                      CALL DTRSM( 'L', 'U', 'N', DIAG, M1, N, ONE, A( M2*M2 ), M2, B, LDB )
*
                  END IF
*
               END IF
*
            END IF
*
         ELSE
*
            // SIDE = 'L' and N is even
*
            IF( NORMALTRANSR ) THEN
*
               // SIDE = 'L', N is even, and TRANSR = 'N'
*
               IF( LOWER ) THEN
*
                  // SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
                     // and TRANS = 'N'
*
                     CALL DTRSM( 'L', 'L', 'N', DIAG, K, N, ALPHA, A( 1 ), M+1, B, LDB )                      CALL DGEMM( 'N', 'N', K, N, K, -ONE, A( K+1 ), M+1, B, LDB, ALPHA, B( K, 0 ), LDB )                      CALL DTRSM( 'L', 'U', 'T', DIAG, K, N, ONE, A( 0 ), M+1, B( K, 0 ), LDB )
*
                  ELSE
*
                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
                     // and TRANS = 'T'
*
                     CALL DTRSM( 'L', 'U', 'N', DIAG, K, N, ALPHA, A( 0 ), M+1, B( K, 0 ), LDB )                      CALL DGEMM( 'T', 'N', K, N, K, -ONE, A( K+1 ), M+1, B( K, 0 ), LDB, ALPHA, B, LDB )                      CALL DTRSM( 'L', 'L', 'T', DIAG, K, N, ONE, A( 1 ), M+1, B, LDB )
*
                  END IF
*
               ELSE
*
                  // SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U'
*
                  IF( .NOT.NOTRANS ) THEN
*
                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'N'
*
                     CALL DTRSM( 'L', 'L', 'N', DIAG, K, N, ALPHA, A( K+1 ), M+1, B, LDB )                      CALL DGEMM( 'T', 'N', K, N, K, -ONE, A( 0 ), M+1, B, LDB, ALPHA, B( K, 0 ), LDB )                      CALL DTRSM( 'L', 'U', 'T', DIAG, K, N, ONE, A( K ), M+1, B( K, 0 ), LDB )
*
                  ELSE
*
                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'T'
                     CALL DTRSM( 'L', 'U', 'N', DIAG, K, N, ALPHA, A( K ), M+1, B( K, 0 ), LDB )                      CALL DGEMM( 'N', 'N', K, N, K, -ONE, A( 0 ), M+1, B( K, 0 ), LDB, ALPHA, B, LDB )                      CALL DTRSM( 'L', 'L', 'T', DIAG, K, N, ONE, A( K+1 ), M+1, B, LDB )
*
                  END IF
*
               END IF
*
            ELSE
*
               // SIDE = 'L', N is even, and TRANSR = 'T'
*
               IF( LOWER ) THEN
*
                  // SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
                     // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L',
                     // and TRANS = 'N'
*
                     CALL DTRSM( 'L', 'U', 'T', DIAG, K, N, ALPHA, A( K ), K, B, LDB )                      CALL DGEMM( 'T', 'N', K, N, K, -ONE, A( K*( K+1 ) ), K, B, LDB, ALPHA, B( K, 0 ), LDB )
                     CALL DTRSM( 'L', 'L', 'N', DIAG, K, N, ONE, A( 0 ), K, B( K, 0 ), LDB )
*
                  ELSE
*
                     // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L',
                     // and TRANS = 'T'
*
                     CALL DTRSM( 'L', 'L', 'T', DIAG, K, N, ALPHA, A( 0 ), K, B( K, 0 ), LDB )                      CALL DGEMM( 'N', 'N', K, N, K, -ONE, A( K*( K+1 ) ), K, B( K, 0 ), LDB, ALPHA, B, LDB )
                     CALL DTRSM( 'L', 'U', 'N', DIAG, K, N, ONE, A( K ), K, B, LDB )
*
                  END IF
*
               ELSE
*
                  // SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'U'
*
                  IF( .NOT.NOTRANS ) THEN
*
                     // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U',
                     // and TRANS = 'N'
*
                     CALL DTRSM( 'L', 'U', 'T', DIAG, K, N, ALPHA, A( K*( K+1 ) ), K, B, LDB )                      CALL DGEMM( 'N', 'N', K, N, K, -ONE, A( 0 ), K, B, LDB, ALPHA, B( K, 0 ), LDB )                      CALL DTRSM( 'L', 'L', 'N', DIAG, K, N, ONE, A( K*K ), K, B( K, 0 ), LDB )
*
                  ELSE
*
                     // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U',
                     // and TRANS = 'T'
*
                     CALL DTRSM( 'L', 'L', 'T', DIAG, K, N, ALPHA, A( K*K ), K, B( K, 0 ), LDB )                      CALL DGEMM( 'T', 'N', K, N, K, -ONE, A( 0 ), K, B( K, 0 ), LDB, ALPHA, B, LDB )                      CALL DTRSM( 'L', 'U', 'N', DIAG, K, N, ONE, A( K*( K+1 ) ), K, B, LDB )
*
                  END IF
*
               END IF
*
            END IF
*
         END IF
*
      ELSE
*
         // SIDE = 'R'
*
         // A is N-by-N.
         // If N is odd, set NISODD = .TRUE., and N1 and N2.
         // If N is even, NISODD = .FALSE., and K.
*
         IF( MOD( N, 2 ).EQ.0 ) THEN
            NISODD = .FALSE.
            K = N / 2
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
            // SIDE = 'R' and N is odd
*
            IF( NORMALTRANSR ) THEN
*
               // SIDE = 'R', N is odd, and TRANSR = 'N'
*
               IF( LOWER ) THEN
*
                  // SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
                     // TRANS = 'N'
*
                     CALL DTRSM( 'R', 'U', 'T', DIAG, M, N2, ALPHA, A( N ), N, B( 0, N1 ), LDB )                      CALL DGEMM( 'N', 'N', M, N1, N2, -ONE, B( 0, N1 ), LDB, A( N1 ), N, ALPHA, B( 0, 0 ), LDB )
                     CALL DTRSM( 'R', 'L', 'N', DIAG, M, N1, ONE, A( 0 ), N, B( 0, 0 ), LDB )
*
                  ELSE
*
                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
                     // TRANS = 'T'
*
                     CALL DTRSM( 'R', 'L', 'T', DIAG, M, N1, ALPHA, A( 0 ), N, B( 0, 0 ), LDB )                      CALL DGEMM( 'N', 'T', M, N2, N1, -ONE, B( 0, 0 ), LDB, A( N1 ), N, ALPHA, B( 0, N1 ), LDB )
                     CALL DTRSM( 'R', 'U', 'N', DIAG, M, N2, ONE, A( N ), N, B( 0, N1 ), LDB )
*
                  END IF
*
               ELSE
*
                  // SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U'
*
                  IF( NOTRANS ) THEN
*
                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'N'
*
                     CALL DTRSM( 'R', 'L', 'T', DIAG, M, N1, ALPHA, A( N2 ), N, B( 0, 0 ), LDB )                      CALL DGEMM( 'N', 'N', M, N2, N1, -ONE, B( 0, 0 ), LDB, A( 0 ), N, ALPHA, B( 0, N1 ), LDB )
                     CALL DTRSM( 'R', 'U', 'N', DIAG, M, N2, ONE, A( N1 ), N, B( 0, N1 ), LDB )
*
                  ELSE
*
                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'T'
*
                     CALL DTRSM( 'R', 'U', 'T', DIAG, M, N2, ALPHA, A( N1 ), N, B( 0, N1 ), LDB )                      CALL DGEMM( 'N', 'T', M, N1, N2, -ONE, B( 0, N1 ), LDB, A( 0 ), N, ALPHA, B( 0, 0 ), LDB )                      CALL DTRSM( 'R', 'L', 'N', DIAG, M, N1, ONE, A( N2 ), N, B( 0, 0 ), LDB )
*
                  END IF
*
               END IF
*
            ELSE
*
               // SIDE = 'R', N is odd, and TRANSR = 'T'
*
               IF( LOWER ) THEN
*
                  // SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
                     // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and
                     // TRANS = 'N'
*
                     CALL DTRSM( 'R', 'L', 'N', DIAG, M, N2, ALPHA, A( 1 ), N1, B( 0, N1 ), LDB )                      CALL DGEMM( 'N', 'T', M, N1, N2, -ONE, B( 0, N1 ), LDB, A( N1*N1 ), N1, ALPHA, B( 0, 0 ), LDB )
                     CALL DTRSM( 'R', 'U', 'T', DIAG, M, N1, ONE, A( 0 ), N1, B( 0, 0 ), LDB )
*
                  ELSE
*
                     // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and
                     // TRANS = 'T'
*
                     CALL DTRSM( 'R', 'U', 'N', DIAG, M, N1, ALPHA, A( 0 ), N1, B( 0, 0 ), LDB )                      CALL DGEMM( 'N', 'N', M, N2, N1, -ONE, B( 0, 0 ), LDB, A( N1*N1 ), N1, ALPHA, B( 0, N1 ), LDB )
                     CALL DTRSM( 'R', 'L', 'T', DIAG, M, N2, ONE, A( 1 ), N1, B( 0, N1 ), LDB )
*
                  END IF
*
               ELSE
*
                  // SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'U'
*
                  IF( NOTRANS ) THEN
*
                     // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and
                     // TRANS = 'N'
*
                     CALL DTRSM( 'R', 'U', 'N', DIAG, M, N1, ALPHA, A( N2*N2 ), N2, B( 0, 0 ), LDB )                      CALL DGEMM( 'N', 'T', M, N2, N1, -ONE, B( 0, 0 ), LDB, A( 0 ), N2, ALPHA, B( 0, N1 ), LDB )
                     CALL DTRSM( 'R', 'L', 'T', DIAG, M, N2, ONE, A( N1*N2 ), N2, B( 0, N1 ), LDB )
*
                  ELSE
*
                     // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and
                     // TRANS = 'T'
*
                     CALL DTRSM( 'R', 'L', 'N', DIAG, M, N2, ALPHA, A( N1*N2 ), N2, B( 0, N1 ), LDB )                      CALL DGEMM( 'N', 'N', M, N1, N2, -ONE, B( 0, N1 ), LDB, A( 0 ), N2, ALPHA, B( 0, 0 ), LDB )
                     CALL DTRSM( 'R', 'U', 'T', DIAG, M, N1, ONE, A( N2*N2 ), N2, B( 0, 0 ), LDB )
*
                  END IF
*
               END IF
*
            END IF
*
         ELSE
*
            // SIDE = 'R' and N is even
*
            IF( NORMALTRANSR ) THEN
*
               // SIDE = 'R', N is even, and TRANSR = 'N'
*
               IF( LOWER ) THEN
*
                  // SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
                     // and TRANS = 'N'
*
                     CALL DTRSM( 'R', 'U', 'T', DIAG, M, K, ALPHA, A( 0 ), N+1, B( 0, K ), LDB )                      CALL DGEMM( 'N', 'N', M, K, K, -ONE, B( 0, K ), LDB, A( K+1 ), N+1, ALPHA, B( 0, 0 ), LDB )
                     CALL DTRSM( 'R', 'L', 'N', DIAG, M, K, ONE, A( 1 ), N+1, B( 0, 0 ), LDB )
*
                  ELSE
*
                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
                     // and TRANS = 'T'
*
                     CALL DTRSM( 'R', 'L', 'T', DIAG, M, K, ALPHA, A( 1 ), N+1, B( 0, 0 ), LDB )                      CALL DGEMM( 'N', 'T', M, K, K, -ONE, B( 0, 0 ), LDB, A( K+1 ), N+1, ALPHA, B( 0, K ), LDB )
                     CALL DTRSM( 'R', 'U', 'N', DIAG, M, K, ONE, A( 0 ), N+1, B( 0, K ), LDB )
*
                  END IF
*
               ELSE
*
                  // SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U'
*
                  IF( NOTRANS ) THEN
*
                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'N'
*
                     CALL DTRSM( 'R', 'L', 'T', DIAG, M, K, ALPHA, A( K+1 ), N+1, B( 0, 0 ), LDB )                      CALL DGEMM( 'N', 'N', M, K, K, -ONE, B( 0, 0 ), LDB, A( 0 ), N+1, ALPHA, B( 0, K ), LDB )
                     CALL DTRSM( 'R', 'U', 'N', DIAG, M, K, ONE, A( K ), N+1, B( 0, K ), LDB )
*
                  ELSE
*
                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'T'
*
                     CALL DTRSM( 'R', 'U', 'T', DIAG, M, K, ALPHA, A( K ), N+1, B( 0, K ), LDB )                      CALL DGEMM( 'N', 'T', M, K, K, -ONE, B( 0, K ), LDB, A( 0 ), N+1, ALPHA, B( 0, 0 ), LDB )
                     CALL DTRSM( 'R', 'L', 'N', DIAG, M, K, ONE, A( K+1 ), N+1, B( 0, 0 ), LDB )
*
                  END IF
*
               END IF
*
            ELSE
*
               // SIDE = 'R', N is even, and TRANSR = 'T'
*
               IF( LOWER ) THEN
*
                  // SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'L'
*
                  IF( NOTRANS ) THEN
*
                     // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L',
                     // and TRANS = 'N'
*
                     CALL DTRSM( 'R', 'L', 'N', DIAG, M, K, ALPHA, A( 0 ), K, B( 0, K ), LDB )                      CALL DGEMM( 'N', 'T', M, K, K, -ONE, B( 0, K ), LDB, A( ( K+1 )*K ), K, ALPHA, B( 0, 0 ), LDB )
                     CALL DTRSM( 'R', 'U', 'T', DIAG, M, K, ONE, A( K ), K, B( 0, 0 ), LDB )
*
                  ELSE
*
                     // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L',
                     // and TRANS = 'T'
*
                     CALL DTRSM( 'R', 'U', 'N', DIAG, M, K, ALPHA, A( K ), K, B( 0, 0 ), LDB )                      CALL DGEMM( 'N', 'N', M, K, K, -ONE, B( 0, 0 ), LDB, A( ( K+1 )*K ), K, ALPHA, B( 0, K ), LDB )
                     CALL DTRSM( 'R', 'L', 'T', DIAG, M, K, ONE, A( 0 ), K, B( 0, K ), LDB )
*
                  END IF
*
               ELSE
*
                  // SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'U'
*
                  IF( NOTRANS ) THEN
*
                     // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U',
                     // and TRANS = 'N'
*
                     CALL DTRSM( 'R', 'U', 'N', DIAG, M, K, ALPHA, A( ( K+1 )*K ), K, B( 0, 0 ), LDB )                      CALL DGEMM( 'N', 'T', M, K, K, -ONE, B( 0, 0 ), LDB, A( 0 ), K, ALPHA, B( 0, K ), LDB )                      CALL DTRSM( 'R', 'L', 'T', DIAG, M, K, ONE, A( K*K ), K, B( 0, K ), LDB )
*
                  ELSE
*
                     // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U',
                     // and TRANS = 'T'
*
                     CALL DTRSM( 'R', 'L', 'N', DIAG, M, K, ALPHA, A( K*K ), K, B( 0, K ), LDB )                      CALL DGEMM( 'N', 'N', M, K, K, -ONE, B( 0, K ), LDB, A( 0 ), K, ALPHA, B( 0, 0 ), LDB )                      CALL DTRSM( 'R', 'U', 'T', DIAG, M, K, ONE, A( ( K+1 )*K ), K, B( 0, 0 ), LDB )
*
                  END IF
*
               END IF
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
      // End of DTFSM
*
      END
