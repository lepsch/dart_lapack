      SUBROUTINE SSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, HALF
      PARAMETER          ( ONE = 1.0, HALF = 0.5 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K, KB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSYGS2, SSYMM, SSYR2K, STRMM, STRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYGST', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'SSYGST', UPLO, N, -1, -1, -1 )
*
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL SSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      ELSE
*
*        Use blocked code
*
         IF( ITYPE.EQ.1 ) THEN
            IF( UPPER ) THEN
*
*              Compute inv(U**T)*A*inv(U)
*
               DO 10 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
*
*                 Update the upper triangle of A(k:n,k:n)
*
                  CALL SSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL STRSM( 'Left', UPLO, 'Transpose', 'Non-unit',
     $                           KB, N-K-KB+1, ONE, B( K, K ), LDB,
     $                           A( K, K+KB ), LDA )
                     CALL SSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF,
     $                           A( K, K ), LDA, B( K, K+KB ), LDB, ONE,
     $                           A( K, K+KB ), LDA )
                     CALL SSYR2K( UPLO, 'Transpose', N-K-KB+1, KB, -ONE,
     $                            A( K, K+KB ), LDA, B( K, K+KB ), LDB,
     $                            ONE, A( K+KB, K+KB ), LDA )
                     CALL SSYMM( 'Left', UPLO, KB, N-K-KB+1, -HALF,
     $                           A( K, K ), LDA, B( K, K+KB ), LDB, ONE,
     $                           A( K, K+KB ), LDA )
                     CALL STRSM( 'Right', UPLO, 'No transpose',
     $                           'Non-unit', KB, N-K-KB+1, ONE,
     $                           B( K+KB, K+KB ), LDB, A( K, K+KB ),
     $                           LDA )
                  END IF
   10          CONTINUE
            ELSE
*
*              Compute inv(L)*A*inv(L**T)
*
               DO 20 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
*
*                 Update the lower triangle of A(k:n,k:n)
*
                  CALL SSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL STRSM( 'Right', UPLO, 'Transpose', 'Non-unit',
     $                           N-K-KB+1, KB, ONE, B( K, K ), LDB,
     $                           A( K+KB, K ), LDA )
                     CALL SSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF,
     $                           A( K, K ), LDA, B( K+KB, K ), LDB, ONE,
     $                           A( K+KB, K ), LDA )
                     CALL SSYR2K( UPLO, 'No transpose', N-K-KB+1, KB,
     $                            -ONE, A( K+KB, K ), LDA, B( K+KB, K ),
     $                            LDB, ONE, A( K+KB, K+KB ), LDA )
                     CALL SSYMM( 'Right', UPLO, N-K-KB+1, KB, -HALF,
     $                           A( K, K ), LDA, B( K+KB, K ), LDB, ONE,
     $                           A( K+KB, K ), LDA )
                     CALL STRSM( 'Left', UPLO, 'No transpose',
     $                           'Non-unit', N-K-KB+1, KB, ONE,
     $                           B( K+KB, K+KB ), LDB, A( K+KB, K ),
     $                           LDA )
                  END IF
   20          CONTINUE
            END IF
         ELSE
            IF( UPPER ) THEN
*
*              Compute U*A*U**T
*
               DO 30 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
*
*                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
*
                  CALL STRMM( 'Left', UPLO, 'No transpose', 'Non-unit',
     $                        K-1, KB, ONE, B, LDB, A( 1, K ), LDA )
                  CALL SSYMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ),
     $                        LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA )
                  CALL SSYR2K( UPLO, 'No transpose', K-1, KB, ONE,
     $                         A( 1, K ), LDA, B( 1, K ), LDB, ONE, A,
     $                         LDA )
                  CALL SSYMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ),
     $                        LDA, B( 1, K ), LDB, ONE, A( 1, K ), LDA )
                  CALL STRMM( 'Right', UPLO, 'Transpose', 'Non-unit',
     $                        K-1, KB, ONE, B( K, K ), LDB, A( 1, K ),
     $                        LDA )
                  CALL SSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
   30          CONTINUE
            ELSE
*
*              Compute L**T*A*L
*
               DO 40 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
*
*                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
*
                  CALL STRMM( 'Right', UPLO, 'No transpose', 'Non-unit',
     $                        KB, K-1, ONE, B, LDB, A( K, 1 ), LDA )
                  CALL SSYMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ),
     $                        LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA )
                  CALL SSYR2K( UPLO, 'Transpose', K-1, KB, ONE,
     $                         A( K, 1 ), LDA, B( K, 1 ), LDB, ONE, A,
     $                         LDA )
                  CALL SSYMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ),
     $                        LDA, B( K, 1 ), LDB, ONE, A( K, 1 ), LDA )
                  CALL STRMM( 'Left', UPLO, 'Transpose', 'Non-unit', KB,
     $                        K-1, ONE, B( K, K ), LDB, A( K, 1 ), LDA )
                  CALL SSYGS2( ITYPE, UPLO, KB, A( K, K ), LDA,
     $                         B( K, K ), LDB, INFO )
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
*
*     End of SSYGST
*
      END