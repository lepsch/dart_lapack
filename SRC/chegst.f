      SUBROUTINE CHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, ITYPE, LDA, LDB, N;
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      COMPLEX            CONE, HALF
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ), HALF = ( 0.5E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      bool               UPPER;
      int                K, KB, NB;
*     ..
*     .. External Subroutines ..
      // EXTERNAL CHEGS2, CHEMM, CHER2K, CTRMM, CTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
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
         CALL XERBLA( 'CHEGST', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'CHEGST', UPLO, N, -1, -1, -1 )
*
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL CHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      ELSE
*
*        Use blocked code
*
         IF( ITYPE.EQ.1 ) THEN
            IF( UPPER ) THEN
*
*              Compute inv(U**H)*A*inv(U)
*
               DO 10 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
*
*                 Update the upper triangle of A(k:n,k:n)
*
                  CALL CHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL CTRSM( 'Left', UPLO, 'Conjugate transpose', 'Non-unit', KB, N-K-KB+1, CONE, B( K, K ), LDB, A( K, K+KB ), LDA )                      CALL CHEMM( 'Left', UPLO, KB, N-K-KB+1, -HALF, A( K, K ), LDA, B( K, K+KB ), LDB, CONE, A( K, K+KB ), LDA )                      CALL CHER2K( UPLO, 'Conjugate transpose', N-K-KB+1, KB, -CONE, A( K, K+KB ), LDA, B( K, K+KB ), LDB, ONE, A( K+KB, K+KB ), LDA )                      CALL CHEMM( 'Left', UPLO, KB, N-K-KB+1, -HALF, A( K, K ), LDA, B( K, K+KB ), LDB, CONE, A( K, K+KB ), LDA )                      CALL CTRSM( 'Right', UPLO, 'No transpose', 'Non-unit', KB, N-K-KB+1, CONE, B( K+KB, K+KB ), LDB, A( K, K+KB ), LDA )
                  END IF
   10          CONTINUE
            ELSE
*
*              Compute inv(L)*A*inv(L**H)
*
               DO 20 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
*
*                 Update the lower triangle of A(k:n,k:n)
*
                  CALL CHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO )
                  IF( K+KB.LE.N ) THEN
                     CALL CTRSM( 'Right', UPLO, 'Conjugate transpose', 'Non-unit', N-K-KB+1, KB, CONE, B( K, K ), LDB, A( K+KB, K ), LDA )                      CALL CHEMM( 'Right', UPLO, N-K-KB+1, KB, -HALF, A( K, K ), LDA, B( K+KB, K ), LDB, CONE, A( K+KB, K ), LDA )                      CALL CHER2K( UPLO, 'No transpose', N-K-KB+1, KB, -CONE, A( K+KB, K ), LDA, B( K+KB, K ), LDB, ONE, A( K+KB, K+KB ), LDA )                      CALL CHEMM( 'Right', UPLO, N-K-KB+1, KB, -HALF, A( K, K ), LDA, B( K+KB, K ), LDB, CONE, A( K+KB, K ), LDA )                      CALL CTRSM( 'Left', UPLO, 'No transpose', 'Non-unit', N-K-KB+1, KB, CONE, B( K+KB, K+KB ), LDB, A( K+KB, K ), LDA )
                  END IF
   20          CONTINUE
            END IF
         ELSE
            IF( UPPER ) THEN
*
*              Compute U*A*U**H
*
               DO 30 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
*
*                 Update the upper triangle of A(1:k+kb-1,1:k+kb-1)
*
                  CALL CTRMM( 'Left', UPLO, 'No transpose', 'Non-unit', K-1, KB, CONE, B, LDB, A( 1, K ), LDA )                   CALL CHEMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ), LDA, B( 1, K ), LDB, CONE, A( 1, K ), LDA )                   CALL CHER2K( UPLO, 'No transpose', K-1, KB, CONE, A( 1, K ), LDA, B( 1, K ), LDB, ONE, A, LDA )                   CALL CHEMM( 'Right', UPLO, K-1, KB, HALF, A( K, K ), LDA, B( 1, K ), LDB, CONE, A( 1, K ), LDA )                   CALL CTRMM( 'Right', UPLO, 'Conjugate transpose', 'Non-unit', K-1, KB, CONE, B( K, K ), LDB, A( 1, K ), LDA )
                  CALL CHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO )
   30          CONTINUE
            ELSE
*
*              Compute L**H*A*L
*
               DO 40 K = 1, N, NB
                  KB = MIN( N-K+1, NB )
*
*                 Update the lower triangle of A(1:k+kb-1,1:k+kb-1)
*
                  CALL CTRMM( 'Right', UPLO, 'No transpose', 'Non-unit', KB, K-1, CONE, B, LDB, A( K, 1 ), LDA )                   CALL CHEMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ), LDA, B( K, 1 ), LDB, CONE, A( K, 1 ), LDA )                   CALL CHER2K( UPLO, 'Conjugate transpose', K-1, KB, CONE, A( K, 1 ), LDA, B( K, 1 ), LDB, ONE, A, LDA )                   CALL CHEMM( 'Left', UPLO, KB, K-1, HALF, A( K, K ), LDA, B( K, 1 ), LDB, CONE, A( K, 1 ), LDA )                   CALL CTRMM( 'Left', UPLO, 'Conjugate transpose', 'Non-unit', KB, K-1, CONE, B( K, K ), LDB, A( K, 1 ), LDA )
                  CALL CHEGS2( ITYPE, UPLO, KB, A( K, K ), LDA, B( K, K ), LDB, INFO )
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
*
*     End of CHEGST
*
      END
