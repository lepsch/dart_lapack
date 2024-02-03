      SUBROUTINE SSYGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
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
      bool               UPPER;
      int                K;
      REAL               AKK, BKK, CT
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SSCAL, SSYR2, STRMV, STRSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX
*     ..
*     .. External Functions ..
      bool               LSAME;
      EXTERNAL           LSAME
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
         CALL XERBLA( 'SSYGS2', -INFO )
         RETURN
      END IF
*
      IF( ITYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
*
*           Compute inv(U**T)*A*inv(U)
*
            DO 10 K = 1, N
*
*              Update the upper triangle of A(k:n,k:n)
*
               AKK = A( K, K )
               BKK = B( K, K )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL SSCAL( N-K, ONE / BKK, A( K, K+1 ), LDA )
                  CT = -HALF*AKK
                  CALL SAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ), LDA )                   CALL SSYR2( UPLO, N-K, -ONE, A( K, K+1 ), LDA, B( K, K+1 ), LDB, A( K+1, K+1 ), LDA )                   CALL SAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ), LDA )                   CALL STRSV( UPLO, 'Transpose', 'Non-unit', N-K, B( K+1, K+1 ), LDB, A( K, K+1 ), LDA )
               END IF
   10       CONTINUE
         ELSE
*
*           Compute inv(L)*A*inv(L**T)
*
            DO 20 K = 1, N
*
*              Update the lower triangle of A(k:n,k:n)
*
               AKK = A( K, K )
               BKK = B( K, K )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL SSCAL( N-K, ONE / BKK, A( K+1, K ), 1 )
                  CT = -HALF*AKK
                  CALL SAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL SSYR2( UPLO, N-K, -ONE, A( K+1, K ), 1, B( K+1, K ), 1, A( K+1, K+1 ), LDA )
                  CALL SAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL STRSV( UPLO, 'No transpose', 'Non-unit', N-K, B( K+1, K+1 ), LDB, A( K+1, K ), 1 )
               END IF
   20       CONTINUE
         END IF
      ELSE
         IF( UPPER ) THEN
*
*           Compute U*A*U**T
*
            DO 30 K = 1, N
*
*              Update the upper triangle of A(1:k,1:k)
*
               AKK = A( K, K )
               BKK = B( K, K )
               CALL STRMV( UPLO, 'No transpose', 'Non-unit', K-1, B, LDB, A( 1, K ), 1 )
               CT = HALF*AKK
               CALL SAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL SSYR2( UPLO, K-1, ONE, A( 1, K ), 1, B( 1, K ), 1, A, LDA )
               CALL SAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL SSCAL( K-1, BKK, A( 1, K ), 1 )
               A( K, K ) = AKK*BKK**2
   30       CONTINUE
         ELSE
*
*           Compute L**T *A*L
*
            DO 40 K = 1, N
*
*              Update the lower triangle of A(1:k,1:k)
*
               AKK = A( K, K )
               BKK = B( K, K )
               CALL STRMV( UPLO, 'Transpose', 'Non-unit', K-1, B, LDB, A( K, 1 ), LDA )
               CT = HALF*AKK
               CALL SAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL SSYR2( UPLO, K-1, ONE, A( K, 1 ), LDA, B( K, 1 ), LDB, A, LDA )
               CALL SAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL SSCAL( K-1, BKK, A( K, 1 ), LDA )
               A( K, K ) = AKK*BKK**2
   40       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of SSYGS2
*
      END
