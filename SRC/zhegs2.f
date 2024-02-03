      SUBROUTINE ZHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, ITYPE, LDA, LDB, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ONE, HALF;
      PARAMETER          ( ONE = 1.0D+0, HALF = 0.5D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                K;
      double             AKK, BKK;
      COMPLEX*16         CT
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZAXPY, ZDSCAL, ZHER2, ZLACGV, ZTRMV, ZTRSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
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
         CALL XERBLA( 'ZHEGS2', -INFO )
         RETURN
      END IF
*
      IF( ITYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
*
            // Compute inv(U**H)*A*inv(U)
*
            DO 10 K = 1, N
*
               // Update the upper triangle of A(k:n,k:n)
*
               AKK = DBLE( A( K, K ) )
               BKK = DBLE( B( K, K ) )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL ZDSCAL( N-K, ONE / BKK, A( K, K+1 ), LDA )
                  CT = -HALF*AKK
                  CALL ZLACGV( N-K, A( K, K+1 ), LDA )
                  CALL ZLACGV( N-K, B( K, K+1 ), LDB )
                  CALL ZAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ), LDA )                   CALL ZHER2( UPLO, N-K, -CONE, A( K, K+1 ), LDA, B( K, K+1 ), LDB, A( K+1, K+1 ), LDA )                   CALL ZAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ), LDA )
                  CALL ZLACGV( N-K, B( K, K+1 ), LDB )
                  CALL ZTRSV( UPLO, 'Conjugate transpose', 'Non-unit', N-K, B( K+1, K+1 ), LDB, A( K, K+1 ), LDA )
                  CALL ZLACGV( N-K, A( K, K+1 ), LDA )
               END IF
   10       CONTINUE
         ELSE
*
            // Compute inv(L)*A*inv(L**H)
*
            DO 20 K = 1, N
*
               // Update the lower triangle of A(k:n,k:n)
*
               AKK = DBLE( A( K, K ) )
               BKK = DBLE( B( K, K ) )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL ZDSCAL( N-K, ONE / BKK, A( K+1, K ), 1 )
                  CT = -HALF*AKK
                  CALL ZAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL ZHER2( UPLO, N-K, -CONE, A( K+1, K ), 1, B( K+1, K ), 1, A( K+1, K+1 ), LDA )
                  CALL ZAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL ZTRSV( UPLO, 'No transpose', 'Non-unit', N-K, B( K+1, K+1 ), LDB, A( K+1, K ), 1 )
               END IF
   20       CONTINUE
         END IF
      ELSE
         IF( UPPER ) THEN
*
            // Compute U*A*U**H
*
            DO 30 K = 1, N
*
               // Update the upper triangle of A(1:k,1:k)
*
               AKK = DBLE( A( K, K ) )
               BKK = DBLE( B( K, K ) )
               CALL ZTRMV( UPLO, 'No transpose', 'Non-unit', K-1, B, LDB, A( 1, K ), 1 )
               CT = HALF*AKK
               CALL ZAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL ZHER2( UPLO, K-1, CONE, A( 1, K ), 1, B( 1, K ), 1, A, LDA )
               CALL ZAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL ZDSCAL( K-1, BKK, A( 1, K ), 1 )
               A( K, K ) = AKK*BKK**2
   30       CONTINUE
         ELSE
*
            // Compute L**H *A*L
*
            DO 40 K = 1, N
*
               // Update the lower triangle of A(1:k,1:k)
*
               AKK = DBLE( A( K, K ) )
               BKK = DBLE( B( K, K ) )
               CALL ZLACGV( K-1, A( K, 1 ), LDA )
               CALL ZTRMV( UPLO, 'Conjugate transpose', 'Non-unit', K-1, B, LDB, A( K, 1 ), LDA )
               CT = HALF*AKK
               CALL ZLACGV( K-1, B( K, 1 ), LDB )
               CALL ZAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL ZHER2( UPLO, K-1, CONE, A( K, 1 ), LDA, B( K, 1 ), LDB, A, LDA )
               CALL ZAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL ZLACGV( K-1, B( K, 1 ), LDB )
               CALL ZDSCAL( K-1, BKK, A( K, 1 ), LDA )
               CALL ZLACGV( K-1, A( K, 1 ), LDA )
               A( K, K ) = AKK*BKK**2
   40       CONTINUE
         END IF
      END IF
      RETURN
*
      // End of ZHEGS2
*
      END
