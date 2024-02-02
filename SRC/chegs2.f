      SUBROUTINE CHEGS2( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
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
      COMPLEX            A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, HALF
      PARAMETER          ( ONE = 1.0E+0, HALF = 0.5E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            K
      REAL               AKK, BKK
      COMPLEX            CT
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CHER2, CLACGV, CSSCAL, CTRMV, CTRSV,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
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
         CALL XERBLA( 'CHEGS2', -INFO )
         RETURN
      END IF
*
      IF( ITYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
*
*           Compute inv(U**H)*A*inv(U)
*
            DO 10 K = 1, N
*
*              Update the upper triangle of A(k:n,k:n)
*
               AKK = REAL( A( K, K ) )
               BKK = REAL( B( K, K ) )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL CSSCAL( N-K, ONE / BKK, A( K, K+1 ), LDA )
                  CT = -HALF*AKK
                  CALL CLACGV( N-K, A( K, K+1 ), LDA )
                  CALL CLACGV( N-K, B( K, K+1 ), LDB )
                  CALL CAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ),
     $                        LDA )
                  CALL CHER2( UPLO, N-K, -CONE, A( K, K+1 ), LDA,
     $                        B( K, K+1 ), LDB, A( K+1, K+1 ), LDA )
                  CALL CAXPY( N-K, CT, B( K, K+1 ), LDB, A( K, K+1 ),
     $                        LDA )
                  CALL CLACGV( N-K, B( K, K+1 ), LDB )
                  CALL CTRSV( UPLO, 'Conjugate transpose', 'Non-unit',
     $                        N-K, B( K+1, K+1 ), LDB, A( K, K+1 ),
     $                        LDA )
                  CALL CLACGV( N-K, A( K, K+1 ), LDA )
               END IF
   10       CONTINUE
         ELSE
*
*           Compute inv(L)*A*inv(L**H)
*
            DO 20 K = 1, N
*
*              Update the lower triangle of A(k:n,k:n)
*
               AKK = REAL( A( K, K ) )
               BKK = REAL( B( K, K ) )
               AKK = AKK / BKK**2
               A( K, K ) = AKK
               IF( K.LT.N ) THEN
                  CALL CSSCAL( N-K, ONE / BKK, A( K+1, K ), 1 )
                  CT = -HALF*AKK
                  CALL CAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL CHER2( UPLO, N-K, -CONE, A( K+1, K ), 1,
     $                        B( K+1, K ), 1, A( K+1, K+1 ), LDA )
                  CALL CAXPY( N-K, CT, B( K+1, K ), 1, A( K+1, K ), 1 )
                  CALL CTRSV( UPLO, 'No transpose', 'Non-unit', N-K,
     $                        B( K+1, K+1 ), LDB, A( K+1, K ), 1 )
               END IF
   20       CONTINUE
         END IF
      ELSE
         IF( UPPER ) THEN
*
*           Compute U*A*U**H
*
            DO 30 K = 1, N
*
*              Update the upper triangle of A(1:k,1:k)
*
               AKK = REAL( A( K, K ) )
               BKK = REAL( B( K, K ) )
               CALL CTRMV( UPLO, 'No transpose', 'Non-unit', K-1, B,
     $                     LDB, A( 1, K ), 1 )
               CT = HALF*AKK
               CALL CAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL CHER2( UPLO, K-1, CONE, A( 1, K ), 1, B( 1, K ), 1,
     $                     A, LDA )
               CALL CAXPY( K-1, CT, B( 1, K ), 1, A( 1, K ), 1 )
               CALL CSSCAL( K-1, BKK, A( 1, K ), 1 )
               A( K, K ) = AKK*BKK**2
   30       CONTINUE
         ELSE
*
*           Compute L**H *A*L
*
            DO 40 K = 1, N
*
*              Update the lower triangle of A(1:k,1:k)
*
               AKK = REAL( A( K, K ) )
               BKK = REAL( B( K, K ) )
               CALL CLACGV( K-1, A( K, 1 ), LDA )
               CALL CTRMV( UPLO, 'Conjugate transpose', 'Non-unit', K-1,
     $                     B, LDB, A( K, 1 ), LDA )
               CT = HALF*AKK
               CALL CLACGV( K-1, B( K, 1 ), LDB )
               CALL CAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL CHER2( UPLO, K-1, CONE, A( K, 1 ), LDA, B( K, 1 ),
     $                     LDB, A, LDA )
               CALL CAXPY( K-1, CT, B( K, 1 ), LDB, A( K, 1 ), LDA )
               CALL CLACGV( K-1, B( K, 1 ), LDB )
               CALL CSSCAL( K-1, BKK, A( K, 1 ), LDA )
               CALL CLACGV( K-1, A( K, 1 ), LDA )
               A( K, K ) = AKK*BKK**2
   40       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of CHEGS2
*
      END