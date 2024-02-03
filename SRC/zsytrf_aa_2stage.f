      SUBROUTINE ZSYTRF_AA_2STAGE( UPLO, N, A, LDA, TB, LTB, IPIV, IPIV2, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                N, LDA, LTB, LWORK, INFO
*     ..
*     .. Array Arguments ..
      int                IPIV( * ), IPIV2( * )
      COMPLEX*16         A( LDA, * ), TB( * ), WORK( * )
*     ..
*
*  =====================================================================
*     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE  = ( 1.0D+0, 0.0D+0 ) )
*
*     .. Local Scalars ..
      LOGICAL            UPPER, TQUERY, WQUERY
      int                I, J, K, I1, I2, TD
      int                LDTB, NB, KB, JB, NT, IINFO
      COMPLEX*16         PIV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      int                ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZCOPY, ZGBTRF, ZGEMM, ZGETRF,   ZLACPY, ZLASET, ZLASWP, ZTRSM, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      WQUERY = ( LWORK.EQ.-1 )
      TQUERY = ( LTB.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF ( LTB .LT. 4*N .AND. .NOT.TQUERY ) THEN
         INFO = -6
      ELSE IF ( LWORK .LT. N .AND. .NOT.WQUERY ) THEN
         INFO = -10
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSYTRF_AA_2STAGE', -INFO )
         RETURN
      END IF
*
*     Answer the query
*
      NB = ILAENV( 1, 'ZSYTRF_AA_2STAGE', UPLO, N, -1, -1, -1 )
      IF( INFO.EQ.0 ) THEN
         IF( TQUERY ) THEN
            TB( 1 ) = (3*NB+1)*N
         END IF
         IF( WQUERY ) THEN
            WORK( 1 ) = N*NB
         END IF
      END IF
      IF( TQUERY .OR. WQUERY ) THEN
         RETURN
      END IF
*
*     Quick return
*
      IF ( N.EQ.0 ) THEN
         RETURN
      ENDIF
*
*     Determine the number of the block size
*
      LDTB = LTB/N
      IF( LDTB .LT. 3*NB+1 ) THEN
         NB = (LDTB-1)/3
      END IF
      IF( LWORK .LT. NB*N ) THEN
         NB = LWORK/N
      END IF
*
*     Determine the number of the block columns
*
      NT = (N+NB-1)/NB
      TD = 2*NB
      KB = MIN(NB, N)
*
*     Initialize vectors/matrices
*
      DO J = 1, KB
         IPIV( J ) = J
      END DO
*
*     Save NB
*
      TB( 1 ) = NB
*
      IF( UPPER ) THEN
*
*        .....................................................
*        Factorize A as U**T*D*U using the upper triangle of A
*        .....................................................
*
         DO J = 0, NT-1
*
*           Generate Jth column of W and H
*
            KB = MIN(NB, N-J*NB)
            DO I = 1, J-1
               IF( I.EQ.1 ) THEN
*                  H(I,J) = T(I,I)*U(I,J) + T(I+1,I)*U(I+1,J)
                  IF( I .EQ. (J-1) ) THEN
                     JB = NB+KB
                  ELSE
                     JB = 2*NB
                  END IF
                  CALL ZGEMM( 'NoTranspose', 'NoTranspose', NB, KB, JB, CONE,  TB( TD+1 + (I*NB)*LDTB ), LDTB-1, A( (I-1)*NB+1, J*NB+1 ), LDA, CZERO, WORK( I*NB+1 ), N )
               ELSE
*                 H(I,J) = T(I,I-1)*U(I-1,J) + T(I,I)*U(I,J) + T(I,I+1)*U(I+1,J)
                  IF( I .EQ. (J-1) ) THEN
                     JB = 2*NB+KB
                  ELSE
                     JB = 3*NB
                  END IF
                  CALL ZGEMM( 'NoTranspose', 'NoTranspose', NB, KB, JB, CONE,  TB( TD+NB+1 + ((I-1)*NB)*LDTB ), LDTB-1, A( (I-2)*NB+1, J*NB+1 ), LDA, CZERO, WORK( I*NB+1 ), N )
               END IF
            END DO
*
*           Compute T(J,J)
*
            CALL ZLACPY( 'Upper', KB, KB, A( J*NB+1, J*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
            IF( J.GT.1 ) THEN
*              T(J,J) = U(1:J,J)'*H(1:J)
               CALL ZGEMM( 'Transpose', 'NoTranspose', KB, KB, (J-1)*NB, -CONE, A( 1, J*NB+1 ), LDA, WORK( NB+1 ), N, CONE, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
*              T(J,J) += U(J,J)'*T(J,J-1)*U(J-1,J)
               CALL ZGEMM( 'Transpose', 'NoTranspose', KB, NB, KB, CONE,  A( (J-1)*NB+1, J*NB+1 ), LDA, TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, CZERO, WORK( 1 ), N )                CALL ZGEMM( 'NoTranspose', 'NoTranspose', KB, KB, NB, -CONE, WORK( 1 ), N, A( (J-2)*NB+1, J*NB+1 ), LDA, CONE, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
            END IF
*
*           Expand T(J,J) into full format
*
            DO I = 1, KB
               DO K = I+1, KB
                  TB( TD+(K-I)+1 + (J*NB+I-1)*LDTB ) = TB( TD-(K-(I+1)) + (J*NB+K-1)*LDTB )
               END DO
            END DO
            IF( J.GT.0 ) THEN
c               CALL CHEGST( 1, 'Upper', KB,
c     $                      TB( TD+1 + (J*NB)*LDTB ), LDTB-1,
c     $                      A( (J-1)*NB+1, J*NB+1 ), LDA, IINFO )
               CALL ZTRSM( 'L', 'U', 'T', 'N', KB, KB, CONE, A( (J-1)*NB+1, J*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )                CALL ZTRSM( 'R', 'U', 'N', 'N', KB, KB, CONE, A( (J-1)*NB+1, J*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
            END IF
*
            IF( J.LT.NT-1 ) THEN
               IF( J.GT.0 ) THEN
*
*                 Compute H(J,J)
*
                  IF( J.EQ.1 ) THEN
                     CALL ZGEMM( 'NoTranspose', 'NoTranspose', KB, KB, KB, CONE,  TB( TD+1 + (J*NB)*LDTB ), LDTB-1, A( (J-1)*NB+1, J*NB+1 ), LDA, CZERO, WORK( J*NB+1 ), N )
                  ELSE
                     CALL ZGEMM( 'NoTranspose', 'NoTranspose', KB, KB, NB+KB, CONE, TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, A( (J-2)*NB+1, J*NB+1 ), LDA, CZERO, WORK( J*NB+1 ), N )
                  END IF
*
*                 Update with the previous column
*
                  CALL ZGEMM( 'Transpose', 'NoTranspose', NB, N-(J+1)*NB, J*NB, -CONE, WORK( NB+1 ), N, A( 1, (J+1)*NB+1 ), LDA, CONE, A( J*NB+1, (J+1)*NB+1 ), LDA )
               END IF
*
*              Copy panel to workspace to call ZGETRF
*
               DO K = 1, NB
                   CALL ZCOPY( N-(J+1)*NB, A( J*NB+K, (J+1)*NB+1 ), LDA, WORK( 1+(K-1)*N ), 1 )
               END DO
*
*              Factorize panel
*
               CALL ZGETRF( N-(J+1)*NB, NB,  WORK, N, IPIV( (J+1)*NB+1 ), IINFO )
c               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN
c                  INFO = IINFO+(J+1)*NB
c               END IF
*
*              Copy panel back
*
               DO K = 1, NB
                   CALL ZCOPY( N-(J+1)*NB, WORK( 1+(K-1)*N ), 1, A( J*NB+K, (J+1)*NB+1 ), LDA )
               END DO
*
*              Compute T(J+1, J), zero out for GEMM update
*
               KB = MIN(NB, N-(J+1)*NB)
               CALL ZLASET( 'Full', KB, NB, CZERO, CZERO,  TB( TD+NB+1 + (J*NB)*LDTB), LDTB-1 )                CALL ZLACPY( 'Upper', KB, NB, WORK, N, TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 )
               IF( J.GT.0 ) THEN
                  CALL ZTRSM( 'R', 'U', 'N', 'U', KB, NB, CONE, A( (J-1)*NB+1, J*NB+1 ), LDA, TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 )
               END IF
*
*              Copy T(J,J+1) into T(J+1, J), both upper/lower for GEMM
*              updates
*
               DO K = 1, NB
                  DO I = 1, KB
                     TB( TD-NB+K-I+1 + (J*NB+NB+I-1)*LDTB ) = TB( TD+NB+I-K+1 + (J*NB+K-1)*LDTB )
                  END DO
               END DO
               CALL ZLASET( 'Lower', KB, NB, CZERO, CONE,  A( J*NB+1, (J+1)*NB+1), LDA )
*
*              Apply pivots to trailing submatrix of A
*
               DO K = 1, KB
*                 > Adjust ipiv
                  IPIV( (J+1)*NB+K ) = IPIV( (J+1)*NB+K ) + (J+1)*NB
*
                  I1 = (J+1)*NB+K
                  I2 = IPIV( (J+1)*NB+K )
                  IF( I1.NE.I2 ) THEN
*                    > Apply pivots to previous columns of L
                     CALL ZSWAP( K-1, A( (J+1)*NB+1, I1 ), 1,  A( (J+1)*NB+1, I2 ), 1 )
*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                     IF( I2.GT.(I1+1) ) CALL ZSWAP( I2-I1-1, A( I1, I1+1 ), LDA, A( I1+1, I2 ), 1 )
*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                     IF( I2.LT.N ) CALL ZSWAP( N-I2, A( I1, I2+1 ), LDA, A( I2, I2+1 ), LDA )
*                    > Swap A(I1, I1) with A(I2, I2)
                     PIV = A( I1, I1 )
                     A( I1, I1 ) = A( I2, I2 )
                     A( I2, I2 ) = PIV
*                    > Apply pivots to previous columns of L
                     IF( J.GT.0 ) THEN
                        CALL ZSWAP( J*NB, A( 1, I1 ), 1, A( 1, I2 ), 1 )
                     END IF
                  ENDIF
               END DO
            END IF
         END DO
      ELSE
*
*        .....................................................
*        Factorize A as L*D*L**T using the lower triangle of A
*        .....................................................
*
         DO J = 0, NT-1
*
*           Generate Jth column of W and H
*
            KB = MIN(NB, N-J*NB)
            DO I = 1, J-1
               IF( I.EQ.1 ) THEN
*                  H(I,J) = T(I,I)*L(J,I)' + T(I+1,I)'*L(J,I+1)'
                  IF( I .EQ. (J-1) ) THEN
                     JB = NB+KB
                  ELSE
                     JB = 2*NB
                  END IF
                  CALL ZGEMM( 'NoTranspose', 'Transpose', NB, KB, JB, CONE, TB( TD+1 + (I*NB)*LDTB ), LDTB-1, A( J*NB+1, (I-1)*NB+1 ), LDA, CZERO, WORK( I*NB+1 ), N )
               ELSE
*                 H(I,J) = T(I,I-1)*L(J,I-1)' + T(I,I)*L(J,I)' + T(I,I+1)*L(J,I+1)'
                  IF( I .EQ. (J-1) ) THEN
                     JB = 2*NB+KB
                  ELSE
                     JB = 3*NB
                  END IF
                  CALL ZGEMM( 'NoTranspose', 'Transpose', NB, KB, JB, CONE,  TB( TD+NB+1 + ((I-1)*NB)*LDTB ), LDTB-1, A( J*NB+1, (I-2)*NB+1 ), LDA, CZERO, WORK( I*NB+1 ), N )
               END IF
            END DO
*
*           Compute T(J,J)
*
            CALL ZLACPY( 'Lower', KB, KB, A( J*NB+1, J*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
            IF( J.GT.1 ) THEN
*              T(J,J) = L(J,1:J)*H(1:J)
               CALL ZGEMM( 'NoTranspose', 'NoTranspose', KB, KB, (J-1)*NB, -CONE, A( J*NB+1, 1 ), LDA, WORK( NB+1 ), N, CONE, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
*              T(J,J) += L(J,J)*T(J,J-1)*L(J,J-1)'
               CALL ZGEMM( 'NoTranspose', 'NoTranspose', KB, NB, KB, CONE,  A( J*NB+1, (J-1)*NB+1 ), LDA, TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, CZERO, WORK( 1 ), N )                CALL ZGEMM( 'NoTranspose', 'Transpose', KB, KB, NB, -CONE, WORK( 1 ), N, A( J*NB+1, (J-2)*NB+1 ), LDA, CONE, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
            END IF
*
*           Expand T(J,J) into full format
*
            DO I = 1, KB
               DO K = I+1, KB
                  TB( TD-(K-(I+1)) + (J*NB+K-1)*LDTB ) = TB( TD+(K-I)+1 + (J*NB+I-1)*LDTB )
               END DO
            END DO
            IF( J.GT.0 ) THEN
c               CALL CHEGST( 1, 'Lower', KB,
c     $                      TB( TD+1 + (J*NB)*LDTB ), LDTB-1,
c     $                      A( J*NB+1, (J-1)*NB+1 ), LDA, IINFO )
               CALL ZTRSM( 'L', 'L', 'N', 'N', KB, KB, CONE, A( J*NB+1, (J-1)*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )                CALL ZTRSM( 'R', 'L', 'T', 'N', KB, KB, CONE, A( J*NB+1, (J-1)*NB+1 ), LDA, TB( TD+1 + (J*NB)*LDTB ), LDTB-1 )
            END IF
*
*           Symmetrize T(J,J)
*
            DO I = 1, KB
               DO K = I+1, KB
                  TB( TD-(K-(I+1)) + (J*NB+K-1)*LDTB ) = TB( TD+(K-I)+1 + (J*NB+I-1)*LDTB )
               END DO
            END DO
*
            IF( J.LT.NT-1 ) THEN
               IF( J.GT.0 ) THEN
*
*                 Compute H(J,J)
*
                  IF( J.EQ.1 ) THEN
                     CALL ZGEMM( 'NoTranspose', 'Transpose', KB, KB, KB, CONE,  TB( TD+1 + (J*NB)*LDTB ), LDTB-1, A( J*NB+1, (J-1)*NB+1 ), LDA, CZERO, WORK( J*NB+1 ), N )
                  ELSE
                     CALL ZGEMM( 'NoTranspose', 'Transpose', KB, KB, NB+KB, CONE, TB( TD+NB+1 + ((J-1)*NB)*LDTB ), LDTB-1, A( J*NB+1, (J-2)*NB+1 ), LDA, CZERO, WORK( J*NB+1 ), N )
                  END IF
*
*                 Update with the previous column
*
                  CALL ZGEMM( 'NoTranspose', 'NoTranspose', N-(J+1)*NB, NB, J*NB, -CONE, A( (J+1)*NB+1, 1 ), LDA, WORK( NB+1 ), N, CONE, A( (J+1)*NB+1, J*NB+1 ), LDA )
               END IF
*
*              Factorize panel
*
               CALL ZGETRF( N-(J+1)*NB, NB,  A( (J+1)*NB+1, J*NB+1 ), LDA, IPIV( (J+1)*NB+1 ), IINFO )
c               IF (IINFO.NE.0 .AND. INFO.EQ.0) THEN
c                  INFO = IINFO+(J+1)*NB
c               END IF
*
*              Compute T(J+1, J), zero out for GEMM update
*
               KB = MIN(NB, N-(J+1)*NB)
               CALL ZLASET( 'Full', KB, NB, CZERO, CZERO,  TB( TD+NB+1 + (J*NB)*LDTB), LDTB-1 )                CALL ZLACPY( 'Upper', KB, NB, A( (J+1)*NB+1, J*NB+1 ), LDA, TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 )
               IF( J.GT.0 ) THEN
                  CALL ZTRSM( 'R', 'L', 'T', 'U', KB, NB, CONE, A( J*NB+1, (J-1)*NB+1 ), LDA, TB( TD+NB+1 + (J*NB)*LDTB ), LDTB-1 )
               END IF
*
*              Copy T(J+1,J) into T(J, J+1), both upper/lower for GEMM
*              updates
*
               DO K = 1, NB
                  DO I = 1, KB
                     TB( TD-NB+K-I+1 + (J*NB+NB+I-1)*LDTB ) = TB( TD+NB+I-K+1 + (J*NB+K-1)*LDTB )
                  END DO
               END DO
               CALL ZLASET( 'Upper', KB, NB, CZERO, CONE,  A( (J+1)*NB+1, J*NB+1 ), LDA )
*
*              Apply pivots to trailing submatrix of A
*
               DO K = 1, KB
*                 > Adjust ipiv
                  IPIV( (J+1)*NB+K ) = IPIV( (J+1)*NB+K ) + (J+1)*NB
*
                  I1 = (J+1)*NB+K
                  I2 = IPIV( (J+1)*NB+K )
                  IF( I1.NE.I2 ) THEN
*                    > Apply pivots to previous columns of L
                     CALL ZSWAP( K-1, A( I1, (J+1)*NB+1 ), LDA,  A( I2, (J+1)*NB+1 ), LDA )
*                    > Swap A(I1+1:M, I1) with A(I2, I1+1:M)
                     IF( I2.GT.(I1+1) ) CALL ZSWAP( I2-I1-1, A( I1+1, I1 ), 1, A( I2, I1+1 ), LDA )
*                    > Swap A(I2+1:M, I1) with A(I2+1:M, I2)
                     IF( I2.LT.N ) CALL ZSWAP( N-I2, A( I2+1, I1 ), 1, A( I2+1, I2 ), 1 )
*                    > Swap A(I1, I1) with A(I2, I2)
                     PIV = A( I1, I1 )
                     A( I1, I1 ) = A( I2, I2 )
                     A( I2, I2 ) = PIV
*                    > Apply pivots to previous columns of L
                     IF( J.GT.0 ) THEN
                        CALL ZSWAP( J*NB, A( I1, 1 ), LDA, A( I2, 1 ), LDA )
                     END IF
                  ENDIF
               END DO
*
*              Apply pivots to previous columns of L
*
c               CALL ZLASWP( J*NB, A( 1, 1 ), LDA,
c     $                     (J+1)*NB+1, (J+1)*NB+KB, IPIV, 1 )
            END IF
         END DO
      END IF
*
*     Factor the band matrix
      CALL ZGBTRF( N, N, NB, NB, TB, LDTB, IPIV2, INFO )
*
      RETURN
*
*     End of ZSYTRF_AA_2STAGE
*
      END
