      SUBROUTINE ZLASYF_AA( UPLO, J1, M, NB, A, LDA, IPIV, H, LDH, WORK )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                M, NB, J1, LDA, LDH
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      COMPLEX*16         A( LDA, * ), H( LDH, * ), WORK( * )
*     ..
*
*  =====================================================================
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*
*     .. Local Scalars ..
      int                J, K, K1, I1, I2, MJ
      COMPLEX*16         PIV, ALPHA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      int                IZAMAX, ILAENV
      EXTERNAL           LSAME, ILAENV, IZAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZAXPY, ZSCAL, ZCOPY, ZSWAP, ZLASET, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      J = 1
*
*     K1 is the first column of the panel to be factorized
*     i.e.,  K1 is 2 for the first block column, and 1 for the rest of the blocks
*
      K1 = (2-J1)+1
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        .....................................................
*        Factorize A as U**T*D*U using the upper triangle of A
*        .....................................................
*
 10      CONTINUE
         IF ( J.GT.MIN(M, NB) ) GO TO 20
*
*        K is the column to be factorized
*         when being called from ZSYTRF_AA,
*         > for the first block column, J1 is 1, hence J1+J-1 is J,
*         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1,
*
         K = J1+J-1
         IF( J.EQ.M ) THEN
*
*            Only need to compute T(J, J)
*
             MJ = 1
         ELSE
             MJ = M-J+1
         END IF
*
*        H(J:M, J) := A(J, J:M) - H(J:M, 1:(J-1)) * L(J1:(J-1), J),
*         where H(J:M, J) has been initialized to be A(J, J:M)
*
         IF( K.GT.2 ) THEN
*
*        K is the column to be factorized
*         > for the first block column, K is J, skipping the first two
*           columns
*         > for the rest of the columns, K is J+1, skipping only the
*           first column
*
            CALL ZGEMV( 'No transpose', MJ, J-K1, -ONE, H( J, K1 ), LDH, A( 1, J ), 1, ONE, H( J, J ), 1 )
         END IF
*
*        Copy H(i:M, i) into WORK
*
         CALL ZCOPY( MJ, H( J, J ), 1, WORK( 1 ), 1 )
*
         IF( J.GT.K1 ) THEN
*
*           Compute WORK := WORK - L(J-1, J:M) * T(J-1,J),
*            where A(J-1, J) stores T(J-1, J) and A(J-2, J:M) stores U(J-1, J:M)
*
            ALPHA = -A( K-1, J )
            CALL ZAXPY( MJ, ALPHA, A( K-2, J ), LDA, WORK( 1 ), 1 )
         END IF
*
*        Set A(J, J) = T(J, J)
*
         A( K, J ) = WORK( 1 )
*
         IF( J.LT.M ) THEN
*
*           Compute WORK(2:M) = T(J, J) L(J, (J+1):M)
*            where A(J, J) stores T(J, J) and A(J-1, (J+1):M) stores U(J, (J+1):M)
*
            IF( K.GT.1 ) THEN
               ALPHA = -A( K, J )
               CALL ZAXPY( M-J, ALPHA, A( K-1, J+1 ), LDA, WORK( 2 ), 1 )
            ENDIF
*
*           Find max(|WORK(2:M)|)
*
            I2 = IZAMAX( M-J, WORK( 2 ), 1 ) + 1
            PIV = WORK( I2 )
*
*           Apply symmetric pivot
*
            IF( (I2.NE.2) .AND. (PIV.NE.0) ) THEN
*
*              Swap WORK(I1) and WORK(I2)
*
               I1 = 2
               WORK( I2 ) = WORK( I1 )
               WORK( I1 ) = PIV
*
*              Swap A(I1, I1+1:M) with A(I1+1:M, I2)
*
               I1 = I1+J-1
               I2 = I2+J-1
               CALL ZSWAP( I2-I1-1, A( J1+I1-1, I1+1 ), LDA, A( J1+I1, I2 ), 1 )
*
*              Swap A(I1, I2+1:M) with A(I2, I2+1:M)
*
               IF( I2.LT.M ) CALL ZSWAP( M-I2, A( J1+I1-1, I2+1 ), LDA, A( J1+I2-1, I2+1 ), LDA )
*
*              Swap A(I1, I1) with A(I2,I2)
*
               PIV = A( I1+J1-1, I1 )
               A( J1+I1-1, I1 ) = A( J1+I2-1, I2 )
               A( J1+I2-1, I2 ) = PIV
*
*              Swap H(I1, 1:J1) with H(I2, 1:J1)
*
               CALL ZSWAP( I1-1, H( I1, 1 ), LDH, H( I2, 1 ), LDH )
               IPIV( I1 ) = I2
*
               IF( I1.GT.(K1-1) ) THEN
*
*                 Swap L(1:I1-1, I1) with L(1:I1-1, I2),
*                  skipping the first column
*
                  CALL ZSWAP( I1-K1+1, A( 1, I1 ), 1, A( 1, I2 ), 1 )
               END IF
            ELSE
               IPIV( J+1 ) = J+1
            ENDIF
*
*           Set A(J, J+1) = T(J, J+1)
*
            A( K, J+1 ) = WORK( 2 )
*
            IF( J.LT.NB ) THEN
*
*              Copy A(J+1:M, J+1) into H(J:M, J),
*
               CALL ZCOPY( M-J, A( K+1, J+1 ), LDA, H( J+1, J+1 ), 1 )
            END IF
*
*           Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1),
*            where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1)
*
            IF( J.LT.(M-1) ) THEN
               IF( A( K, J+1 ).NE.ZERO ) THEN
                  ALPHA = ONE / A( K, J+1 )
                  CALL ZCOPY( M-J-1, WORK( 3 ), 1, A( K, J+2 ), LDA )
                  CALL ZSCAL( M-J-1, ALPHA, A( K, J+2 ), LDA )
               ELSE
                  CALL ZLASET( 'Full', 1, M-J-1, ZERO, ZERO, A( K, J+2 ), LDA)
               END IF
            END IF
         END IF
         J = J + 1
         GO TO 10
 20      CONTINUE
*
      ELSE
*
*        .....................................................
*        Factorize A as L*D*L**T using the lower triangle of A
*        .....................................................
*
 30      CONTINUE
         IF( J.GT.MIN( M, NB ) ) GO TO 40
*
*        K is the column to be factorized
*         when being called from ZSYTRF_AA,
*         > for the first block column, J1 is 1, hence J1+J-1 is J,
*         > for the rest of the columns, J1 is 2, and J1+J-1 is J+1,
*
         K = J1+J-1
         IF( J.EQ.M ) THEN
*
*            Only need to compute T(J, J)
*
             MJ = 1
         ELSE
             MJ = M-J+1
         END IF
*
*        H(J:M, J) := A(J:M, J) - H(J:M, 1:(J-1)) * L(J, J1:(J-1))^T,
*         where H(J:M, J) has been initialized to be A(J:M, J)
*
         IF( K.GT.2 ) THEN
*
*        K is the column to be factorized
*         > for the first block column, K is J, skipping the first two
*           columns
*         > for the rest of the columns, K is J+1, skipping only the
*           first column
*
            CALL ZGEMV( 'No transpose', MJ, J-K1, -ONE, H( J, K1 ), LDH, A( J, 1 ), LDA, ONE, H( J, J ), 1 )
         END IF
*
*        Copy H(J:M, J) into WORK
*
         CALL ZCOPY( MJ, H( J, J ), 1, WORK( 1 ), 1 )
*
         IF( J.GT.K1 ) THEN
*
*           Compute WORK := WORK - L(J:M, J-1) * T(J-1,J),
*            where A(J-1, J) = T(J-1, J) and A(J, J-2) = L(J, J-1)
*
            ALPHA = -A( J, K-1 )
            CALL ZAXPY( MJ, ALPHA, A( J, K-2 ), 1, WORK( 1 ), 1 )
         END IF
*
*        Set A(J, J) = T(J, J)
*
         A( J, K ) = WORK( 1 )
*
         IF( J.LT.M ) THEN
*
*           Compute WORK(2:M) = T(J, J) L((J+1):M, J)
*            where A(J, J) = T(J, J) and A((J+1):M, J-1) = L((J+1):M, J)
*
            IF( K.GT.1 ) THEN
               ALPHA = -A( J, K )
               CALL ZAXPY( M-J, ALPHA, A( J+1, K-1 ), 1, WORK( 2 ), 1 )
            ENDIF
*
*           Find max(|WORK(2:M)|)
*
            I2 = IZAMAX( M-J, WORK( 2 ), 1 ) + 1
            PIV = WORK( I2 )
*
*           Apply symmetric pivot
*
            IF( (I2.NE.2) .AND. (PIV.NE.0) ) THEN
*
*              Swap WORK(I1) and WORK(I2)
*
               I1 = 2
               WORK( I2 ) = WORK( I1 )
               WORK( I1 ) = PIV
*
*              Swap A(I1+1:M, I1) with A(I2, I1+1:M)
*
               I1 = I1+J-1
               I2 = I2+J-1
               CALL ZSWAP( I2-I1-1, A( I1+1, J1+I1-1 ), 1, A( I2, J1+I1 ), LDA )
*
*              Swap A(I2+1:M, I1) with A(I2+1:M, I2)
*
               IF( I2.LT.M ) CALL ZSWAP( M-I2, A( I2+1, J1+I1-1 ), 1, A( I2+1, J1+I2-1 ), 1 )
*
*              Swap A(I1, I1) with A(I2, I2)
*
               PIV = A( I1, J1+I1-1 )
               A( I1, J1+I1-1 ) = A( I2, J1+I2-1 )
               A( I2, J1+I2-1 ) = PIV
*
*              Swap H(I1, I1:J1) with H(I2, I2:J1)
*
               CALL ZSWAP( I1-1, H( I1, 1 ), LDH, H( I2, 1 ), LDH )
               IPIV( I1 ) = I2
*
               IF( I1.GT.(K1-1) ) THEN
*
*                 Swap L(1:I1-1, I1) with L(1:I1-1, I2),
*                  skipping the first column
*
                  CALL ZSWAP( I1-K1+1, A( I1, 1 ), LDA, A( I2, 1 ), LDA )
               END IF
            ELSE
               IPIV( J+1 ) = J+1
            ENDIF
*
*           Set A(J+1, J) = T(J+1, J)
*
            A( J+1, K ) = WORK( 2 )
*
            IF( J.LT.NB ) THEN
*
*              Copy A(J+1:M, J+1) into H(J+1:M, J),
*
               CALL ZCOPY( M-J, A( J+1, K+1 ), 1, H( J+1, J+1 ), 1 )
            END IF
*
*           Compute L(J+2, J+1) = WORK( 3:M ) / T(J, J+1),
*            where A(J, J+1) = T(J, J+1) and A(J+2:M, J) = L(J+2:M, J+1)
*
            IF( J.LT.(M-1) ) THEN
               IF( A( J+1, K ).NE.ZERO ) THEN
                  ALPHA = ONE / A( J+1, K )
                  CALL ZCOPY( M-J-1, WORK( 3 ), 1, A( J+2, K ), 1 )
                  CALL ZSCAL( M-J-1, ALPHA, A( J+2, K ), 1 )
               ELSE
                  CALL ZLASET( 'Full', M-J-1, 1, ZERO, ZERO, A( J+2, K ), LDA )
               END IF
            END IF
         END IF
         J = J + 1
         GO TO 30
 40      CONTINUE
      END IF
      RETURN
*
*     End of ZLASYF_AA
*
      END
