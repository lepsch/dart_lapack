      SUBROUTINE SSYTRF_RK( UPLO, N, A, LDA, E, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      int                IPIV( * )
      REAL               A( LDA, * ), E( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      bool               LQUERY, UPPER;
      int                I, IINFO, IP, IWS, K, KB, LDWORK, LWKOPT, NB, NBMIN
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV
      REAL               SROUNDUP_LWORK
      EXTERNAL           LSAME, ILAENV, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASYF_RK, SSYTF2_RK, SSWAP,  XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size
*
         NB = ILAENV( 1, 'SSYTRF_RK', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, N*NB )
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYTRF_RK', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = LDWORK*NB
         IF( LWORK.LT.IWS ) THEN
            NB = MAX( LWORK / LDWORK, 1 )
            NBMIN = MAX( 2, ILAENV( 2, 'SSYTRF_RK', UPLO, N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = 1
      END IF
      IF( NB.LT.NBMIN ) NB = N
*
      IF( UPPER ) THEN
*
*        Factorize A as U*D*U**T using the upper triangle of A
*
*        K is the main loop index, decreasing from N to 1 in steps of
*        KB, where KB is the number of columns factorized by SLASYF_RK;
*        KB is either NB or NB-1, or K for the last block
*
         K = N
   10    CONTINUE
*
*        If K < 1, exit from loop
*
         IF( K.LT.1 ) GO TO 15
*
         IF( K.GT.NB ) THEN
*
*           Factorize columns k-kb+1:k of A and use blocked code to
*           update columns 1:k-kb
*
            CALL SLASYF_RK( UPLO, K, NB, KB, A, LDA, E, IPIV, WORK, LDWORK, IINFO )
         ELSE
*
*           Use unblocked code to factorize columns 1:k of A
*
            CALL SSYTF2_RK( UPLO, K, A, LDA, E, IPIV, IINFO )
            KB = K
         END IF
*
*        Set INFO on the first occurrence of a zero pivot
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO
*
*        No need to adjust IPIV
*
*
*        Apply permutations to the leading panel 1:k-1
*
*        Read IPIV from the last block factored, i.e.
*        indices  k-kb+1:k and apply row permutations to the
*        last k+1 colunms k+1:N after that block
*        (We can do the simple loop over IPIV with decrement -1,
*        since the ABS value of IPIV( I ) represents the row index
*        of the interchange with row i in both 1x1 and 2x2 pivot cases)
*
         IF( K.LT.N ) THEN
            DO I = K, ( K - KB + 1 ), -1
               IP = ABS( IPIV( I ) )
               IF( IP.NE.I ) THEN
                  CALL SSWAP( N-K, A( I, K+1 ), LDA, A( IP, K+1 ), LDA )
               END IF
            END DO
         END IF
*
*        Decrease K and return to the start of the main loop
*
         K = K - KB
         GO TO 10
*
*        This label is the exit from main loop over K decreasing
*        from N to 1 in steps of KB
*
   15    CONTINUE
*
      ELSE
*
*        Factorize A as L*D*L**T using the lower triangle of A
*
*        K is the main loop index, increasing from 1 to N in steps of
*        KB, where KB is the number of columns factorized by SLASYF_RK;
*        KB is either NB or NB-1, or N-K+1 for the last block
*
         K = 1
   20    CONTINUE
*
*        If K > N, exit from loop
*
         IF( K.GT.N ) GO TO 35
*
         IF( K.LE.N-NB ) THEN
*
*           Factorize columns k:k+kb-1 of A and use blocked code to
*           update columns k+kb:n
*
            CALL SLASYF_RK( UPLO, N-K+1, NB, KB, A( K, K ), LDA, E( K ), IPIV( K ), WORK, LDWORK, IINFO )


         ELSE
*
*           Use unblocked code to factorize columns k:n of A
*
            CALL SSYTF2_RK( UPLO, N-K+1, A( K, K ), LDA, E( K ), IPIV( K ), IINFO )
            KB = N - K + 1
*
         END IF
*
*        Set INFO on the first occurrence of a zero pivot
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + K - 1
*
*        Adjust IPIV
*
         DO I = K, K + KB - 1
            IF( IPIV( I ).GT.0 ) THEN
               IPIV( I ) = IPIV( I ) + K - 1
            ELSE
               IPIV( I ) = IPIV( I ) - K + 1
            END IF
         END DO
*
*        Apply permutations to the leading panel 1:k-1
*
*        Read IPIV from the last block factored, i.e.
*        indices  k:k+kb-1 and apply row permutations to the
*        first k-1 colunms 1:k-1 before that block
*        (We can do the simple loop over IPIV with increment 1,
*        since the ABS value of IPIV( I ) represents the row index
*        of the interchange with row i in both 1x1 and 2x2 pivot cases)
*
         IF( K.GT.1 ) THEN
            DO I = K, ( K + KB - 1 ), 1
               IP = ABS( IPIV( I ) )
               IF( IP.NE.I ) THEN
                  CALL SSWAP( K-1, A( I, 1 ), LDA, A( IP, 1 ), LDA )
               END IF
            END DO
         END IF
*
*        Increase K and return to the start of the main loop
*
         K = K + KB
         GO TO 20
*
*        This label is the exit from main loop over K increasing
*        from 1 to N in steps of KB
*
   35    CONTINUE
*
*     End Lower
*
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      RETURN
*
*     End of SSYTRF_RK
*
      END
