      SUBROUTINE CSYTRF_ROOK( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY, UPPER;
      int                IINFO, IWS, J, K, KB, LDWORK, LWKOPT, NB, NBMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASYF_ROOK, CSYTF2_ROOK, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      } else if ( LWORK.LT.1 .AND. .NOT.LQUERY ) {
         INFO = -7
      }

      if ( INFO.EQ.0 ) {

         // Determine the block size

         NB = ILAENV( 1, 'CSYTRF_ROOK', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, N*NB )
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CSYTRF_ROOK', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      NBMIN = 2
      LDWORK = N
      if ( NB.GT.1 .AND. NB.LT.N ) {
         IWS = LDWORK*NB
         if ( LWORK.LT.IWS ) {
            NB = MAX( LWORK / LDWORK, 1 )
            NBMIN = MAX( 2, ILAENV( 2, 'CSYTRF_ROOK', UPLO, N, -1, -1, -1 ) )
         }
      } else {
         IWS = 1
      }
      IF( NB.LT.NBMIN ) NB = N

      if ( UPPER ) {

         // Factorize A as U*D*U**T using the upper triangle of A

         // K is the main loop index, decreasing from N to 1 in steps of
         // KB, where KB is the number of columns factorized by CLASYF_ROOK;
         // KB is either NB or NB-1, or K for the last block

         K = N
   10    CONTINUE

         // If K < 1, exit from loop

         IF( K.LT.1 ) GO TO 40

         if ( K.GT.NB ) {

            // Factorize columns k-kb+1:k of A and use blocked code to
            // update columns 1:k-kb

            CALL CLASYF_ROOK( UPLO, K, NB, KB, A, LDA, IPIV, WORK, LDWORK, IINFO )
         } else {

            // Use unblocked code to factorize columns 1:k of A

            CALL CSYTF2_ROOK( UPLO, K, A, LDA, IPIV, IINFO )
            KB = K
         }

         // Set INFO on the first occurrence of a zero pivot

         IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO

         // No need to adjust IPIV

         // Decrease K and return to the start of the main loop

         K = K - KB
         GO TO 10

      } else {

         // Factorize A as L*D*L**T using the lower triangle of A

         // K is the main loop index, increasing from 1 to N in steps of
         // KB, where KB is the number of columns factorized by CLASYF_ROOK;
         // KB is either NB or NB-1, or N-K+1 for the last block

         K = 1
   20    CONTINUE

         // If K > N, exit from loop

         IF( K.GT.N ) GO TO 40

         if ( K.LE.N-NB ) {

            // Factorize columns k:k+kb-1 of A and use blocked code to
            // update columns k+kb:n

            CALL CLASYF_ROOK( UPLO, N-K+1, NB, KB, A( K, K ), LDA, IPIV( K ), WORK, LDWORK, IINFO )
         } else {

            // Use unblocked code to factorize columns k:n of A

            CALL CSYTF2_ROOK( UPLO, N-K+1, A( K, K ), LDA, IPIV( K ), IINFO )
            KB = N - K + 1
         }

         // Set INFO on the first occurrence of a zero pivot

         IF( INFO.EQ.0 .AND. IINFO.GT.0 ) INFO = IINFO + K - 1

         // Adjust IPIV

         DO 30 J = K, K + KB - 1
            if ( IPIV( J ).GT.0 ) {
               IPIV( J ) = IPIV( J ) + K - 1
            } else {
               IPIV( J ) = IPIV( J ) - K + 1
            }
   30    CONTINUE

         // Increase K and return to the start of the main loop

         K = K + KB
         GO TO 20

      }

   40 CONTINUE
      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
      RETURN

      // End of CSYTRF_ROOK

      }
