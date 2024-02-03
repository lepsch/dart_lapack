      SUBROUTINE DSYGV_2STAGE( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDA, LDB, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), W( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER, WANTZ;
      String             TRANS;
      int                NEIG, LWMIN, LHTRD, LWTRD, KD, IB;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      // EXTERNAL LSAME, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DPOTRF, DSYGST, DTRMM, DTRSM, XERBLA, DSYEV_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 )

      INFO = 0
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1
      } else if ( .NOT.( LSAME( JOBZ, 'N' ) ) ) {
         INFO = -2
      } else if ( .NOT.( UPPER || LSAME( UPLO, 'L' ) ) ) {
         INFO = -3
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -8
      }

      if ( INFO == 0 ) {
         KD    = ILAENV2STAGE( 1, 'DSYTRD_2STAGE', JOBZ, N, -1, -1, -1 )
         IB    = ILAENV2STAGE( 2, 'DSYTRD_2STAGE', JOBZ, N, KD, -1, -1 )
         LHTRD = ILAENV2STAGE( 3, 'DSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )
         LWTRD = ILAENV2STAGE( 4, 'DSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )
         LWMIN = 2*N + LHTRD + LWTRD
         WORK( 1 )  = LWMIN

         if ( LWORK < LWMIN && .NOT.LQUERY ) {
            INFO = -11
         }
      }

      if ( INFO != 0 ) {
         xerbla('DSYGV_2STAGE ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Form a Cholesky factorization of B.

      dpotrf(UPLO, N, B, LDB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem and solve.

      dsygst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO );
      dsyev_2stage(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO );

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         NEIG = N
         if (INFO > 0) NEIG = INFO - 1;
         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N'
            } else {
               TRANS = 'T'
            }

            dtrsm('Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, B, LDB, A, LDA );

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**T*y

            if ( UPPER ) {
               TRANS = 'T'
            } else {
               TRANS = 'N'
            }

            dtrmm('Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, B, LDB, A, LDA );
         }
      }

      WORK( 1 ) = LWMIN
      RETURN

      // End of DSYGV_2STAGE

      }
