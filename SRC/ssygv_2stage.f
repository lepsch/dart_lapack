      SUBROUTINE SSYGV_2STAGE( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDA, LDB, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER, WANTZ;
      String             TRANS;
      int                NEIG, LWMIN, LHTRD, LWTRD, KD, IB;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SPOTRF, SSYGST, STRMM, STRSM, XERBLA, SSYEV_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )

      INFO = 0
      if ( ITYPE.LT.1 .OR. ITYPE.GT.3 ) {
         INFO = -1
      } else if ( .NOT.( LSAME( JOBZ, 'N' ) ) ) {
         INFO = -2
      } else if ( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -8
      }

      if ( INFO.EQ.0 ) {
         KD    = ILAENV2STAGE( 1, 'SSYTRD_2STAGE', JOBZ, N, -1, -1, -1 )
         IB    = ILAENV2STAGE( 2, 'SSYTRD_2STAGE', JOBZ, N, KD, -1, -1 )
         LHTRD = ILAENV2STAGE( 3, 'SSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )
         LWTRD = ILAENV2STAGE( 4, 'SSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )
         LWMIN = 2*N + LHTRD + LWTRD
         WORK( 1 )  = LWMIN

         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -11
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('SSYGV_2STAGE ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Form a Cholesky factorization of B.

      spotrf(UPLO, N, B, LDB, INFO );
      if ( INFO.NE.0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem and solve.

      ssygst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO );
      ssyev_2stage(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO );

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         NEIG = N
         IF( INFO.GT.0 ) NEIG = INFO - 1
         if ( ITYPE.EQ.1 .OR. ITYPE.EQ.2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N'
            } else {
               TRANS = 'T'
            }

            strsm('Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, B, LDB, A, LDA );

         } else if ( ITYPE.EQ.3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**T*y

            if ( UPPER ) {
               TRANS = 'T'
            } else {
               TRANS = 'N'
            }

            strmm('Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, B, LDB, A, LDA );
         }
      }

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      RETURN

      // End of SSYGV_2STAGE

      }
