      SUBROUTINE ZHEGV_2STAGE( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, RWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDA, LDB, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             RWORK( * ), W( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      const              ONE = ( 1.0D+0, 0.0D+0 ) ;
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
      // EXTERNAL XERBLA, ZHEGST, ZPOTRF, ZTRMM, ZTRSM, ZHEEV_2STAGE
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
         KD    = ILAENV2STAGE( 1, 'ZHETRD_2STAGE', JOBZ, N, -1, -1, -1 )
         IB    = ILAENV2STAGE( 2, 'ZHETRD_2STAGE', JOBZ, N, KD, -1, -1 )
         LHTRD = ILAENV2STAGE( 3, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 )
         LWTRD = ILAENV2STAGE( 4, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 )
         LWMIN = N + LHTRD + LWTRD
         WORK( 1 )  = LWMIN

         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -11
         }
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZHEGV_2STAGE ', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Form a Cholesky factorization of B.

      CALL ZPOTRF( UPLO, N, B, LDB, INFO )
      if ( INFO.NE.0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem and solve.

      CALL ZHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      CALL ZHEEV_2STAGE( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         NEIG = N
         IF( INFO.GT.0 ) NEIG = INFO - 1
         if ( ITYPE.EQ.1 .OR. ITYPE.EQ.2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N'
            } else {
               TRANS = 'C'
            }

            CALL ZTRSM( 'Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, B, LDB, A, LDA )

         } else if ( ITYPE.EQ.3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**H *y

            if ( UPPER ) {
               TRANS = 'C'
            } else {
               TRANS = 'N'
            }

            CALL ZTRMM( 'Left', UPLO, TRANS, 'Non-unit', N, NEIG, ONE, B, LDB, A, LDA )
         }
      }

      WORK( 1 ) = LWMIN

      RETURN

      // End of ZHEGV_2STAGE

      }
