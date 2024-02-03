      SUBROUTINE ZHEGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDA, LDB, LIWORK, LRWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             RWORK( * ), W( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         CONE
      const              CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER, WANTZ;
      String             TRANS;
      int                LIOPT, LIWMIN, LOPT, LROPT, LRWMIN, LWMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZHEEVD, ZHEGST, ZPOTRF, ZTRMM, ZTRSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 || LRWORK == -1 || LIWORK == -1 )

      INFO = 0
      if ( N <= 1 ) {
         LWMIN = 1
         LRWMIN = 1
         LIWMIN = 1
      } else if ( WANTZ ) {
         LWMIN = 2*N + N*N
         LRWMIN = 1 + 5*N + 2*N*N
         LIWMIN = 3 + 5*N
      } else {
         LWMIN = N + 1
         LRWMIN = N
         LIWMIN = 1
      }
      LOPT = LWMIN
      LROPT = LRWMIN
      LIOPT = LIWMIN
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1
      } else if ( .NOT.( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
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
         WORK( 1 ) = LOPT
         RWORK( 1 ) = LROPT
         IWORK( 1 ) = LIOPT

         if ( LWORK < LWMIN && .NOT.LQUERY ) {
            INFO = -11
         } else if ( LRWORK < LRWMIN && .NOT.LQUERY ) {
            INFO = -13
         } else if ( LIWORK < LIWMIN && .NOT.LQUERY ) {
            INFO = -15
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZHEGVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Form a Cholesky factorization of B.

      zpotrf(UPLO, N, B, LDB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem and solve.

      zhegst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO );
      zheevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO );
      LOPT = INT( MAX( DBLE( LOPT ), DBLE( WORK( 1 ) ) ) )
      LROPT = INT( MAX( DBLE( LROPT ), DBLE( RWORK( 1 ) ) ) )
      LIOPT = INT( MAX( DBLE( LIOPT ), DBLE( IWORK( 1 ) ) ) )

      if ( WANTZ && INFO == 0 ) {

         // Backtransform eigenvectors to the original problem.

         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N'
            } else {
               TRANS = 'C'
            }

            ztrsm('Left', UPLO, TRANS, 'Non-unit', N, N, CONE, B, LDB, A, LDA );

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**H *y

            if ( UPPER ) {
               TRANS = 'C'
            } else {
               TRANS = 'N'
            }

            ztrmm('Left', UPLO, TRANS, 'Non-unit', N, N, CONE, B, LDB, A, LDA );
         }
      }

      WORK( 1 ) = LOPT
      RWORK( 1 ) = LROPT
      IWORK( 1 ) = LIOPT

      RETURN

      // End of ZHEGVD

      }
