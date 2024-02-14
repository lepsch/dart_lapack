      void chegvd(final int ITYPE, final int JOBZ, final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int W, final Array<double> WORK_, final int LWORK, final Array<int> RWORK_, final int LRWORK, final Array<int> IWORK_, final int LIWORK, final Box<int> INFO,) {
  final A = A_.dim();
  final B = B_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDA, LDB, LIWORK, LRWORK, LWORK, N;
      int                IWORK( * );
      double               RWORK( * ), W( * );
      Complex            A( LDA, * ), B( LDB, * ), WORK( * );
      // ..

      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      bool               LQUERY, UPPER, WANTZ;
      String             TRANS;
      int                LIOPT, LIWMIN, LOPT, LROPT, LRWMIN, LWMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHEEVD, CHEGST, CPOTRF, CTRMM, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 || LRWORK == -1 || LIWORK == -1 );

      INFO = 0;
      if ( N <= 1 ) {
         LWMIN = 1;
         LRWMIN = 1;
         LIWMIN = 1;
      } else if ( WANTZ ) {
         LWMIN = 2*N + N*N;
         LRWMIN = 1 + 5*N + 2*N*N;
         LIWMIN = 3 + 5*N;
      } else {
         LWMIN = N + 1;
         LRWMIN = N;
         LIWMIN = 1;
      }
      LOPT = LWMIN;
      LROPT = LRWMIN;
      LIOPT = LIWMIN;
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1;
      } else if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( UPPER || lsame( UPLO, 'L' ) ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      }

      if ( INFO == 0 ) {
         WORK[1] = SROUNDUP_LWORK(LOPT);
         RWORK[1] = LROPT;
         IWORK[1] = LIOPT;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -11;
         } else if ( LRWORK < LRWMIN && !LQUERY ) {
            INFO = -13;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -15;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CHEGVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Form a Cholesky factorization of B.

      cpotrf(UPLO, N, B, LDB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem and solve.

      chegst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO );
      cheevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO );
      LOPT = INT( max( REAL( LOPT ), double( WORK( 1 ) ) ) );
      LROPT = INT( max( REAL( LROPT ), double( RWORK( 1 ) ) ) );
      LIOPT = INT( max( REAL( LIOPT ), double( IWORK( 1 ) ) ) );

      if ( WANTZ && INFO == 0 ) {

         // Backtransform eigenvectors to the original problem.

         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N';
            } else {
               TRANS = 'C';
            }

            ctrsm('Left', UPLO, TRANS, 'Non-unit', N, N, CONE, B, LDB, A, LDA );

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**H *y

            if ( UPPER ) {
               TRANS = 'C';
            } else {
               TRANS = 'N';
            }

            ctrmm('Left', UPLO, TRANS, 'Non-unit', N, N, CONE, B, LDB, A, LDA );
         }
      }

      WORK[1] = SROUNDUP_LWORK(LOPT);
      RWORK[1] = LROPT;
      IWORK[1] = LIOPT;

      }
