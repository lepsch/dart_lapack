      void chegvx(final int ITYPE, final int JOBZ, final int RANGE, final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int VL, final int VU, final int IL, final int IU, final int ABSTOL, final int M, final int W, final Matrix<double> Z_, final int LDZ, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final Array<int> IWORK_, final int IFAIL, final Box<int> INFO,) {
  final A = A_.dim();
  final B = B_.dim();
  final Z = Z_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N;
      double               ABSTOL, VL, VU;
      int                IFAIL( * ), IWORK( * );
      double               RWORK( * ), W( * );
      Complex            A( LDA, * ), B( LDB, * ), WORK( * ), Z( LDZ, * );
      // ..

      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      bool               ALLEIG, INDEIG, LQUERY, UPPER, VALEIG, WANTZ;
      String             TRANS;
      int                LWKOPT, NB;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHEEVX, CHEGST, CPOTRF, CTRMM, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      UPPER = lsame( UPLO, 'U' );
      ALLEIG = lsame( RANGE, 'A' );
      VALEIG = lsame( RANGE, 'V' );
      INDEIG = lsame( RANGE, 'I' );
      LQUERY = ( LWORK == -1 );

      INFO = 0;
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1;
      } else if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -3;
      } else if ( !( UPPER || lsame( UPLO, 'L' ) ) ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else {
         if ( VALEIG ) {
            if (N > 0 && VU <= VL) INFO = -11;
         } else if ( INDEIG ) {
            if ( IL < 1 || IL > max( 1, N ) ) {
               INFO = -12;
            } else if ( IU < min( N, IL ) || IU > N ) {
               INFO = -13;
            }
         }
      }
      if (INFO == 0) {
         if (LDZ < 1 || (WANTZ && LDZ < N)) {
            INFO = -18;
         }
      }

      if ( INFO == 0 ) {
         NB = ilaenv( 1, 'CHETRD', UPLO, N, -1, -1, -1 );
         LWKOPT = max( 1, ( NB + 1 )*N );
         WORK[1] = SROUNDUP_LWORK(LWKOPT);

         if ( LWORK < max( 1, 2*N ) && !LQUERY ) {
            INFO = -20;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CHEGVX', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      M = 0;
      if ( N == 0 ) {
         return;
      }

      // Form a Cholesky factorization of B.

      cpotrf(UPLO, N, B, LDB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem and solve.

      chegst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO );
      cheevx(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL, INFO );

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         if (INFO > 0) M = INFO - 1;
         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**H*y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N';
            } else {
               TRANS = 'C';
            }

            ctrsm('Left', UPLO, TRANS, 'Non-unit', N, M, CONE, B, LDB, Z, LDZ );

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**H*y

            if ( UPPER ) {
               TRANS = 'C';
            } else {
               TRANS = 'N';
            }

            ctrmm('Left', UPLO, TRANS, 'Non-unit', N, M, CONE, B, LDB, Z, LDZ );
         }
      }

      // Set WORK(1) to optimal complex workspace size.

      WORK[1] = SROUNDUP_LWORK(LWKOPT);

      }
