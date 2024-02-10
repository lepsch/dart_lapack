      void zhpgvd(final int ITYPE, final int JOBZ, final int UPLO, final int N, final int AP, final int BP, final int W, final Matrix<double> Z, final int LDZ, final Array<double> WORK, final int LWORK, final Array<int> RWORK, final int LRWORK, final Array<int> IWORK, final int LIWORK, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDZ, LIWORK, LRWORK, LWORK, N;
      int                IWORK( * );
      double             RWORK( * ), W( * );
      Complex         AP( * ), BP( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY, UPPER, WANTZ;
      String             TRANS;
      int                J, LIWMIN, LRWMIN, LWMIN, NEIG;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZHPEVD, ZHPGST, ZPPTRF, ZTPMV, ZTPSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 || LRWORK == -1 || LIWORK == -1 );

      INFO = 0;
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1;
      } else if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( UPPER || lsame( UPLO, 'L' ) ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -9;
      }

      if ( INFO == 0 ) {
         if ( N <= 1 ) {
            LWMIN = 1;
            LIWMIN = 1;
            LRWMIN = 1;
         } else {
            if ( WANTZ ) {
               LWMIN = 2*N;
               LRWMIN = 1 + 5*N + 2*N**2;
               LIWMIN = 3 + 5*N;
            } else {
               LWMIN = N;
               LRWMIN = N;
               LIWMIN = 1;
            }
         }

         WORK[1] = LWMIN;
         RWORK[1] = LRWMIN;
         IWORK[1] = LIWMIN;
         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -11;
         } else if ( LRWORK < LRWMIN && !LQUERY ) {
            INFO = -13;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -15;
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZHPGVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Form a Cholesky factorization of B.

      zpptrf(UPLO, N, BP, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem and solve.

      zhpgst(ITYPE, UPLO, N, AP, BP, INFO );
      zhpevd(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO );
      LWMIN = INT( max( LWMIN.toDouble(), (WORK( 1 )).toDouble() ) );
      LRWMIN = INT( max( LRWMIN.toDouble(), (RWORK( 1 )).toDouble() ) );
      LIWMIN = INT( max( LIWMIN.toDouble(), (IWORK( 1 )).toDouble() ) );

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         NEIG = N;
         if (INFO > 0) NEIG = INFO - 1;
         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N';
            } else {
               TRANS = 'C';
            }

            for (J = 1; J <= NEIG; J++) { // 10
               ztpsv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
            } // 10

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**H *y

            if ( UPPER ) {
               TRANS = 'C';
            } else {
               TRANS = 'N';
            }

            for (J = 1; J <= NEIG; J++) { // 20
               ztpmv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
            } // 20
         }
      }

      WORK[1] = LWMIN;
      RWORK[1] = LRWMIN;
      IWORK[1] = LIWMIN;
      }
