      void zhbgvd(JOBZ, UPLO, N, KA, KB, final Matrix<double> AB, final int LDAB, final Matrix<double> BB, final int LDBB, W, final Matrix<double> Z, final int LDZ, final Array<double> WORK, final int LWORK, final Array<int> RWORK, final int LRWORK, final Array<int> IWORK, final int LIWORK, Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LRWORK, LWORK, N;
      int                IWORK( * );
      double             RWORK( * ), W( * );
      Complex         AB( LDAB, * ), BB( LDBB, * ), WORK( * ), Z( LDZ, * );
      // ..

      Complex         CONE, CZERO;
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      bool               LQUERY, UPPER, WANTZ;
      String             VECT;
      int                IINFO, INDE, INDWK2, INDWRK, LIWMIN, LLRWK, LLWK2, LRWMIN, LWMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSTERF, XERBLA, ZGEMM, ZHBGST, ZHBTRD, ZLACPY, ZPBSTF, ZSTEDC

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 || LRWORK == -1 || LIWORK == -1 );

      INFO = 0;
      if ( N <= 1 ) {
         LWMIN = 1+N;
         LRWMIN = 1+N;
         LIWMIN = 1;
      } else if ( WANTZ ) {
         LWMIN = 2*N**2;
         LRWMIN = 1 + 5*N + 2*N**2;
         LIWMIN = 3 + 5*N;
      } else {
         LWMIN = N;
         LRWMIN = N;
         LIWMIN = 1;
      }
      if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( UPPER || lsame( UPLO, 'L' ) ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( KA < 0 ) {
         INFO = -4;
      } else if ( KB < 0 || KB > KA ) {
         INFO = -5;
      } else if ( LDAB < KA+1 ) {
         INFO = -7;
      } else if ( LDBB < KB+1 ) {
         INFO = -9;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -12;
      }

      if ( INFO == 0 ) {
         WORK[1] = LWMIN;
         RWORK[1] = LRWMIN;
         IWORK[1] = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -14;
         } else if ( LRWORK < LRWMIN && !LQUERY ) {
            INFO = -16;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -18;
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZHBGVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Form a split Cholesky factorization of B.

      zpbstf(UPLO, N, KB, BB, LDBB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem.

      INDE = 1;
      INDWRK = INDE + N;
      INDWK2 = 1 + N*N;
      LLWK2 = LWORK - INDWK2 + 2;
      LLRWK = LRWORK - INDWRK + 2;
      zhbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK, RWORK, IINFO );

      // Reduce Hermitian band matrix to tridiagonal form.

      if ( WANTZ ) {
         VECT = 'U';
      } else {
         VECT = 'N';
      }
      zhbtrd(VECT, UPLO, N, KA, AB, LDAB, W, RWORK( INDE ), Z, LDZ, WORK, IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEDC.

      if ( !WANTZ ) {
         dsterf(N, W, RWORK( INDE ), INFO );
      } else {
         zstedc('I', N, W, RWORK( INDE ), WORK, N, WORK( INDWK2 ), LLWK2, RWORK( INDWRK ), LLRWK, IWORK, LIWORK, INFO );
         zgemm('N', 'N', N, N, N, CONE, Z, LDZ, WORK, N, CZERO, WORK( INDWK2 ), N );
         zlacpy('A', N, N, WORK( INDWK2 ), N, Z, LDZ );
      }

      WORK[1] = LWMIN;
      RWORK[1] = LRWMIN;
      IWORK[1] = LIWMIN;
      }
