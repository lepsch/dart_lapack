      void ssbgvd(final int JOBZ, final int UPLO, final int N, final int KA, final int KB, final Matrix<double> AB_, final int LDAB, final Matrix<double> BB_, final int LDBB, final int W, final Matrix<double> Z_, final int LDZ, final Array<double> WORK_, final int LWORK, final Array<int> IWORK_, final int LIWORK, final Box<int> INFO,) {
  final AB = AB_.dim();
  final BB = BB_.dim();
  final Z = Z_.dim();
  final WORK = WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LWORK, N;
      int                IWORK( * );
      double               AB( LDAB, * ), BB( LDBB, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      bool               LQUERY, UPPER, WANTZ;
      String             VECT;
      int                IINFO, INDE, INDWK2, INDWRK, LIWMIN, LLWRK2, LWMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL lsame, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SPBSTF, SSBGST, SSBTRD, SSTEDC, SSTERF, XERBLA

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      UPPER = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 || LIWORK == -1 );

      INFO = 0;
      if ( N <= 1 ) {
         LIWMIN = 1;
         LWMIN = 1;
      } else if ( WANTZ ) {
         LIWMIN = 3 + 5*N;
         LWMIN = 1 + 5*N + 2*N**2;
      } else {
         LIWMIN = 1;
         LWMIN = 2*N;
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
         WORK[1] = SROUNDUP_LWORK(LWMIN);
         IWORK[1] = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -14;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -16;
         }
      }

      if ( INFO != 0 ) {
         xerbla('SSBGVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Form a split Cholesky factorization of B.

      spbstf(UPLO, N, KB, BB, LDBB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem.

      INDE = 1;
      INDWRK = INDE + N;
      INDWK2 = INDWRK + N*N;
      LLWRK2 = LWORK - INDWK2 + 1;
      ssbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK, IINFO );

      // Reduce to tridiagonal form.

      if ( WANTZ ) {
         VECT = 'U';
      } else {
         VECT = 'N';
      }
      ssbtrd(VECT, UPLO, N, KA, AB, LDAB, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), IINFO );

      // For eigenvalues only, call SSTERF. For eigenvectors, call SSTEDC.

      if ( !WANTZ ) {
         ssterf(N, W, WORK( INDE ), INFO );
      } else {
         sstedc('I', N, W, WORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO );
         sgemm('N', 'N', N, N, N, ONE, Z, LDZ, WORK( INDWRK ), N, ZERO, WORK( INDWK2 ), N );
         slacpy('A', N, N, WORK( INDWK2 ), N, Z, LDZ );
      }

      WORK[1] = SROUNDUP_LWORK(LWMIN);
      IWORK[1] = LIWMIN;

      }
