      void dsbgvd(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             AB( LDAB, * ), BB( LDBB, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER, WANTZ;
      String             VECT;
      int                IINFO, INDE, INDWK2, INDWRK, LIWMIN, LLWRK2, LWMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DPBSTF, DSBGST, DSBTRD, DSTEDC, DSTERF, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' );
      UPPER = LSAME( UPLO, 'U' );
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

      if ( !( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( UPPER || LSAME( UPLO, 'L' ) ) ) {
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
         WORK( 1 ) = LWMIN;
         IWORK( 1 ) = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -14;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -16;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DSBGVD', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Form a split Cholesky factorization of B.

      dpbstf(UPLO, N, KB, BB, LDBB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem.

      INDE = 1;
      INDWRK = INDE + N;
      INDWK2 = INDWRK + N*N;
      LLWRK2 = LWORK - INDWK2 + 1;
      dsbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK, IINFO );

      // Reduce to tridiagonal form.

      if ( WANTZ ) {
         VECT = 'U';
      } else {
         VECT = 'N';
      }
      dsbtrd(VECT, UPLO, N, KA, AB, LDAB, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), IINFO );

      // For eigenvalues only, call DSTERF. For eigenvectors, call SSTEDC.

      if ( !WANTZ ) {
         dsterf(N, W, WORK( INDE ), INFO );
      } else {
         dstedc('I', N, W, WORK( INDE ), WORK( INDWRK ), N, WORK( INDWK2 ), LLWRK2, IWORK, LIWORK, INFO );
         dgemm('N', 'N', N, N, N, ONE, Z, LDZ, WORK( INDWRK ), N, ZERO, WORK( INDWK2 ), N );
         dlacpy('A', N, N, WORK( INDWK2 ), N, Z, LDZ );
      }

      WORK( 1 ) = LWMIN;
      IWORK( 1 ) = LIWMIN;

      return;

      // End of DSBGVD

      }
