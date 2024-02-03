      SUBROUTINE CHBGVD( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KA, KB, LDAB, LDBB, LDZ, LIWORK, LRWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               RWORK( * ), W( * )
      COMPLEX            AB( LDAB, * ), BB( LDBB, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            CONE, CZERO
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER, WANTZ;
      String             VECT;
      int                IINFO, INDE, INDWK2, INDWRK, LIWMIN, LLRWK, LLWK2, LRWMIN, LWMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSTERF, XERBLA, CGEMM, CHBGST, CHBTRD, CLACPY, CPBSTF, CSTEDC
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 || LRWORK == -1 || LIWORK == -1 )

      INFO = 0
      if ( N <= 1 ) {
         LWMIN = 1+N
         LRWMIN = 1+N
         LIWMIN = 1
      } else if ( WANTZ ) {
         LWMIN = 2*N**2
         LRWMIN = 1 + 5*N + 2*N**2
         LIWMIN = 3 + 5*N
      } else {
         LWMIN = N
         LRWMIN = N
         LIWMIN = 1
      }
      if ( !( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( !( UPPER || LSAME( UPLO, 'L' ) ) ) {
         INFO = -2
      } else if ( N < 0 ) {
         INFO = -3
      } else if ( KA < 0 ) {
         INFO = -4
      } else if ( KB < 0 || KB > KA ) {
         INFO = -5
      } else if ( LDAB < KA+1 ) {
         INFO = -7
      } else if ( LDBB < KB+1 ) {
         INFO = -9
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -12
      }

      if ( INFO == 0 ) {
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -14
         } else if ( LRWORK < LRWMIN && !LQUERY ) {
            INFO = -16
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -18
         }
      }

      if ( INFO != 0 ) {
         xerbla('CHBGVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Form a split Cholesky factorization of B.

      cpbstf(UPLO, N, KB, BB, LDBB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem.

      INDE = 1
      INDWRK = INDE + N
      INDWK2 = 1 + N*N
      LLWK2 = LWORK - INDWK2 + 2
      LLRWK = LRWORK - INDWRK + 2
      chbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK, RWORK, IINFO );

      // Reduce Hermitian band matrix to tridiagonal form.

      if ( WANTZ ) {
         VECT = 'U'
      } else {
         VECT = 'N'
      }
      chbtrd(VECT, UPLO, N, KA, AB, LDAB, W, RWORK( INDE ), Z, LDZ, WORK, IINFO );

      // For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEDC.

      if ( !WANTZ ) {
         ssterf(N, W, RWORK( INDE ), INFO );
      } else {
         cstedc('I', N, W, RWORK( INDE ), WORK, N, WORK( INDWK2 ), LLWK2, RWORK( INDWRK ), LLRWK, IWORK, LIWORK, INFO );
         cgemm('N', 'N', N, N, N, CONE, Z, LDZ, WORK, N, CZERO, WORK( INDWK2 ), N );
         clacpy('A', N, N, WORK( INDWK2 ), N, Z, LDZ );
      }

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
      RETURN

      // End of CHBGVD

      }
