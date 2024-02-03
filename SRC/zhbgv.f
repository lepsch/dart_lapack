      SUBROUTINE ZHBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, LDZ, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KA, KB, LDAB, LDBB, LDZ, N;
      // ..
      // .. Array Arguments ..
      double             RWORK( * ), W( * );
      COMPLEX*16         AB( LDAB, * ), BB( LDBB, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               UPPER, WANTZ;
      String             VECT;
      int                IINFO, INDE, INDWRK;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSTERF, XERBLA, ZHBGST, ZHBTRD, ZPBSTF, ZSTEQR
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )

      INFO = 0
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( KA.LT.0 ) {
         INFO = -4
      } else if ( KB.LT.0 .OR. KB.GT.KA ) {
         INFO = -5
      } else if ( LDAB.LT.KA+1 ) {
         INFO = -7
      } else if ( LDBB.LT.KB+1 ) {
         INFO = -9
      } else if ( LDZ.LT.1 .OR. ( WANTZ && LDZ.LT.N ) ) {
         INFO = -12
      }
      if ( INFO != 0 ) {
         xerbla('ZHBGV', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Form a split Cholesky factorization of B.

      zpbstf(UPLO, N, KB, BB, LDBB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem.

      INDE = 1
      INDWRK = INDE + N
      zhbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK, RWORK( INDWRK ), IINFO );

      // Reduce to tridiagonal form.

      if ( WANTZ ) {
         VECT = 'U'
      } else {
         VECT = 'N'
      }
      zhbtrd(VECT, UPLO, N, KA, AB, LDAB, W, RWORK( INDE ), Z, LDZ, WORK, IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEQR.

      if ( .NOT.WANTZ ) {
         dsterf(N, W, RWORK( INDE ), INFO );
      } else {
         zsteqr(JOBZ, N, W, RWORK( INDE ), Z, LDZ, RWORK( INDWRK ), INFO );
      }
      RETURN

      // End of ZHBGV

      }
