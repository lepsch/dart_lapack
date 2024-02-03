      SUBROUTINE DSBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, LDZ, WORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KA, KB, LDAB, LDBB, LDZ, N;
      // ..
      // .. Array Arguments ..
      double             AB( LDAB, * ), BB( LDBB, * ), W( * ), WORK( * ), Z( LDZ, * );
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
      // EXTERNAL DPBSTF, DSBGST, DSBTRD, DSTEQR, DSTERF, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )

      INFO = 0
      if ( .NOT.( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( UPPER || LSAME( UPLO, 'L' ) ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( KA.LT.0 ) {
         INFO = -4
      } else if ( KB.LT.0 || KB.GT.KA ) {
         INFO = -5
      } else if ( LDAB.LT.KA+1 ) {
         INFO = -7
      } else if ( LDBB.LT.KB+1 ) {
         INFO = -9
      } else if ( LDZ.LT.1 || ( WANTZ && LDZ.LT.N ) ) {
         INFO = -12
      }
      if ( INFO != 0 ) {
         xerbla('DSBGV', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Form a split Cholesky factorization of B.

      dpbstf(UPLO, N, KB, BB, LDBB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem.

      INDE = 1
      INDWRK = INDE + N
      dsbgst(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK( INDWRK ), IINFO );

      // Reduce to tridiagonal form.

      if ( WANTZ ) {
         VECT = 'U'
      } else {
         VECT = 'N'
      }
      dsbtrd(VECT, UPLO, N, KA, AB, LDAB, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), IINFO );

      // For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR.

      if ( .NOT.WANTZ ) {
         dsterf(N, W, WORK( INDE ), INFO );
      } else {
         dsteqr(JOBZ, N, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), INFO );
      }
      RETURN

      // End of DSBGV

      }
