      SUBROUTINE SSBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, LDZ, WORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KA, KB, LDAB, LDBB, LDZ, N;
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * ), BB( LDBB, * ), W( * ), WORK( * ), Z( LDZ, * )
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
      // EXTERNAL SPBSTF, SSBGST, SSBTRD, SSTEQR, SSTERF, XERBLA
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
      } else if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
         INFO = -12
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SSBGV', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Form a split Cholesky factorization of B.

      CALL SPBSTF( UPLO, N, KB, BB, LDBB, INFO )
      if ( INFO.NE.0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem.

      INDE = 1
      INDWRK = INDE + N
      CALL SSBGST( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK( INDWRK ), IINFO )

      // Reduce to tridiagonal form.

      if ( WANTZ ) {
         VECT = 'U'
      } else {
         VECT = 'N'
      }
      CALL SSBTRD( VECT, UPLO, N, KA, AB, LDAB, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), IINFO )

      // For eigenvalues only, call SSTERF.  For eigenvectors, call SSTEQR.

      if ( .NOT.WANTZ ) {
         CALL SSTERF( N, W, WORK( INDE ), INFO )
      } else {
         CALL SSTEQR( JOBZ, N, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), INFO )
      }
      RETURN

      // End of SSBGV

      }
