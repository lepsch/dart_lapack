      SUBROUTINE CHBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, LDZ, WORK, RWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KA, KB, LDAB, LDBB, LDZ, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * ), W( * )
      COMPLEX            AB( LDAB, * ), BB( LDBB, * ), WORK( * ), Z( LDZ, * )
      // ..
*
*  =====================================================================
*
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
      // EXTERNAL CHBGST, CHBTRD, CPBSTF, CSTEQR, SSTERF, XERBLA
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KA.LT.0 ) THEN
         INFO = -4
      ELSE IF( KB.LT.0 .OR. KB.GT.KA ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.KA+1 ) THEN
         INFO = -7
      ELSE IF( LDBB.LT.KB+1 ) THEN
         INFO = -9
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHBGV', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      // Form a split Cholesky factorization of B.
*
      CALL CPBSTF( UPLO, N, KB, BB, LDBB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
*
      // Transform problem to standard eigenvalue problem.
*
      INDE = 1
      INDWRK = INDE + N
      CALL CHBGST( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK, RWORK( INDWRK ), IINFO )
*
      // Reduce to tridiagonal form.
*
      IF( WANTZ ) THEN
         VECT = 'U'
      ELSE
         VECT = 'N'
      END IF
      CALL CHBTRD( VECT, UPLO, N, KA, AB, LDAB, W, RWORK( INDE ), Z, LDZ, WORK, IINFO )
*
      // For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL SSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL CSTEQR( JOBZ, N, W, RWORK( INDE ), Z, LDZ, RWORK( INDWRK ), INFO )
      END IF
      RETURN
*
      // End of CHBGV
*
      END
