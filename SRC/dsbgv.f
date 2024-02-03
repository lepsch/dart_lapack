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
         CALL XERBLA( 'DSBGV', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Form a split Cholesky factorization of B.

      CALL DPBSTF( UPLO, N, KB, BB, LDBB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF

      // Transform problem to standard eigenvalue problem.

      INDE = 1
      INDWRK = INDE + N
      CALL DSBGST( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK( INDWRK ), IINFO )

      // Reduce to tridiagonal form.

      IF( WANTZ ) THEN
         VECT = 'U'
      } else {
         VECT = 'N'
      END IF
      CALL DSBTRD( VECT, UPLO, N, KA, AB, LDAB, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), IINFO )

      // For eigenvalues only, call DSTERF.  For eigenvectors, call SSTEQR.

      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      } else {
         CALL DSTEQR( JOBZ, N, W, WORK( INDE ), Z, LDZ, WORK( INDWRK ), INFO )
      END IF
      RETURN

      // End of DSBGV

      }
