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
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ), CZERO = ( 0.0E+0, 0.0E+0 ) )
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
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

      INFO = 0
      IF( N.LE.1 ) THEN
         LWMIN = 1+N
         LRWMIN = 1+N
         LIWMIN = 1
      ELSE IF( WANTZ ) THEN
         LWMIN = 2*N**2
         LRWMIN = 1 + 5*N + 2*N**2
         LIWMIN = 3 + 5*N
      ELSE
         LWMIN = N
         LRWMIN = N
         LIWMIN = 1
      END IF
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

      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN

         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -14
         ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -16
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -18
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHBGVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Form a split Cholesky factorization of B.

      CALL CPBSTF( UPLO, N, KB, BB, LDBB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF

      // Transform problem to standard eigenvalue problem.

      INDE = 1
      INDWRK = INDE + N
      INDWK2 = 1 + N*N
      LLWK2 = LWORK - INDWK2 + 2
      LLRWK = LRWORK - INDWRK + 2
      CALL CHBGST( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, Z, LDZ, WORK, RWORK, IINFO )

      // Reduce Hermitian band matrix to tridiagonal form.

      IF( WANTZ ) THEN
         VECT = 'U'
      ELSE
         VECT = 'N'
      END IF
      CALL CHBTRD( VECT, UPLO, N, KA, AB, LDAB, W, RWORK( INDE ), Z, LDZ, WORK, IINFO )

      // For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEDC.

      IF( .NOT.WANTZ ) THEN
         CALL SSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL CSTEDC( 'I', N, W, RWORK( INDE ), WORK, N, WORK( INDWK2 ), LLWK2, RWORK( INDWRK ), LLRWK, IWORK, LIWORK, INFO )
         CALL CGEMM( 'N', 'N', N, N, N, CONE, Z, LDZ, WORK, N, CZERO, WORK( INDWK2 ), N )
         CALL CLACPY( 'A', N, N, WORK( INDWK2 ), N, Z, LDZ )
      END IF

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
      RETURN

      // End of CHBGVD

      END
