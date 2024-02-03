      SUBROUTINE CHBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, KD, LDAB, LDZ, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * ), W( * )
      COMPLEX            AB( LDAB, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, WANTZ;
      int                IINFO, IMAX, INDE, INDRWK, ISCALE;
      REAL               ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHB, SLAMCH
      // EXTERNAL LSAME, CLANHB, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHBTRD, CLASCL, CSTEQR, SSCAL, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )

      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KD.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -9
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHBEV ', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      IF( N.EQ.1 ) THEN
         IF( LOWER ) THEN
            W( 1 ) = REAL( AB( 1, 1 ) )
         } else {
            W( 1 ) = REAL( AB( KD+1, 1 ) )
         END IF
         IF( WANTZ ) Z( 1, 1 ) = ONE
         RETURN
      END IF

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )

      // Scale matrix to allowable range, if necessary.

      ANRM = CLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            CALL CLASCL( 'B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         } else {
            CALL CLASCL( 'Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         END IF
      END IF

      // Call CHBTRD to reduce Hermitian band matrix to tridiagonal form.

      INDE = 1
      CALL CHBTRD( JOBZ, UPLO, N, KD, AB, LDAB, W, RWORK( INDE ), Z, LDZ, WORK, IINFO )

      // For eigenvalues only, call SSTERF.  For eigenvectors, call CSTEQR.

      IF( .NOT.WANTZ ) THEN
         CALL SSTERF( N, W, RWORK( INDE ), INFO )
      } else {
         INDRWK = INDE + N
         CALL CSTEQR( JOBZ, N, W, RWORK( INDE ), Z, LDZ, RWORK( INDRWK ), INFO )
      END IF

      // If matrix was scaled, then rescale eigenvalues appropriately.

      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         } else {
            IMAX = INFO - 1
         END IF
         CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF

      RETURN

      // End of CHBEV

      }
