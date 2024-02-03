      SUBROUTINE ZHPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDZ, M, N;
      double             ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      double             RWORK( * ), W( * );
      COMPLEX*16         AP( * ), WORK( * ), Z( LDZ, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, TEST, VALEIG, WANTZ;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDISP, INDIWK, INDRWK, INDTAU, INDWRK, ISCALE, ITMP1, J, JJ, NSPLIT;
      double             ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANHP;
      // EXTERNAL LSAME, DLAMCH, ZLANHP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSCAL, DSTEBZ, DSTERF, XERBLA, ZDSCAL, ZHPTRD, ZSTEIN, ZSTEQR, ZSWAP, ZUPGTR, ZUPMTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( LSAME( UPLO, 'L' ) .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE
         IF( VALEIG ) THEN
            IF( N.GT.0 .AND. VU.LE.VL ) INFO = -7
         ELSE IF( INDEIG ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) THEN
               INFO = -8
            ELSE IF( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) THEN
               INFO = -9
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) INFO = -14
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHPEVX', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      M = 0
      IF( N.EQ.0 ) RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = DBLE( AP( 1 ) )
         ELSE
            IF( VL.LT.DBLE( AP( 1 ) ) .AND. VU.GE.DBLE( AP( 1 ) ) ) THEN
               M = 1
               W( 1 ) = DBLE( AP( 1 ) )
            END IF
         END IF
         IF( WANTZ ) Z( 1, 1 ) = CONE
         RETURN
      END IF
*
      // Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
      // Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
      ABSTLL = ABSTOL
      IF( VALEIG ) THEN
         VLL = VL
         VUU = VU
      ELSE
         VLL = ZERO
         VUU = ZERO
      END IF
      ANRM = ZLANHP( 'M', UPLO, N, AP, RWORK )
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         CALL ZDSCAL( ( N*( N+1 ) ) / 2, SIGMA, AP, 1 )
         IF( ABSTOL.GT.0 ) ABSTLL = ABSTOL*SIGMA
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF
*
      // Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form.
*
      INDD = 1
      INDE = INDD + N
      INDRWK = INDE + N
      INDTAU = 1
      INDWRK = INDTAU + N
      CALL ZHPTRD( UPLO, N, AP, RWORK( INDD ), RWORK( INDE ), WORK( INDTAU ), IINFO )
*
      // If all eigenvalues are desired and ABSTOL is less than or equal
     t // o zero, then call DSTERF or ZUPGTR and ZSTEQR.  If this fails
      // for some eigenvalue, then try DSTEBZ.
*
      TEST = .FALSE.
      IF (INDEIG) THEN
         IF (IL.EQ.1 .AND. IU.EQ.N) THEN
            TEST = .TRUE.
         END IF
      END IF
      IF ((ALLEIG .OR. TEST) .AND. (ABSTOL.LE.ZERO)) THEN
         CALL DCOPY( N, RWORK( INDD ), 1, W, 1 )
         INDEE = INDRWK + 2*N
         IF( .NOT.WANTZ ) THEN
            CALL DCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
            CALL DSTERF( N, W, RWORK( INDEE ), INFO )
         ELSE
            CALL ZUPGTR( UPLO, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO )
            CALL DCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
            CALL ZSTEQR( JOBZ, N, W, RWORK( INDEE ), Z, LDZ, RWORK( INDRWK ), INFO )
            IF( INFO.EQ.0 ) THEN
               DO 10 I = 1, N
                  IFAIL( I ) = 0
   10          CONTINUE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 20
         END IF
         INFO = 0
      END IF
*
      // Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.
*
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDISP = 1 + N
      INDIWK = INDISP + N
      CALL DSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, RWORK( INDD ), RWORK( INDE ), M, NSPLIT, W, IWORK( 1 ), IWORK( INDISP ), RWORK( INDRWK ), IWORK( INDIWK ), INFO )
*
      IF( WANTZ ) THEN
         CALL ZSTEIN( N, RWORK( INDD ), RWORK( INDE ), M, W, IWORK( 1 ), IWORK( INDISP ), Z, LDZ, RWORK( INDRWK ), IWORK( INDIWK ), IFAIL, INFO )
*
         // Apply unitary matrix used in reduction to tridiagonal
         // form to eigenvectors returned by ZSTEIN.
*
         INDWRK = INDTAU + N
         CALL ZUPMTR( 'L', UPLO, 'N', N, M, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO )
      END IF
*
      // If matrix was scaled, then rescale eigenvalues appropriately.
*
   20 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = M
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.
*
      IF( WANTZ ) THEN
         DO 40 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 30 JJ = J + 1, M
               IF( W( JJ ).LT.TMP1 ) THEN
                  I = JJ
                  TMP1 = W( JJ )
               END IF
   30       CONTINUE
*
            IF( I.NE.0 ) THEN
               ITMP1 = IWORK( 1 + I-1 )
               W( I ) = W( J )
               IWORK( 1 + I-1 ) = IWORK( 1 + J-1 )
               W( J ) = TMP1
               IWORK( 1 + J-1 ) = ITMP1
               CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
               IF( INFO.NE.0 ) THEN
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               END IF
            END IF
   40    CONTINUE
      END IF
*
      RETURN
*
      // End of ZHPEVX
*
      END
