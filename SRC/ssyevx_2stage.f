      SUBROUTINE SSYEVX_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDA, LDZ, LWORK, M, N;
      REAL               ABSTOL, VL, VU
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      REAL               A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LOWER, LQUERY, TEST, VALEIG, WANTZ;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, INDISP, INDIWO, INDTAU, INDWKN, INDWRK, ISCALE, ITMP1, J, JJ, LLWORK, LLWRKN, NSPLIT, LWMIN, LHTRD, LWTRD, KD, IB, INDHOUS;
      REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      REAL               SLAMCH, SLANSY, SROUNDUP_LWORK
      // EXTERNAL LSAME, SLAMCH, SLANSY, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLACPY, SORGTR, SORMTR, SSCAL, SSTEBZ, SSTEIN, SSTEQR, SSTERF, SSWAP, XERBLA, SSYTRD_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
      LQUERY = ( LWORK.EQ.-1 )

      INFO = 0
      IF( .NOT.( LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE
         IF( VALEIG ) THEN
            IF( N.GT.0 .AND. VU.LE.VL ) INFO = -8
         ELSE IF( INDEIG ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) THEN
               INFO = -9
            ELSE IF( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) THEN
               INFO = -10
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
            INFO = -15
         END IF
      END IF

      IF( INFO.EQ.0 ) THEN
         IF( N.LE.1 ) THEN
            LWMIN = 1
            WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         ELSE
            KD    = ILAENV2STAGE( 1, 'SSYTRD_2STAGE', JOBZ, N, -1, -1, -1 )             IB    = ILAENV2STAGE( 2, 'SSYTRD_2STAGE', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'SSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'SSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )
            LWMIN = MAX( 8*N, 3*N + LHTRD + LWTRD )
            WORK( 1 )  = LWMIN
         END IF

         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) INFO = -17
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSYEVX_2STAGE', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      M = 0
      IF( N.EQ.0 ) THEN
         RETURN
      END IF

      IF( N.EQ.1 ) THEN
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = A( 1, 1 )
         ELSE
            IF( VL.LT.A( 1, 1 ) .AND. VU.GE.A( 1, 1 ) ) THEN
               M = 1
               W( 1 ) = A( 1, 1 )
            END IF
         END IF
         IF( WANTZ ) Z( 1, 1 ) = ONE
         RETURN
      END IF

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS    = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )

      // Scale matrix to allowable range, if necessary.

      ISCALE = 0
      ABSTLL = ABSTOL
      IF( VALEIG ) THEN
         VLL = VL
         VUU = VU
      END IF
      ANRM = SLANSY( 'M', UPLO, N, A, LDA, WORK )
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            DO 10 J = 1, N
               CALL SSCAL( N-J+1, SIGMA, A( J, J ), 1 )
   10       CONTINUE
         ELSE
            DO 20 J = 1, N
               CALL SSCAL( J, SIGMA, A( 1, J ), 1 )
   20       CONTINUE
         END IF
         IF( ABSTOL.GT.0 ) ABSTLL = ABSTOL*SIGMA
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF

      // Call SSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form.

      INDTAU  = 1
      INDE    = INDTAU + N
      INDD    = INDE + N
      INDHOUS = INDD + N
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1

      CALL SSYTRD_2STAGE( JOBZ, UPLO, N, A, LDA, WORK( INDD ), WORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO )

      // If all eigenvalues are desired and ABSTOL is less than or equal to
      // zero, then call SSTERF or SORGTR and SSTEQR.  If this fails for
      // some eigenvalue, then try SSTEBZ.

      TEST = .FALSE.
      IF( INDEIG ) THEN
         IF( IL.EQ.1 .AND. IU.EQ.N ) THEN
            TEST = .TRUE.
         END IF
      END IF
      IF( ( ALLEIG .OR. TEST ) .AND. ( ABSTOL.LE.ZERO ) ) THEN
         CALL SCOPY( N, WORK( INDD ), 1, W, 1 )
         INDEE = INDWRK + 2*N
         IF( .NOT.WANTZ ) THEN
            CALL SCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL SSTERF( N, W, WORK( INDEE ), INFO )
         ELSE
            CALL SLACPY( 'A', N, N, A, LDA, Z, LDZ )
            CALL SORGTR( UPLO, N, Z, LDZ, WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO )
            CALL SCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL SSTEQR( JOBZ, N, W, WORK( INDEE ), Z, LDZ, WORK( INDWRK ), INFO )
            IF( INFO.EQ.0 ) THEN
               DO 30 I = 1, N
                  IFAIL( I ) = 0
   30          CONTINUE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 40
         END IF
         INFO = 0
      END IF

      // Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN.

      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWO = INDISP + N
      CALL SSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, WORK( INDD ), WORK( INDE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWRK ), IWORK( INDIWO ), INFO )

      IF( WANTZ ) THEN
         CALL SSTEIN( N, WORK( INDD ), WORK( INDE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO )

         // Apply orthogonal matrix used in reduction to tridiagonal
         // form to eigenvectors returned by SSTEIN.

         INDWKN = INDE
         LLWRKN = LWORK - INDWKN + 1
         CALL SORMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO )
      END IF

      // If matrix was scaled, then rescale eigenvalues appropriately.

   40 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = M
         ELSE
            IMAX = INFO - 1
         END IF
         CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF

      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.

      IF( WANTZ ) THEN
         DO 60 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 50 JJ = J + 1, M
               IF( W( JJ ).LT.TMP1 ) THEN
                  I = JJ
                  TMP1 = W( JJ )
               END IF
   50       CONTINUE

            IF( I.NE.0 ) THEN
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL SSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
               IF( INFO.NE.0 ) THEN
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               END IF
            END IF
   60    CONTINUE
      END IF

      // Set WORK(1) to optimal workspace size.

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)

      RETURN

      // End of SSYEVX_2STAGE

      }
