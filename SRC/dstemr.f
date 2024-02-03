      SUBROUTINE DSTEMR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC, WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             JOBZ, RANGE;
      bool               TRYRAC;
      int                IL, INFO, IU, LDZ, NZC, LIWORK, LWORK, M, N
      DOUBLE PRECISION VL, VU
*     ..
*     .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
      DOUBLE PRECISION   Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, FOUR, MINRGP
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0, MINRGP = 1.0D-3 )
*     ..
*     .. Local Scalars ..
      bool               ALLEIG, INDEIG, LQUERY, VALEIG, WANTZ, ZQUERY, LAESWAP       int                I, IBEGIN, IEND, IFIRST, IIL, IINDBL, IINDW, IINDWK, IINFO, IINSPL, IIU, ILAST, IN, INDD, INDE2, INDERR, INDGP, INDGRS, INDWRK, ITMP, ITMP2, J, JBLK, JJ, LIWMIN, LWMIN, NSPLIT, NZCMIN, OFFSET, WBEGIN, WEND;
      DOUBLE PRECISION   BIGNUM, CS, EPS, PIVMIN, R1, R2, RMAX, RMIN, RTOL1, RTOL2, SAFMIN, SCALE, SMLNUM, SN, THRESH, TMP, TNRM, WL, WU
*     ..
*     ..
*     .. External Functions ..
      bool               LSAME;
      DOUBLE PRECISION   DLAMCH, DLANST
      EXTERNAL           LSAME, DLAMCH, DLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLAE2, DLAEV2, DLARRC, DLARRE, DLARRJ, DLARRR, DLARRV, DLASRT, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT


*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
*
      LQUERY = ( ( LWORK.EQ.-1 ).OR.( LIWORK.EQ.-1 ) )
      ZQUERY = ( NZC.EQ.-1 )
      LAESWAP = .FALSE.

*     DSTEMR needs WORK of size 6*N, IWORK of size 3*N.
*     In addition, DLARRE needs WORK of size 6*N, IWORK of size 5*N.
*     Furthermore, DLARRV needs WORK of size 12*N, IWORK of size 7*N.
      IF( WANTZ ) THEN
         LWMIN = 18*N
         LIWMIN = 10*N
      ELSE
*        need less workspace if only the eigenvalues are wanted
         LWMIN = 12*N
         LIWMIN = 8*N
      ENDIF

      WL = ZERO
      WU = ZERO
      IIL = 0
      IIU = 0
      NSPLIT = 0

      IF( VALEIG ) THEN
*        We do not reference VL, VU in the cases RANGE = 'I','A'
*        The interval (WL, WU] contains all the wanted eigenvalues.
*        It is either given by the user or computed in DLARRE.
         WL = VL
         WU = VU
      ELSEIF( INDEIG ) THEN
*        We do not reference IL, IU in the cases RANGE = 'V','A'
         IIL = IL
         IIU = IU
      ENDIF
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( VALEIG .AND. N.GT.0 .AND. WU.LE.WL ) THEN
         INFO = -7
      ELSE IF( INDEIG .AND. ( IIL.LT.1 .OR. IIL.GT.N ) ) THEN
         INFO = -8
      ELSE IF( INDEIG .AND. ( IIU.LT.IIL .OR. IIU.GT.N ) ) THEN
         INFO = -9
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -13
      ELSE IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -17
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -19
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
*
         IF( WANTZ .AND. ALLEIG ) THEN
            NZCMIN = N
         ELSE IF( WANTZ .AND. VALEIG ) THEN
            CALL DLARRC( 'T', N, VL, VU, D, E, SAFMIN, NZCMIN, ITMP, ITMP2, INFO )
         ELSE IF( WANTZ .AND. INDEIG ) THEN
            NZCMIN = IIU-IIL+1
         ELSE
*           WANTZ .EQ. FALSE.
            NZCMIN = 0
         ENDIF
         IF( ZQUERY .AND. INFO.EQ.0 ) THEN
            Z( 1,1 ) = NZCMIN
         ELSE IF( NZC.LT.NZCMIN .AND. .NOT.ZQUERY ) THEN
            INFO = -14
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
*
         CALL XERBLA( 'DSTEMR', -INFO )
*
         RETURN
      ELSE IF( LQUERY .OR. ZQUERY ) THEN
         RETURN
      END IF
*
*     Handle N = 0, 1, and 2 cases immediately
*
      M = 0
      IF( N.EQ.0 ) RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = D( 1 )
         ELSE
            IF( WL.LT.D( 1 ) .AND. WU.GE.D( 1 ) ) THEN
               M = 1
               W( 1 ) = D( 1 )
            END IF
         END IF
         IF( WANTZ.AND.(.NOT.ZQUERY) ) THEN
            Z( 1, 1 ) = ONE
            ISUPPZ(1) = 1
            ISUPPZ(2) = 1
         END IF
         RETURN
      END IF
*
      IF( N.EQ.2 ) THEN
         IF( .NOT.WANTZ ) THEN
            CALL DLAE2( D(1), E(1), D(2), R1, R2 )
         ELSE IF( WANTZ.AND.(.NOT.ZQUERY) ) THEN
            CALL DLAEV2( D(1), E(1), D(2), R1, R2, CS, SN )
         END IF
*        D/S/LAE2 and D/S/LAEV2 outputs satisfy |R1| >= |R2|. However,
*        the following code requires R1 >= R2. Hence, we correct
*        the order of R1, R2, CS, SN if R1 < R2 before further processing.
         IF( R1.LT.R2 ) THEN
            E(2) = R1
            R1 = R2
            R2 = E(2)
            LAESWAP = .TRUE.
         ENDIF
         IF( ALLEIG.OR. (VALEIG.AND.(R2.GT.WL).AND. (R2.LE.WU)).OR. (INDEIG.AND.(IIL.EQ.1)) ) THEN
            M = M+1
            W( M ) = R2
            IF( WANTZ.AND.(.NOT.ZQUERY) ) THEN
               IF( LAESWAP ) THEN
                  Z( 1, M ) = CS
                  Z( 2, M ) = SN
               ELSE
                  Z( 1, M ) = -SN
                  Z( 2, M ) = CS
               ENDIF
*              Note: At most one of SN and CS can be zero.
               IF (SN.NE.ZERO) THEN
                  IF (CS.NE.ZERO) THEN
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 2
                  ELSE
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 1
                  END IF
               ELSE
                  ISUPPZ(2*M-1) = 2
                  ISUPPZ(2*M) = 2
               END IF
            ENDIF
         ENDIF
         IF( ALLEIG.OR. (VALEIG.AND.(R1.GT.WL).AND. (R1.LE.WU)).OR. (INDEIG.AND.(IIU.EQ.2)) ) THEN
            M = M+1
            W( M ) = R1
            IF( WANTZ.AND.(.NOT.ZQUERY) ) THEN
               IF( LAESWAP ) THEN
                  Z( 1, M ) = -SN
                  Z( 2, M ) = CS
               ELSE
                  Z( 1, M ) = CS
                  Z( 2, M ) = SN
               ENDIF
*              Note: At most one of SN and CS can be zero.
               IF (SN.NE.ZERO) THEN
                  IF (CS.NE.ZERO) THEN
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 2
                  ELSE
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 1
                  END IF
               ELSE
                  ISUPPZ(2*M-1) = 2
                  ISUPPZ(2*M) = 2
               END IF
            ENDIF
         ENDIF

      ELSE

*     Continue with general N

         INDGRS = 1
         INDERR = 2*N + 1
         INDGP = 3*N + 1
         INDD = 4*N + 1
         INDE2 = 5*N + 1
         INDWRK = 6*N + 1
*
         IINSPL = 1
         IINDBL = N + 1
         IINDW = 2*N + 1
         IINDWK = 3*N + 1
*
*        Scale matrix to allowable range, if necessary.
*        The allowable range is related to the PIVMIN parameter; see the
*        comments in DLARRD.  The preference for scaling small values
*        up is heuristic; we expect users' matrices not to be close to the
*        RMAX threshold.
*
         SCALE = ONE
         TNRM = DLANST( 'M', N, D, E )
         IF( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) THEN
            SCALE = RMIN / TNRM
         ELSE IF( TNRM.GT.RMAX ) THEN
            SCALE = RMAX / TNRM
         END IF
         IF( SCALE.NE.ONE ) THEN
            CALL DSCAL( N, SCALE, D, 1 )
            CALL DSCAL( N-1, SCALE, E, 1 )
            TNRM = TNRM*SCALE
            IF( VALEIG ) THEN
*              If eigenvalues in interval have to be found,
*              scale (WL, WU] accordingly
               WL = WL*SCALE
               WU = WU*SCALE
            ENDIF
         END IF
*
*        Compute the desired eigenvalues of the tridiagonal after splitting
*        into smaller subblocks if the corresponding off-diagonal elements
*        are small
*        THRESH is the splitting parameter for DLARRE
*        A negative THRESH forces the old splitting criterion based on the
*        size of the off-diagonal. A positive THRESH switches to splitting
*        which preserves relative accuracy.
*
         IF( TRYRAC ) THEN
*           Test whether the matrix warrants the more expensive relative approach.
            CALL DLARRR( N, D, E, IINFO )
         ELSE
*           The user does not care about relative accurately eigenvalues
            IINFO = -1
         ENDIF
*        Set the splitting criterion
         IF (IINFO.EQ.0) THEN
            THRESH = EPS
         ELSE
            THRESH = -EPS
*           relative accuracy is desired but T does not guarantee it
            TRYRAC = .FALSE.
         ENDIF
*
         IF( TRYRAC ) THEN
*           Copy original diagonal, needed to guarantee relative accuracy
            CALL DCOPY(N,D,1,WORK(INDD),1)
         ENDIF
*        Store the squares of the offdiagonal values of T
         DO 5 J = 1, N-1
            WORK( INDE2+J-1 ) = E(J)**2
 5       CONTINUE

*        Set the tolerance parameters for bisection
         IF( .NOT.WANTZ ) THEN
*           DLARRE computes the eigenvalues to full precision.
            RTOL1 = FOUR * EPS
            RTOL2 = FOUR * EPS
         ELSE
*           DLARRE computes the eigenvalues to less than full precision.
*           DLARRV will refine the eigenvalue approximations, and we can
*           need less accurate initial bisection in DLARRE.
*           Note: these settings do only affect the subset case and DLARRE
            RTOL1 = SQRT(EPS)
            RTOL2 = MAX( SQRT(EPS)*5.0D-3, FOUR * EPS )
         ENDIF
         CALL DLARRE( RANGE, N, WL, WU, IIL, IIU, D, E, WORK(INDE2), RTOL1, RTOL2, THRESH, NSPLIT, IWORK( IINSPL ), M, W, WORK( INDERR ), WORK( INDGP ), IWORK( IINDBL ), IWORK( IINDW ), WORK( INDGRS ), PIVMIN, WORK( INDWRK ), IWORK( IINDWK ), IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 10 + ABS( IINFO )
            RETURN
         END IF
*        Note that if RANGE .NE. 'V', DLARRE computes bounds on the desired
*        part of the spectrum. All desired eigenvalues are contained in
*        (WL,WU]


         IF( WANTZ ) THEN
*
*           Compute the desired eigenvectors corresponding to the computed
*           eigenvalues
*
            CALL DLARRV( N, WL, WU, D, E, PIVMIN, IWORK( IINSPL ), M, 1, M, MINRGP, RTOL1, RTOL2, W, WORK( INDERR ), WORK( INDGP ), IWORK( IINDBL ), IWORK( IINDW ), WORK( INDGRS ), Z, LDZ, ISUPPZ, WORK( INDWRK ), IWORK( IINDWK ), IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = 20 + ABS( IINFO )
               RETURN
            END IF
         ELSE
*           DLARRE computes eigenvalues of the (shifted) root representation
*           DLARRV returns the eigenvalues of the unshifted matrix.
*           However, if the eigenvectors are not desired by the user, we need
*           to apply the corresponding shifts from DLARRE to obtain the
*           eigenvalues of the original matrix.
            DO 20 J = 1, M
               ITMP = IWORK( IINDBL+J-1 )
               W( J ) = W( J ) + E( IWORK( IINSPL+ITMP-1 ) )
 20         CONTINUE
         END IF
*

         IF ( TRYRAC ) THEN
*           Refine computed eigenvalues so that they are relatively accurate
*           with respect to the original matrix T.
            IBEGIN = 1
            WBEGIN = 1
            DO 39  JBLK = 1, IWORK( IINDBL+M-1 )
               IEND = IWORK( IINSPL+JBLK-1 )
               IN = IEND - IBEGIN + 1
               WEND = WBEGIN - 1
*              check if any eigenvalues have to be refined in this block
 36            CONTINUE
               IF( WEND.LT.M ) THEN
                  IF( IWORK( IINDBL+WEND ).EQ.JBLK ) THEN
                     WEND = WEND + 1
                     GO TO 36
                  END IF
               END IF
               IF( WEND.LT.WBEGIN ) THEN
                  IBEGIN = IEND + 1
                  GO TO 39
               END IF

               OFFSET = IWORK(IINDW+WBEGIN-1)-1
               IFIRST = IWORK(IINDW+WBEGIN-1)
               ILAST = IWORK(IINDW+WEND-1)
               RTOL2 = FOUR * EPS
               CALL DLARRJ( IN, WORK(INDD+IBEGIN-1), WORK(INDE2+IBEGIN-1), IFIRST, ILAST, RTOL2, OFFSET, W(WBEGIN), WORK( INDERR+WBEGIN-1 ), WORK( INDWRK ), IWORK( IINDWK ), PIVMIN, TNRM, IINFO )
               IBEGIN = IEND + 1
               WBEGIN = WEND + 1
 39         CONTINUE
         ENDIF
*
*        If matrix was scaled, then rescale eigenvalues appropriately.
*
         IF( SCALE.NE.ONE ) THEN
            CALL DSCAL( M, ONE / SCALE, W, 1 )
         END IF

      END IF

*
*     If eigenvalues are not in increasing order, then sort them,
*     possibly along with eigenvectors.
*
      IF( NSPLIT.GT.1 .OR. N.EQ.2 ) THEN
         IF( .NOT. WANTZ ) THEN
            CALL DLASRT( 'I', M, W, IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = 3
               RETURN
            END IF
         ELSE
            DO 60 J = 1, M - 1
               I = 0
               TMP = W( J )
               DO 50 JJ = J + 1, M
                  IF( W( JJ ).LT.TMP ) THEN
                     I = JJ
                     TMP = W( JJ )
                  END IF
 50            CONTINUE
               IF( I.NE.0 ) THEN
                  W( I ) = W( J )
                  W( J ) = TMP
                  IF( WANTZ ) THEN
                     CALL DSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
                     ITMP = ISUPPZ( 2*I-1 )
                     ISUPPZ( 2*I-1 ) = ISUPPZ( 2*J-1 )
                     ISUPPZ( 2*J-1 ) = ITMP
                     ITMP = ISUPPZ( 2*I )
                     ISUPPZ( 2*I ) = ISUPPZ( 2*J )
                     ISUPPZ( 2*J ) = ITMP
                  END IF
               END IF
 60         CONTINUE
         END IF
      ENDIF
*
*
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
      RETURN
*
*     End of DSTEMR
*
      END
