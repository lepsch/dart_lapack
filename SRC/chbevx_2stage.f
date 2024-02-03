      SUBROUTINE CHBEVX_2STAGE( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N, LWORK;
      REAL               ABSTOL, VL, VU
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      REAL               RWORK( * ), W( * )
      COMPLEX            AB( LDAB, * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E0, 0.0E0 ), CONE = ( 1.0E0, 0.0E0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LOWER, TEST, VALEIG, WANTZ, LQUERY;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, INDISP, INDIWK, INDRWK, INDWRK, ISCALE, ITMP1, LLWORK, LWMIN, LHTRD, LWTRD, IB, INDHOUS, J, JJ, NSPLIT;
      REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU
      COMPLEX            CTMP1
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      REAL               SLAMCH, CLANHB, SROUNDUP_LWORK
      // EXTERNAL LSAME, SLAMCH, CLANHB, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SSCAL, SSTEBZ, SSTERF, XERBLA, CCOPY, CGEMV, CLACPY, CLASCL, CSTEIN, CSTEQR, CSWAP, CHETRD_HB2ST
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )

      INFO = 0
      if ( .NOT.( LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) {
         INFO = -2
      } else if ( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( KD.LT.0 ) {
         INFO = -5
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -7
      } else if ( WANTZ .AND. LDQ.LT.MAX( 1, N ) ) {
         INFO = -9
      } else {
         if ( VALEIG ) {
            IF( N.GT.0 .AND. VU.LE.VL ) INFO = -11
         } else if ( INDEIG ) {
            if ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) {
               INFO = -12
            } else if ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) {
               INFO = -13
            }
         }
      }
      if ( INFO.EQ.0 ) {
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) INFO = -18
      }

      if ( INFO.EQ.0 ) {
         if ( N.LE.1 ) {
            LWMIN = 1
            WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         } else {
            IB    = ILAENV2STAGE( 2, 'CHETRD_HB2ST', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'CHETRD_HB2ST', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'CHETRD_HB2ST', JOBZ, N, KD, IB, -1 )
            LWMIN = LHTRD + LWTRD
            WORK( 1 )  = SROUNDUP_LWORK(LWMIN)
         ENDIF

         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) INFO = -20
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CHBEVX_2STAGE', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      M = 0
      IF( N.EQ.0 ) RETURN

      if ( N.EQ.1 ) {
         M = 1
         if ( LOWER ) {
            CTMP1 = AB( 1, 1 )
         } else {
            CTMP1 = AB( KD+1, 1 )
         }
         TMP1 = REAL( CTMP1 )
         if ( VALEIG ) {
            IF( .NOT.( VL.LT.TMP1 .AND. VU.GE.TMP1 ) ) M = 0
         }
         if ( M.EQ.1 ) {
            W( 1 ) = REAL( CTMP1 )
            IF( WANTZ ) Z( 1, 1 ) = CONE
         }
         RETURN
      }

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
      if ( VALEIG ) {
         VLL = VL
         VUU = VU
      } else {
         VLL = ZERO
         VUU = ZERO
      }
      ANRM = CLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK )
      if ( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE.EQ.1 ) {
         if ( LOWER ) {
            CALL CLASCL( 'B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         } else {
            CALL CLASCL( 'Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         }
         IF( ABSTOL.GT.0 ) ABSTLL = ABSTOL*SIGMA
         if ( VALEIG ) {
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         }
      }

      // Call CHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form.

      INDD = 1
      INDE = INDD + N
      INDRWK = INDE + N

      INDHOUS = 1
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1

      CALL CHETRD_HB2ST( 'N', JOBZ, UPLO, N, KD, AB, LDAB, RWORK( INDD ), RWORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO )

      // If all eigenvalues are desired and ABSTOL is less than or equal
      // to zero, then call SSTERF or CSTEQR.  If this fails for some
      // eigenvalue, then try SSTEBZ.

      TEST = .FALSE.
      if (INDEIG) {
         if (IL.EQ.1 .AND. IU.EQ.N) {
            TEST = .TRUE.
         }
      }
      if ((ALLEIG .OR. TEST) .AND. (ABSTOL.LE.ZERO)) {
         CALL SCOPY( N, RWORK( INDD ), 1, W, 1 )
         INDEE = INDRWK + 2*N
         if ( .NOT.WANTZ ) {
            CALL SCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
            CALL SSTERF( N, W, RWORK( INDEE ), INFO )
         } else {
            CALL CLACPY( 'A', N, N, Q, LDQ, Z, LDZ )
            CALL SCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
            CALL CSTEQR( JOBZ, N, W, RWORK( INDEE ), Z, LDZ, RWORK( INDRWK ), INFO )
            if ( INFO.EQ.0 ) {
               DO 10 I = 1, N
                  IFAIL( I ) = 0
   10          CONTINUE
            }
         }
         if ( INFO.EQ.0 ) {
            M = N
            GO TO 30
         }
         INFO = 0
      }

      // Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN.

      if ( WANTZ ) {
         ORDER = 'B'
      } else {
         ORDER = 'E'
      }
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWK = INDISP + N
      CALL SSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, RWORK( INDD ), RWORK( INDE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ), IWORK( INDIWK ), INFO )

      if ( WANTZ ) {
         CALL CSTEIN( N, RWORK( INDD ), RWORK( INDE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, RWORK( INDRWK ), IWORK( INDIWK ), IFAIL, INFO )

         // Apply unitary matrix used in reduction to tridiagonal
         // form to eigenvectors returned by CSTEIN.

         DO 20 J = 1, M
            CALL CCOPY( N, Z( 1, J ), 1, WORK( 1 ), 1 )
            CALL CGEMV( 'N', N, N, CONE, Q, LDQ, WORK, 1, CZERO, Z( 1, J ), 1 )
   20    CONTINUE
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

   30 CONTINUE
      if ( ISCALE.EQ.1 ) {
         if ( INFO.EQ.0 ) {
            IMAX = M
         } else {
            IMAX = INFO - 1
         }
         CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
      }

      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.

      if ( WANTZ ) {
         DO 50 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 40 JJ = J + 1, M
               if ( W( JJ ).LT.TMP1 ) {
                  I = JJ
                  TMP1 = W( JJ )
               }
   40       CONTINUE

            if ( I.NE.0 ) {
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL CSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
               if ( INFO.NE.0 ) {
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               }
            }
   50    CONTINUE
      }

      // Set WORK(1) to optimal workspace size.

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)

      RETURN

      // End of CHBEVX_2STAGE

      }
