      SUBROUTINE SSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N;
      REAL               ABSTOL, VL, VU
      // ..
      // .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * );
      REAL               A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LOWER, LQUERY, TEST, VALEIG, WANTZ, TRYRAC;
      String             ORDER;
      int                I, IEEEOK, IINFO, IMAX, INDD, INDDD, INDE, INDEE, INDIBL, INDIFL, INDISP, INDIWO, INDTAU, INDWK, INDWKN, ISCALE, J, JJ, LIWMIN, LLWORK, LLWRKN, LWKOPT, LWMIN, NB, NSPLIT;
      REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANSY, SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SLAMCH, SLANSY, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SORMTR, SSCAL, SSTEBZ, SSTEMR, SSTEIN, SSTERF, SSWAP, SSYTRD, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      IEEEOK = ILAENV( 10, 'SSYEVR', 'N', 1, 2, 3, 4 )

      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )

      LQUERY = ( ( LWORK.EQ.-1 ) .OR. ( LIWORK.EQ.-1 ) )

      if ( N.LE.1 ) {
         LWMIN  = 1
         LIWMIN = 1
      } else {
         LWMIN  = 26*N
         LIWMIN = 10*N
      }

      INFO = 0
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) {
         INFO = -2
      } else if ( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      } else {
         if ( VALEIG ) {
            IF( N.GT.0 .AND. VU.LE.VL ) INFO = -8
         } else if ( INDEIG ) {
            if ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) {
               INFO = -9
            } else if ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) {
               INFO = -10
            }
         }
      }
      if ( INFO.EQ.0 ) {
         if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
            INFO = -15
         }
      }

      if ( INFO.EQ.0 ) {
         NB = ILAENV( 1, 'SSYTRD', UPLO, N, -1, -1, -1 )
         NB = MAX( NB, ILAENV( 1, 'SORMTR', UPLO, N, -1, -1, -1 ) )
         LWKOPT = MAX( ( NB+1 )*N, LWMIN )
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -18
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
            INFO = -20
         }
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SSYEVR', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      M = 0
      if ( N.EQ.0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      if ( N.EQ.1 ) {
         WORK( 1 ) = 26
         if ( ALLEIG .OR. INDEIG ) {
            M = 1
            W( 1 ) = A( 1, 1 )
         } else {
            if ( VL.LT.A( 1, 1 ) .AND. VU.GE.A( 1, 1 ) ) {
               M = 1
               W( 1 ) = A( 1, 1 )
            }
         }
         if ( WANTZ ) {
            Z( 1, 1 ) = ONE
            ISUPPZ( 1 ) = 1
            ISUPPZ( 2 ) = 1
         }
         RETURN
      }

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )

      // Scale matrix to allowable range, if necessary.

      ISCALE = 0
      ABSTLL = ABSTOL
      if (VALEIG) {
         VLL = VL
         VUU = VU
      }
      ANRM = SLANSY( 'M', UPLO, N, A, LDA, WORK )
      if ( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE.EQ.1 ) {
         if ( LOWER ) {
            DO 10 J = 1, N
               CALL SSCAL( N-J+1, SIGMA, A( J, J ), 1 )
   10       CONTINUE
         } else {
            DO 20 J = 1, N
               CALL SSCAL( J, SIGMA, A( 1, J ), 1 )
   20       CONTINUE
         }
         IF( ABSTOL.GT.0 ) ABSTLL = ABSTOL*SIGMA
         if ( VALEIG ) {
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         }
      }

      // Initialize indices into workspaces.  Note: The IWORK indices are
      // used only if SSTERF or SSTEMR fail.

      // WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the
      // elementary reflectors used in SSYTRD.
      INDTAU = 1
      // WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries.
      INDD = INDTAU + N
      // WORK(INDE:INDE+N-1) stores the off-diagonal entries of the
      // tridiagonal matrix from SSYTRD.
      INDE = INDD + N
      // WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over
      // -written by SSTEMR (the SSTERF path copies the diagonal to W).
      INDDD = INDE + N
      // WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over
      // -written while computing the eigenvalues in SSTERF and SSTEMR.
      INDEE = INDDD + N
      // INDWK is the starting offset of the left-over workspace, and
      // LLWORK is the remaining workspace size.
      INDWK = INDEE + N
      LLWORK = LWORK - INDWK + 1

      // IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in SSTEBZ and
      // stores the block indices of each of the M<=N eigenvalues.
      INDIBL = 1
      // IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in SSTEBZ and
      // stores the starting and finishing indices of each block.
      INDISP = INDIBL + N
      // IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
      // that corresponding to eigenvectors that fail to converge in
      // SSTEIN.  This information is discarded; if any fail, the driver
      // returns INFO > 0.
      INDIFL = INDISP + N
      // INDIWO is the offset of the remaining integer workspace.
      INDIWO = INDIFL + N


      // Call SSYTRD to reduce symmetric matrix to tridiagonal form.

      CALL SSYTRD( UPLO, N, A, LDA, WORK( INDD ), WORK( INDE ), WORK( INDTAU ), WORK( INDWK ), LLWORK, IINFO )

      // If all eigenvalues are desired
      // then call SSTERF or SSTEMR and SORMTR.

      TEST = .FALSE.
      if ( INDEIG ) {
         if ( IL.EQ.1 .AND. IU.EQ.N ) {
            TEST = .TRUE.
         }
      }
      if ( ( ALLEIG.OR.TEST ) .AND. ( IEEEOK.EQ.1 ) ) {
         if ( .NOT.WANTZ ) {
            CALL SCOPY( N, WORK( INDD ), 1, W, 1 )
            CALL SCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL SSTERF( N, W, WORK( INDEE ), INFO )
         } else {
            CALL SCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL SCOPY( N, WORK( INDD ), 1, WORK( INDDD ), 1 )

            if (ABSTOL .LE. TWO*N*EPS) {
               TRYRAC = .TRUE.
            } else {
               TRYRAC = .FALSE.
            }
            CALL SSTEMR( JOBZ, 'A', N, WORK( INDDD ), WORK( INDEE ), VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK( INDWK ), LWORK, IWORK, LIWORK, INFO )



         // Apply orthogonal matrix used in reduction to tridiagonal
         // form to eigenvectors returned by SSTEMR.

            if ( WANTZ .AND. INFO.EQ.0 ) {
               INDWKN = INDE
               LLWRKN = LWORK - INDWKN + 1
               CALL SORMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO )
            }
         }


         if ( INFO.EQ.0 ) {
            // Everything worked.  Skip SSTEBZ/SSTEIN.  IWORK(:) are
            // undefined.
            M = N
            GO TO 30
         }
         INFO = 0
      }

      // Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN.
      // Also call SSTEBZ and SSTEIN if SSTEMR fails.

      if ( WANTZ ) {
         ORDER = 'B'
      } else {
         ORDER = 'E'
      }
       CALL SSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, WORK( INDD ), WORK( INDE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWK ), IWORK( INDIWO ), INFO )

      if ( WANTZ ) {
         CALL SSTEIN( N, WORK( INDD ), WORK( INDE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, WORK( INDWK ), IWORK( INDIWO ), IWORK( INDIFL ), INFO )

         // Apply orthogonal matrix used in reduction to tridiagonal
         // form to eigenvectors returned by SSTEIN.

         INDWKN = INDE
         LLWRKN = LWORK - INDWKN + 1
         CALL SORMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO )
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

*  Jump here if SSTEMR/SSTEIN succeeded.
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
      // eigenvectors.  Note: We do not sort the IFAIL portion of IWORK.
      // It may not be initialized (if SSTEMR/SSTEIN succeeded), and we do
      // not return this detailed information to the user.

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
               W( I ) = W( J )
               W( J ) = TMP1
               CALL SSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
            }
   50    CONTINUE
      }

      // Set WORK(1) to optimal workspace size.

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of SSYEVR

      }
