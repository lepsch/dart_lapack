      SUBROUTINE CHEEVR_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDA, LDZ, LIWORK, LRWORK, LWORK, M, N;
      REAL               ABSTOL, VL, VU
      // ..
      // .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * );
      REAL               RWORK( * ), W( * )
      COMPLEX            A( LDA, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LOWER, LQUERY, TEST, VALEIG, WANTZ, TRYRAC;
      String             ORDER;
      int                I, IEEEOK, IINFO, IMAX, INDIBL, INDIFL, INDISP, INDIWO, INDRD, INDRDD, INDRE, INDREE, INDRWK, INDTAU, INDWK, INDWKN, ISCALE, ITMP1, J, JJ, LIWMIN, LLWORK, LLRWORK, LLWRKN, LRWMIN, LWMIN, NSPLIT, LHTRD, LWTRD, KD, IB, INDHOUS;
      REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV, ILAENV2STAGE;
      REAL               SLAMCH, CLANSY, SROUNDUP_LWORK
      // EXTERNAL LSAME, SLAMCH, CLANSY, ILAENV, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SSCAL, SSTEBZ, SSTERF, XERBLA, CSSCAL, CHETRD_2STAGE, CSTEMR, CSTEIN, CSWAP, CUNMTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      IEEEOK = ILAENV( 10, 'CHEEVR', 'N', 1, 2, 3, 4 )

      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )

      LQUERY = ( ( LWORK == -1 ) || ( LRWORK == -1 ) || ( LIWORK == -1 ) )

      KD     = ILAENV2STAGE( 1, 'CHETRD_2STAGE', JOBZ, N, -1, -1, -1 )
      IB     = ILAENV2STAGE( 2, 'CHETRD_2STAGE', JOBZ, N, KD, -1, -1 )
      LHTRD  = ILAENV2STAGE( 3, 'CHETRD_2STAGE', JOBZ, N, KD, IB, -1 )
      LWTRD  = ILAENV2STAGE( 4, 'CHETRD_2STAGE', JOBZ, N, KD, IB, -1 )

      if ( N.LE.1 ) {
         LWMIN  = 1
         LRWMIN = 1
         LIWMIN = 1
      } else {
         LWMIN  = N + LHTRD + LWTRD
         LRWMIN = 24*N
         LIWMIN = 10*N
      }

      INFO = 0
      if ( .NOT.( LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -2
      } else if ( .NOT.( LOWER || LSAME( UPLO, 'U' ) ) ) {
         INFO = -3
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -6
      } else {
         if ( VALEIG ) {
            if (N > 0 && VU.LE.VL) INFO = -8;
         } else if ( INDEIG ) {
            if ( IL < 1 || IL > MAX( 1, N ) ) {
               INFO = -9
            } else if ( IU < MIN( N, IL ) || IU > N ) {
               INFO = -10
            }
         }
      }
      if ( INFO == 0 ) {
         if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
            INFO = -15
         }
      }

      if ( INFO == 0 ) {
         WORK( 1 )  = SROUNDUP_LWORK( LWMIN )
         RWORK( 1 ) = SROUNDUP_LWORK( LRWMIN )
         IWORK( 1 ) = LIWMIN

         if ( LWORK < LWMIN && .NOT.LQUERY ) {
            INFO = -18
         } else if ( LRWORK < LRWMIN && .NOT.LQUERY ) {
            INFO = -20
         } else if ( LIWORK < LIWMIN && .NOT.LQUERY ) {
            INFO = -22
         }
      }

      if ( INFO != 0 ) {
         xerbla('CHEEVR_2STAGE', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      M = 0
      if ( N == 0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      if ( N == 1 ) {
         WORK( 1 ) = 1
         if ( ALLEIG || INDEIG ) {
            M = 1
            W( 1 ) = REAL( A( 1, 1 ) )
         } else {
            if ( VL < REAL( A( 1, 1 ) ) && VU.GE.REAL( A( 1, 1 ) ) ) {
               M = 1
               W( 1 ) = REAL( A( 1, 1 ) )
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
      EPS    = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )

      // Scale matrix to allowable range, if necessary.

      ISCALE = 0
      ABSTLL = ABSTOL
      if (VALEIG) {
         VLL = VL
         VUU = VU
      }
      ANRM = CLANSY( 'M', UPLO, N, A, LDA, RWORK )
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM > RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE == 1 ) {
         if ( LOWER ) {
            for (J = 1; J <= N; J++) { // 10
               csscal(N-J+1, SIGMA, A( J, J ), 1 );
            } // 10
         } else {
            for (J = 1; J <= N; J++) { // 20
               csscal(J, SIGMA, A( 1, J ), 1 );
            } // 20
         }
         if (ABSTOL > 0) ABSTLL = ABSTOL*SIGMA;
         if ( VALEIG ) {
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         }
      }

      // Initialize indices into workspaces.  Note: The IWORK indices are
      // used only if SSTERF or CSTEMR fail.

      // WORK(INDTAU:INDTAU+N-1) stores the complex scalar factors of the
      // elementary reflectors used in CHETRD.
      INDTAU = 1
      // INDWK is the starting offset of the remaining complex workspace,
      // and LLWORK is the remaining complex workspace size.
      INDHOUS = INDTAU + N
      INDWK   = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWK + 1

      // RWORK(INDRD:INDRD+N-1) stores the real tridiagonal's diagonal
      // entries.
      INDRD = 1
      // RWORK(INDRE:INDRE+N-1) stores the off-diagonal entries of the
      // tridiagonal matrix from CHETRD.
      INDRE = INDRD + N
      // RWORK(INDRDD:INDRDD+N-1) is a copy of the diagonal entries over
      // -written by CSTEMR (the SSTERF path copies the diagonal to W).
      INDRDD = INDRE + N
      // RWORK(INDREE:INDREE+N-1) is a copy of the off-diagonal entries over
      // -written while computing the eigenvalues in SSTERF and CSTEMR.
      INDREE = INDRDD + N
      // INDRWK is the starting offset of the left-over real workspace, and
      // LLRWORK is the remaining workspace size.
      INDRWK = INDREE + N
      LLRWORK = LRWORK - INDRWK + 1

      // IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in SSTEBZ and
      // stores the block indices of each of the M<=N eigenvalues.
      INDIBL = 1
      // IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in SSTEBZ and
      // stores the starting and finishing indices of each block.
      INDISP = INDIBL + N
      // IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
      // that corresponding to eigenvectors that fail to converge in
      // CSTEIN.  This information is discarded; if any fail, the driver
      // returns INFO > 0.
      INDIFL = INDISP + N
      // INDIWO is the offset of the remaining integer workspace.
      INDIWO = INDIFL + N


      // Call CHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.

      chetrd_2stage(JOBZ, UPLO, N, A, LDA, RWORK( INDRD ), RWORK( INDRE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWK ), LLWORK, IINFO );

      // If all eigenvalues are desired
      // then call SSTERF or CSTEMR and CUNMTR.

      TEST = false;
      if ( INDEIG ) {
         if ( IL == 1 && IU == N ) {
            TEST = true;
         }
      }
      if ( ( ALLEIG || TEST ) && ( IEEEOK == 1 ) ) {
         if ( .NOT.WANTZ ) {
            scopy(N, RWORK( INDRD ), 1, W, 1 );
            scopy(N-1, RWORK( INDRE ), 1, RWORK( INDREE ), 1 );
            ssterf(N, W, RWORK( INDREE ), INFO );
         } else {
            scopy(N-1, RWORK( INDRE ), 1, RWORK( INDREE ), 1 );
            scopy(N, RWORK( INDRD ), 1, RWORK( INDRDD ), 1 );

            if ( ABSTOL .LE. TWO*N*EPS ) {
               TRYRAC = true;
            } else {
               TRYRAC = false;
            }
            cstemr(JOBZ, 'A', N, RWORK( INDRDD ), RWORK( INDREE ), VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, RWORK( INDRWK ), LLRWORK, IWORK, LIWORK, INFO );

            // Apply unitary matrix used in reduction to tridiagonal
            // form to eigenvectors returned by CSTEMR.

            if ( WANTZ && INFO == 0 ) {
               INDWKN = INDWK
               LLWRKN = LWORK - INDWKN + 1
               cunmtr('L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO );
            }
         }


         if ( INFO == 0 ) {
            M = N
            GO TO 30
         }
         INFO = 0
      }

      // Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN.
      // Also call SSTEBZ and CSTEIN if CSTEMR fails.

      if ( WANTZ ) {
         ORDER = 'B'
      } else {
         ORDER = 'E'
      }
       sstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, RWORK( INDRD ), RWORK( INDRE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ), IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         cstein(N, RWORK( INDRD ), RWORK( INDRE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, RWORK( INDRWK ), IWORK( INDIWO ), IWORK( INDIFL ), INFO );

         // Apply unitary matrix used in reduction to tridiagonal
         // form to eigenvectors returned by CSTEIN.

         INDWKN = INDWK
         LLWRKN = LWORK - INDWKN + 1
         cunmtr('L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      } // 30
      if ( ISCALE == 1 ) {
         if ( INFO == 0 ) {
            IMAX = M
         } else {
            IMAX = INFO - 1
         }
         sscal(IMAX, ONE / SIGMA, W, 1 );
      }

      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.

      if ( WANTZ ) {
         for (J = 1; J <= M - 1; J++) { // 50
            I = 0
            TMP1 = W( J )
            for (JJ = J + 1; JJ <= M; JJ++) { // 40
               if ( W( JJ ) < TMP1 ) {
                  I = JJ
                  TMP1 = W( JJ )
               }
            } // 40

            if ( I != 0 ) {
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               cswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
            }
         } // 50
      }

      // Set WORK(1) to optimal workspace size.

      WORK( 1 )  = SROUNDUP_LWORK( LWMIN )
      RWORK( 1 ) = SROUNDUP_LWORK( LRWMIN )
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of CHEEVR_2STAGE

      }
