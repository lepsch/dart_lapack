      void zheevr(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDA, LDZ, LIWORK, LRWORK, LWORK, M, N;
      double             ABSTOL, VL, VU;
      int                ISUPPZ( * ), IWORK( * );
      double             RWORK( * ), W( * );
      Complex         A( LDA, * ), WORK( * ), Z( LDZ, * );
      // ..

      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      bool               ALLEIG, INDEIG, LOWER, LQUERY, TEST, VALEIG, WANTZ, TRYRAC;
      String             ORDER;
      int                I, IEEEOK, IINFO, IMAX, INDIBL, INDIFL, INDISP, INDIWO, INDRD, INDRDD, INDRE, INDREE, INDRWK, INDTAU, INDWK, INDWKN, ISCALE, ITMP1, J, JJ, LIWMIN, LLWORK, LLRWORK, LLWRKN, LRWMIN, LWKOPT, LWMIN, NB, NSPLIT;
      double             ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- double             DLAMCH, ZLANSY;
      // EXTERNAL lsame, ILAENV, DLAMCH, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSCAL, DSTEBZ, DSTERF, XERBLA, ZDSCAL, ZHETRD, ZSTEMR, ZSTEIN, ZSWAP, ZUNMTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, SQRT

      // Test the input parameters.

      IEEEOK = ilaenv( 10, 'ZHEEVR', 'N', 1, 2, 3, 4 );

      LOWER = lsame( UPLO, 'L' );
      WANTZ = lsame( JOBZ, 'V' );
      ALLEIG = lsame( RANGE, 'A' );
      VALEIG = lsame( RANGE, 'V' );
      INDEIG = lsame( RANGE, 'I' );

      LQUERY = ( ( LWORK == -1 ) || ( LRWORK == -1 ) || ( LIWORK == -1 ) );

      if ( N <= 1 ) {
         LWMIN  = 1;
         LRWMIN = 1;
         LIWMIN = 1;
      } else {
         LWMIN  = 2*N;
         LRWMIN = 24*N;
         LIWMIN = 10*N;
      }

      INFO = 0;
      if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -2;
      } else if ( !( LOWER || lsame( UPLO, 'U' ) ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -6;
      } else {
         if ( VALEIG ) {
            if (N > 0 && VU <= VL) INFO = -8;
         } else if ( INDEIG ) {
            if ( IL < 1 || IL > max( 1, N ) ) {
               INFO = -9;
            } else if ( IU < min( N, IL ) || IU > N ) {
               INFO = -10;
            }
         }
      }
      if ( INFO == 0 ) {
         if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
            INFO = -15;
         }
      }

      if ( INFO == 0 ) {
         NB = ilaenv( 1, 'ZHETRD', UPLO, N, -1, -1, -1 );
         NB = max( NB, ilaenv( 1, 'ZUNMTR', UPLO, N, -1, -1, -1 ) );
         LWKOPT = max( ( NB+1 )*N, LWMIN );
         WORK[1] = LWKOPT;
         RWORK[1] = LRWMIN;
         IWORK[1] = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -18;
         } else if ( LRWORK < LRWMIN && !LQUERY ) {
            INFO = -20;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -22;
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZHEEVR', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      M = 0;
      if ( N == 0 ) {
         WORK[1] = 1;
         return;
      }

      if ( N == 1 ) {
         WORK[1] = 1;
         if ( ALLEIG || INDEIG ) {
            M = 1;
            W[1] = (A( 1, 1 )).toDouble();
         } else {
            if ( VL < (A( 1, 1 )).toDouble() && VU >= (A( 1, 1 )).toDouble() ) {
               M = 1;
               W[1] = (A( 1, 1 )).toDouble();
            }
         }
         if ( WANTZ ) {
            Z[1][1] = ONE;
            ISUPPZ[1] = 1;
            ISUPPZ[2] = 1;
         }
         return;
      }

      // Get machine constants.

      SAFMIN = dlamch( 'Safe minimum' );
      EPS = dlamch( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN = sqrt( SMLNUM );
      RMAX = min( sqrt( BIGNUM ), ONE / sqrt( sqrt( SAFMIN ) ) );

      // Scale matrix to allowable range, if necessary.

      ISCALE = 0;
      ABSTLL = ABSTOL;
      if (VALEIG) {
         VLL = VL;
         VUU = VU;
      }
      ANRM = ZLANSY( 'M', UPLO, N, A, LDA, RWORK );
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if ( ISCALE == 1 ) {
         if ( LOWER ) {
            for (J = 1; J <= N; J++) { // 10
               zdscal(N-J+1, SIGMA, A( J, J ), 1 );
            } // 10
         } else {
            for (J = 1; J <= N; J++) { // 20
               zdscal(J, SIGMA, A( 1, J ), 1 );
            } // 20
         }
         if (ABSTOL > 0) ABSTLL = ABSTOL*SIGMA;
         if ( VALEIG ) {
            VLL = VL*SIGMA;
            VUU = VU*SIGMA;
         }
      }

      // Initialize indices into workspaces.  Note: The IWORK indices are
      // used only if DSTERF or ZSTEMR fail.

      // WORK(INDTAU:INDTAU+N-1) stores the complex scalar factors of the
      // elementary reflectors used in ZHETRD.
      INDTAU = 1;
      // INDWK is the starting offset of the remaining complex workspace,
      // and LLWORK is the remaining complex workspace size.
      INDWK = INDTAU + N;
      LLWORK = LWORK - INDWK + 1;

      // RWORK(INDRD:INDRD+N-1) stores the real tridiagonal's diagonal
      // entries.
      INDRD = 1;
      // RWORK(INDRE:INDRE+N-1) stores the off-diagonal entries of the
      // tridiagonal matrix from ZHETRD.
      INDRE = INDRD + N;
      // RWORK(INDRDD:INDRDD+N-1) is a copy of the diagonal entries over
      // -written by ZSTEMR (the DSTERF path copies the diagonal to W).
      INDRDD = INDRE + N;
      // RWORK(INDREE:INDREE+N-1) is a copy of the off-diagonal entries over
      // -written while computing the eigenvalues in DSTERF and ZSTEMR.
      INDREE = INDRDD + N;
      // INDRWK is the starting offset of the left-over real workspace, and
      // LLRWORK is the remaining workspace size.
      INDRWK = INDREE + N;
      LLRWORK = LRWORK - INDRWK + 1;

      // IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and
      // stores the block indices of each of the M<=N eigenvalues.
      INDIBL = 1;
      // IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and
      // stores the starting and finishing indices of each block.
      INDISP = INDIBL + N;
      // IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
      // that corresponding to eigenvectors that fail to converge in
      // DSTEIN.  This information is discarded; if any fail, the driver
      // returns INFO > 0.
      INDIFL = INDISP + N;
      // INDIWO is the offset of the remaining integer workspace.
      INDIWO = INDIFL + N;


      // Call ZHETRD to reduce Hermitian matrix to tridiagonal form.

      zhetrd(UPLO, N, A, LDA, RWORK( INDRD ), RWORK( INDRE ), WORK( INDTAU ), WORK( INDWK ), LLWORK, IINFO );

      // If all eigenvalues are desired
      // then call DSTERF or ZSTEMR and ZUNMTR.

      TEST = false;
      if ( INDEIG ) {
         if ( IL == 1 && IU == N ) {
            TEST = true;
         }
      }
      if ( ( ALLEIG || TEST ) && ( IEEEOK == 1 ) ) {
         if ( !WANTZ ) {
            dcopy(N, RWORK( INDRD ), 1, W, 1 );
            dcopy(N-1, RWORK( INDRE ), 1, RWORK( INDREE ), 1 );
            dsterf(N, W, RWORK( INDREE ), INFO );
         } else {
            dcopy(N-1, RWORK( INDRE ), 1, RWORK( INDREE ), 1 );
            dcopy(N, RWORK( INDRD ), 1, RWORK( INDRDD ), 1 );

            if (ABSTOL <= TWO*N*EPS) {
               TRYRAC = true;
            } else {
               TRYRAC = false;
            }
            zstemr(JOBZ, 'A', N, RWORK( INDRDD ), RWORK( INDREE ), VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, RWORK( INDRWK ), LLRWORK, IWORK, LIWORK, INFO );

            // Apply unitary matrix used in reduction to tridiagonal
            // form to eigenvectors returned by ZSTEMR.

            if ( WANTZ && INFO == 0 ) {
               INDWKN = INDWK;
               LLWRKN = LWORK - INDWKN + 1;
               zunmtr('L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO );
            }
         }


         if ( INFO == 0 ) {
            M = N;
            GO TO 30;
         }
         INFO = 0;
      }

      // Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.
      // Also call DSTEBZ and ZSTEIN if ZSTEMR fails.

      if ( WANTZ ) {
         ORDER = 'B';
      } else {
         ORDER = 'E';
      }
       dstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, RWORK( INDRD ), RWORK( INDRE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ), IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         zstein(N, RWORK( INDRD ), RWORK( INDRE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, RWORK( INDRWK ), IWORK( INDIWO ), IWORK( INDIFL ), INFO );

         // Apply unitary matrix used in reduction to tridiagonal
         // form to eigenvectors returned by ZSTEIN.

         INDWKN = INDWK;
         LLWRKN = LWORK - INDWKN + 1;
         zunmtr('L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      } // 30
      if ( ISCALE == 1 ) {
         if ( INFO == 0 ) {
            IMAX = M;
         } else {
            IMAX = INFO - 1;
         }
         dscal(IMAX, ONE / SIGMA, W, 1 );
      }

      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.

      if ( WANTZ ) {
         for (J = 1; J <= M - 1; J++) { // 50
            I = 0;
            TMP1 = W( J );
            for (JJ = J + 1; JJ <= M; JJ++) { // 40
               if ( W( JJ ) < TMP1 ) {
                  I = JJ;
                  TMP1 = W( JJ );
               }
            } // 40

            if ( I != 0 ) {
               ITMP1 = IWORK( INDIBL+I-1 );
               W[I] = W( J );
               IWORK[INDIBL+I-1] = IWORK( INDIBL+J-1 );
               W[J] = TMP1;
               IWORK[INDIBL+J-1] = ITMP1;
               zswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
            }
         } // 50
      }

      // Set WORK(1) to optimal workspace size.

      WORK[1] = LWKOPT;
      RWORK[1] = LRWMIN;
      IWORK[1] = LIWMIN;

      }
