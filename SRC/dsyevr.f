      void dsyevr(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N;
      double             ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * );
      double             A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LOWER, LQUERY, VALEIG, WANTZ, TRYRAC;
      String             ORDER;
      int                I, IEEEOK, IINFO, IMAX, INDD, INDDD, INDE, INDEE, INDIBL, INDIFL, INDISP, INDIWO, INDTAU, INDWK, INDWKN, ISCALE, J, JJ, LIWMIN, LLWORK, LLWRKN, LWKOPT, LWMIN, NB, NSPLIT;
      double             ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- int                ILAENV;
      //- double             DLAMCH, DLANSY;
      // EXTERNAL LSAME, ILAENV, DLAMCH, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DORMTR, DSCAL, DSTEBZ, DSTEMR, DSTEIN, DSTERF, DSWAP, DSYTRD, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      IEEEOK = ILAENV( 10, 'DSYEVR', 'N', 1, 2, 3, 4 );

      LOWER = LSAME( UPLO, 'L' );
      WANTZ = LSAME( JOBZ, 'V' );
      ALLEIG = LSAME( RANGE, 'A' );
      VALEIG = LSAME( RANGE, 'V' );
      INDEIG = LSAME( RANGE, 'I' );

      LQUERY = ( ( LWORK == -1 ) || ( LIWORK == -1 ) );

      if ( N <= 1 ) {
         LWMIN  = 1;
         LIWMIN = 1;
      } else {
         LWMIN  = 26*N;
         LIWMIN = 10*N;
      }

      INFO = 0;
      if ( !( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -2;
      } else if ( !( LOWER || LSAME( UPLO, 'U' ) ) ) {
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
         } else if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -18;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -20;
         }
      }

      if ( INFO == 0 ) {
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 );
         NB = max( NB, ILAENV( 1, 'DORMTR', UPLO, N, -1, -1, -1 ) );
         LWKOPT = max( ( NB+1 )*N, LWMIN );
         WORK( 1 ) = LWKOPT;
         IWORK( 1 ) = LIWMIN;
      }

      if ( INFO != 0 ) {
         xerbla('DSYEVR', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      M = 0;
      if ( N == 0 ) {
         WORK( 1 ) = 1;
         return;
      }

      if ( N == 1 ) {
         WORK( 1 ) = 1;
         if ( ALLEIG || INDEIG ) {
            M = 1;
            W( 1 ) = A( 1, 1 );
         } else {
            if ( VL < A( 1, 1 ) && VU >= A( 1, 1 ) ) {
               M = 1;
               W( 1 ) = A( 1, 1 );
            }
         }
         if ( WANTZ ) {
            Z( 1, 1 ) = ONE;
            ISUPPZ( 1 ) = 1;
            ISUPPZ( 2 ) = 1;
         }
         return;
      }

      // Get machine constants.

      SAFMIN = DLAMCH( 'Safe minimum' );
      EPS = DLAMCH( 'Precision' );
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
      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK );
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
               dscal(N-J+1, SIGMA, A( J, J ), 1 );
            } // 10
         } else {
            for (J = 1; J <= N; J++) { // 20
               dscal(J, SIGMA, A( 1, J ), 1 );
            } // 20
         }
         if (ABSTOL > 0) ABSTLL = ABSTOL*SIGMA;
         if ( VALEIG ) {
            VLL = VL*SIGMA;
            VUU = VU*SIGMA;
         }
      }

      // Initialize indices into workspaces.  Note: The IWORK indices are
      // used only if DSTERF or DSTEMR fail.

      // WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the
      // elementary reflectors used in DSYTRD.
      INDTAU = 1;
      // WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries.
      INDD = INDTAU + N;
      // WORK(INDE:INDE+N-1) stores the off-diagonal entries of the
      // tridiagonal matrix from DSYTRD.
      INDE = INDD + N;
      // WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over
      // -written by DSTEMR (the DSTERF path copies the diagonal to W).
      INDDD = INDE + N;
      // WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over
      // -written while computing the eigenvalues in DSTERF and DSTEMR.
      INDEE = INDDD + N;
      // INDWK is the starting offset of the left-over workspace, and
      // LLWORK is the remaining workspace size.
      INDWK = INDEE + N;
      LLWORK = LWORK - INDWK + 1;

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


      // Call DSYTRD to reduce symmetric matrix to tridiagonal form.

      dsytrd(UPLO, N, A, LDA, WORK( INDD ), WORK( INDE ), WORK( INDTAU ), WORK( INDWK ), LLWORK, IINFO );

      // If all eigenvalues are desired
      // then call DSTERF or DSTEMR and DORMTR.

      if ( ( ALLEIG || ( INDEIG && IL == 1 && IU == N ) ) && IEEEOK == 1 ) {
         if ( !WANTZ ) {
            dcopy(N, WORK( INDD ), 1, W, 1 );
            dcopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            dsterf(N, W, WORK( INDEE ), INFO );
         } else {
            dcopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            dcopy(N, WORK( INDD ), 1, WORK( INDDD ), 1 );

            if (ABSTOL <= TWO*N*EPS) {
               TRYRAC = true;
            } else {
               TRYRAC = false;
            }
            dstemr(JOBZ, 'A', N, WORK( INDDD ), WORK( INDEE ), VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK( INDWK ), LWORK, IWORK, LIWORK, INFO );



         // Apply orthogonal matrix used in reduction to tridiagonal
         // form to eigenvectors returned by DSTEMR.

            if ( WANTZ && INFO == 0 ) {
               INDWKN = INDE;
               LLWRKN = LWORK - INDWKN + 1;
               dormtr('L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO );
            }
         }


         if ( INFO == 0 ) {
            // Everything worked.  Skip DSTEBZ/DSTEIN.  IWORK(:) are
            // undefined.
            M = N;
            GO TO 30;
         }
         INFO = 0;
      }

      // Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN.
      // Also call DSTEBZ and DSTEIN if DSTEMR fails.

      if ( WANTZ ) {
         ORDER = 'B';
      } else {
         ORDER = 'E';
      }
       dstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, WORK( INDD ), WORK( INDE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWK ), IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         dstein(N, WORK( INDD ), WORK( INDE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, WORK( INDWK ), IWORK( INDIWO ), IWORK( INDIFL ), INFO );

         // Apply orthogonal matrix used in reduction to tridiagonal
         // form to eigenvectors returned by DSTEIN.

         INDWKN = INDE;
         LLWRKN = LWORK - INDWKN + 1;
         dormtr('L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

// Jump here if DSTEMR/DSTEIN succeeded.
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
      // eigenvectors.  Note: We do not sort the IFAIL portion of IWORK.
      // It may not be initialized (if DSTEMR/DSTEIN succeeded), and we do
      // not return this detailed information to the user.

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
               W( I ) = W( J );
               W( J ) = TMP1;
               dswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
            }
         } // 50
      }

      // Set WORK(1) to optimal workspace size.

      WORK( 1 ) = LWKOPT;
      IWORK( 1 ) = LIWMIN;

      return;
      }
