      void zheevx_2stage(final int JOBZ, final int RANGE, final int UPLO, final int N, final Matrix<double> A_, final int LDA, final int VL, final int VU, final int IL, final int IU, final int ABSTOL, final int M, final int W, final Matrix<double> Z_, final int LDZ, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final Array<int> IWORK_, final int IFAIL, final Box<int> INFO,) {
  final A = A_.dim();
  final Z = Z_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDA, LDZ, LWORK, M, N;
      double             ABSTOL, VL, VU;
      int                IFAIL( * ), IWORK( * );
      double             RWORK( * ), W( * );
      Complex         A( LDA, * ), WORK( * ), Z( LDZ, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      bool               ALLEIG, INDEIG, LOWER, LQUERY, TEST, VALEIG, WANTZ;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, INDISP, INDIWK, INDRWK, INDTAU, INDWRK, ISCALE, ITMP1, J, JJ, LLWORK, NSPLIT, LWMIN, LHTRD, LWTRD, KD, IB, INDHOUS;
      double             ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV2STAGE;
      //- double             DLAMCH, ZLANHE;
      // EXTERNAL lsame, DLAMCH, ZLANHE, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSCAL, DSTEBZ, DSTERF, XERBLA, ZDSCAL, ZLACPY, ZSTEIN, ZSTEQR, ZSWAP, ZUNGTR, ZUNMTR, ZHETRD_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, SQRT

      // Test the input parameters.

      LOWER = lsame( UPLO, 'L' );
      WANTZ = lsame( JOBZ, 'V' );
      ALLEIG = lsame( RANGE, 'A' );
      VALEIG = lsame( RANGE, 'V' );
      INDEIG = lsame( RANGE, 'I' );
      LQUERY = ( LWORK == -1 );

      INFO = 0;
      if ( !( lsame( JOBZ, 'N' ) ) ) {
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
         if ( N <= 1 ) {
            LWMIN = 1;
            WORK[1] = LWMIN;
         } else {
            KD    = ILAENV2STAGE( 1, 'ZHETRD_2STAGE', JOBZ, N, -1, -1, -1 )             IB    = ILAENV2STAGE( 2, 'ZHETRD_2STAGE', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 );
            LWMIN = N + LHTRD + LWTRD;
            WORK[1] = LWMIN;
         }

         if (LWORK < LWMIN && !LQUERY) INFO = -17;
      }

      if ( INFO != 0 ) {
         xerbla('ZHEEVX_2STAGE', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      M = 0;
      if ( N == 0 ) {
         return;
      }

      if ( N == 1 ) {
         if ( ALLEIG || INDEIG ) {
            M = 1;
         W[1] = (A( 1, 1 )).toDouble();
         } else if ( VALEIG ) {
            if ( VL < (A( 1, 1 )).toDouble() && VU >= (A( 1, 1 )).toDouble() ) {
               M = 1;
               W[1] = (A( 1, 1 )).toDouble();
            }
         }
         if (WANTZ) Z( 1, 1 ) = CONE;
         return;
      }

      // Get machine constants.

      SAFMIN = dlamch( 'Safe minimum' );
      EPS    = dlamch( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN   = sqrt( SMLNUM );
      RMAX   = min( sqrt( BIGNUM ), ONE / sqrt( sqrt( SAFMIN ) ) );

      // Scale matrix to allowable range, if necessary.

      ISCALE = 0;
      ABSTLL = ABSTOL;
      if ( VALEIG ) {
         VLL = VL;
         VUU = VU;
      }
      ANRM = ZLANHE( 'M', UPLO, N, A, LDA, RWORK );
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

      // Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.

      INDD    = 1;
      INDE    = INDD + N;
      INDRWK  = INDE + N;
      INDTAU  = 1;
      INDHOUS = INDTAU + N;
      INDWRK  = INDHOUS + LHTRD;
      LLWORK  = LWORK - INDWRK + 1;

      zhetrd_2stage(JOBZ, UPLO, N, A, LDA, RWORK( INDD ), RWORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // If all eigenvalues are desired and ABSTOL is less than or equal to
      // zero, then call DSTERF or ZUNGTR and ZSTEQR.  If this fails for
      // some eigenvalue, then try DSTEBZ.

      TEST = false;
      if ( INDEIG ) {
         if ( IL == 1 && IU == N ) {
            TEST = true;
         }
      }
      if ( ( ALLEIG || TEST ) && ( ABSTOL <= ZERO ) ) {
         dcopy(N, RWORK( INDD ), 1, W, 1 );
         INDEE = INDRWK + 2*N;
         if ( !WANTZ ) {
            dcopy(N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 );
            dsterf(N, W, RWORK( INDEE ), INFO );
         } else {
            zlacpy('A', N, N, A, LDA, Z, LDZ );
            zungtr(UPLO, N, Z, LDZ, WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );
            dcopy(N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 );
            zsteqr(JOBZ, N, W, RWORK( INDEE ), Z, LDZ, RWORK( INDRWK ), INFO );
            if ( INFO == 0 ) {
               for (I = 1; I <= N; I++) { // 30
                  IFAIL[I] = 0;
               } // 30
            }
         }
         if ( INFO == 0 ) {
            M = N;
            GO TO 40;
         }
         INFO = 0;
      }

      // Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.

      if ( WANTZ ) {
         ORDER = 'B';
      } else {
         ORDER = 'E';
      }
      INDIBL = 1;
      INDISP = INDIBL + N;
      INDIWK = INDISP + N;
      dstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, RWORK( INDD ), RWORK( INDE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ), IWORK( INDIWK ), INFO );

      if ( WANTZ ) {
         zstein(N, RWORK( INDD ), RWORK( INDE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, RWORK( INDRWK ), IWORK( INDIWK ), IFAIL, INFO );

         // Apply unitary matrix used in reduction to tridiagonal
         // form to eigenvectors returned by ZSTEIN.

         zunmtr('L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), LLWORK, IINFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      } // 40
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
         for (J = 1; J <= M - 1; J++) { // 60
            I = 0;
            TMP1 = W( J );
            for (JJ = J + 1; JJ <= M; JJ++) { // 50
               if ( W( JJ ) < TMP1 ) {
                  I = JJ;
                  TMP1 = W( JJ );
               }
            } // 50

            if ( I != 0 ) {
               ITMP1 = IWORK( INDIBL+I-1 );
               W[I] = W( J );
               IWORK[INDIBL+I-1] = IWORK( INDIBL+J-1 );
               W[J] = TMP1;
               IWORK[INDIBL+J-1] = ITMP1;
               zswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
               if ( INFO != 0 ) {
                  ITMP1 = IFAIL( I );
                  IFAIL[I] = IFAIL( J );
                  IFAIL[J] = ITMP1;
               }
            }
         } // 60
      }

      // Set WORK(1) to optimal complex workspace size.

      WORK[1] = LWMIN;

      }
