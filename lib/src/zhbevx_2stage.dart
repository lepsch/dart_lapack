      void zhbevx_2stage(final int JOBZ, final int RANGE, final int UPLO, final int N, final int KD, final Matrix<double> AB, final int LDAB, final Matrix<double> Q, final int LDQ, final int VL, final int VU, final int IL, final int IU, final int ABSTOL, final int M, final int W, final Matrix<double> Z, final int LDZ, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final Array<int> IWORK, final int IFAIL, final Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N, LWORK;
      double             ABSTOL, VL, VU;
      int                IFAIL( * ), IWORK( * );
      double             RWORK( * ), W( * );
      Complex         AB( LDAB, * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      bool               ALLEIG, INDEIG, LOWER, TEST, VALEIG, WANTZ, LQUERY;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, INDISP, INDIWK, INDRWK, INDWRK, ISCALE, ITMP1, LLWORK, LWMIN, LHTRD, LWTRD, IB, INDHOUS, J, JJ, NSPLIT;
      double             ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      Complex         CTMP1;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV2STAGE;
      //- double             DLAMCH, ZLANHB;
      // EXTERNAL lsame, DLAMCH, ZLANHB, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSCAL, DSTEBZ, DSTERF, XERBLA, ZCOPY, ZGEMV, ZLACPY, ZLASCL, ZSTEIN, ZSTEQR, ZSWAP, ZHETRD_HB2ST
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, SQRT

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      ALLEIG = lsame( RANGE, 'A' );
      VALEIG = lsame( RANGE, 'V' );
      INDEIG = lsame( RANGE, 'I' );
      LOWER = lsame( UPLO, 'L' );
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
      } else if ( KD < 0 ) {
         INFO = -5;
      } else if ( LDAB < KD+1 ) {
         INFO = -7;
      } else if ( WANTZ && LDQ < max( 1, N ) ) {
         INFO = -9;
      } else {
         if ( VALEIG ) {
            if (N > 0 && VU <= VL) INFO = -11;
         } else if ( INDEIG ) {
            if ( IL < 1 || IL > max( 1, N ) ) {
               INFO = -12;
            } else if ( IU < min( N, IL ) || IU > N ) {
               INFO = -13;
            }
         }
      }
      if ( INFO == 0 ) {
         if( LDZ < 1 || ( WANTZ && LDZ < N ) ) INFO = -18;
      }

      if ( INFO == 0 ) {
         if ( N <= 1 ) {
            LWMIN = 1;
            WORK[1] = LWMIN;
         } else {
            IB    = ILAENV2STAGE( 2, 'ZHETRD_HB2ST', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'ZHETRD_HB2ST', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'ZHETRD_HB2ST', JOBZ, N, KD, IB, -1 );
            LWMIN = LHTRD + LWTRD;
            WORK[1] = LWMIN;
         }

         if (LWORK < LWMIN && !LQUERY) INFO = -20;
      }

      if ( INFO != 0 ) {
         xerbla('ZHBEVX_2STAGE', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      M = 0;
      if (N == 0) return;

      if ( N == 1 ) {
         M = 1;
         if ( LOWER ) {
            CTMP1 = AB( 1, 1 );
         } else {
            CTMP1 = AB( KD+1, 1 );
         }
         TMP1 = CTMP1.toDouble();
         if ( VALEIG ) {
            if( !( VL < TMP1 && VU >= TMP1 ) ) M = 0;
         }
         if ( M == 1 ) {
            W[1] = CTMP1.toDouble();
            if (WANTZ) Z( 1, 1 ) = CONE;
         }
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
      } else {
         VLL = ZERO;
         VUU = ZERO;
      }
      ANRM = ZLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK );
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if ( ISCALE == 1 ) {
         if ( LOWER ) {
            zlascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            zlascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
         if (ABSTOL > 0) ABSTLL = ABSTOL*SIGMA;
         if ( VALEIG ) {
            VLL = VL*SIGMA;
            VUU = VU*SIGMA;
         }
      }

      // Call ZHBTRD_HB2ST to reduce Hermitian band matrix to tridiagonal form.

      INDD = 1;
      INDE = INDD + N;
      INDRWK = INDE + N;

      INDHOUS = 1;
      INDWRK  = INDHOUS + LHTRD;
      LLWORK  = LWORK - INDWRK + 1;

      zhetrd_hb2st('N', JOBZ, UPLO, N, KD, AB, LDAB, RWORK( INDD ), RWORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // If all eigenvalues are desired and ABSTOL is less than or equal
      // to zero, then call DSTERF or ZSTEQR.  If this fails for some
      // eigenvalue, then try DSTEBZ.

      TEST = false;
      if (INDEIG) {
         if (IL == 1 && IU == N) {
            TEST = true;
         }
      }
      if ((ALLEIG || TEST) && (ABSTOL <= ZERO)) {
         dcopy(N, RWORK( INDD ), 1, W, 1 );
         INDEE = INDRWK + 2*N;
         if ( !WANTZ ) {
            dcopy(N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 );
            dsterf(N, W, RWORK( INDEE ), INFO );
         } else {
            zlacpy('A', N, N, Q, LDQ, Z, LDZ );
            dcopy(N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 );
            zsteqr(JOBZ, N, W, RWORK( INDEE ), Z, LDZ, RWORK( INDRWK ), INFO );
            if ( INFO == 0 ) {
               for (I = 1; I <= N; I++) { // 10
                  IFAIL[I] = 0;
               } // 10
            }
         }
         if ( INFO == 0 ) {
            M = N;
            GO TO 30;
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

         for (J = 1; J <= M; J++) { // 20
            zcopy(N, Z( 1, J ), 1, WORK( 1 ), 1 );
            zgemv('N', N, N, CONE, Q, LDQ, WORK, 1, CZERO, Z( 1, J ), 1 );
         } // 20
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
               if ( INFO != 0 ) {
                  ITMP1 = IFAIL( I );
                  IFAIL[I] = IFAIL( J );
                  IFAIL[J] = ITMP1;
               }
            }
         } // 50
      }

      // Set WORK(1) to optimal workspace size.

      WORK[1] = LWMIN;

      }
