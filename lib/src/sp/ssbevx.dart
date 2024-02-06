      void ssbevx(JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N;
      double               ABSTOL, VL, VU;
      int                IFAIL( * ), IWORK( * );
      double               AB( LDAB, * ), Q( LDQ, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               ALLEIG, INDEIG, LOWER, TEST, VALEIG, WANTZ;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, INDISP, INDIWO, INDWRK, ISCALE, ITMP1, J, JJ, NSPLIT;
      double               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, SLANSB;
      // EXTERNAL lsame, SLAMCH, SLANSB
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMV, SLACPY, SLASCL, SSBTRD, SSCAL, SSTEBZ, SSTEIN, SSTEQR, SSTERF, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      ALLEIG = lsame( RANGE, 'A' );
      VALEIG = lsame( RANGE, 'V' );
      INDEIG = lsame( RANGE, 'I' );
      LOWER = lsame( UPLO, 'L' );

      INFO = 0;
      if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
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

      if ( INFO != 0 ) {
         xerbla('SSBEVX', -INFO );
         return;
      }

      // Quick return if possible

      M = 0;
      if (N == 0) return;

      if ( N == 1 ) {
         M = 1;
         if ( LOWER ) {
            TMP1 = AB( 1, 1 );
         } else {
            TMP1 = AB( KD+1, 1 );
         }
         if ( VALEIG ) {
            if( !( VL < TMP1 && VU >= TMP1 ) ) M = 0;
         }
         if ( M == 1 ) {
            W[1] = TMP1;
            if (WANTZ) Z( 1, 1 ) = ONE;
         }
         return;
      }

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' );
      EPS = SLAMCH( 'Precision' );
      SMLNUM = SAFMIN / EPS;
      BIGNUM = ONE / SMLNUM;
      RMIN = sqrt( SMLNUM );
      RMAX = min( sqrt( BIGNUM ), ONE / sqrt( sqrt( SAFMIN ) ) );

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
      ANRM = SLANSB( 'M', UPLO, N, KD, AB, LDAB, WORK );
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if ( ISCALE == 1 ) {
         if ( LOWER ) {
            slascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            slascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
         if (ABSTOL > 0) ABSTLL = ABSTOL*SIGMA;
         if ( VALEIG ) {
            VLL = VL*SIGMA;
            VUU = VU*SIGMA;
         }
      }

      // Call SSBTRD to reduce symmetric band matrix to tridiagonal form.

      INDD = 1;
      INDE = INDD + N;
      INDWRK = INDE + N;
      ssbtrd(JOBZ, UPLO, N, KD, AB, LDAB, WORK( INDD ), WORK( INDE ), Q, LDQ, WORK( INDWRK ), IINFO );

      // If all eigenvalues are desired and ABSTOL is less than or equal
      // to zero, then call SSTERF or SSTEQR.  If this fails for some
      // eigenvalue, then try SSTEBZ.

      TEST = false;
      if (INDEIG) {
         if (IL == 1 && IU == N) {
            TEST = true;
         }
      }
      if ((ALLEIG || TEST) && (ABSTOL <= ZERO)) {
         scopy(N, WORK( INDD ), 1, W, 1 );
         INDEE = INDWRK + 2*N;
         if ( !WANTZ ) {
            scopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            ssterf(N, W, WORK( INDEE ), INFO );
         } else {
            slacpy('A', N, N, Q, LDQ, Z, LDZ );
            scopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            ssteqr(JOBZ, N, W, WORK( INDEE ), Z, LDZ, WORK( INDWRK ), INFO );
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

      // Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN.

      if ( WANTZ ) {
         ORDER = 'B';
      } else {
         ORDER = 'E';
      }
      INDIBL = 1;
      INDISP = INDIBL + N;
      INDIWO = INDISP + N;
      sstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, WORK( INDD ), WORK( INDE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWRK ), IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         sstein(N, WORK( INDD ), WORK( INDE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO );

         // Apply orthogonal matrix used in reduction to tridiagonal
         // form to eigenvectors returned by SSTEIN.

         for (J = 1; J <= M; J++) { // 20
            scopy(N, Z( 1, J ), 1, WORK( 1 ), 1 );
            sgemv('N', N, N, ONE, Q, LDQ, WORK, 1, ZERO, Z( 1, J ), 1 );
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
         sscal(IMAX, ONE / SIGMA, W, 1 );
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
               sswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
               if ( INFO != 0 ) {
                  ITMP1 = IFAIL( I );
                  IFAIL[I] = IFAIL( J );
                  IFAIL[J] = ITMP1;
               }
            }
         } // 50
      }

      return;
      }
