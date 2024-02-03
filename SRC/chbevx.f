      SUBROUTINE CHBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO );

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N;
      REAL               ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      REAL               RWORK( * ), W( * );
      COMPLEX            AB( LDAB, * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LOWER, TEST, VALEIG, WANTZ;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, INDISP, INDIWK, INDRWK, INDWRK, ISCALE, ITMP1, J, JJ, NSPLIT;
      REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      COMPLEX            CTMP1;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANHB, SLAMCH;
      // EXTERNAL LSAME, CLANHB, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMV, CHBTRD, CLACPY, CLASCL, CSTEIN, CSTEQR, CSWAP, SCOPY, SSCAL, SSTEBZ, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' );
      ALLEIG = LSAME( RANGE, 'A' );
      VALEIG = LSAME( RANGE, 'V' );
      INDEIG = LSAME( RANGE, 'I' );
      LOWER = LSAME( UPLO, 'L' );

      INFO = 0;
      if ( !( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -2;
      } else if ( !( LOWER || LSAME( UPLO, 'U' ) ) ) {
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
         xerbla('CHBEVX', -INFO );
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
         TMP1 = REAL( CTMP1 );
         if ( VALEIG ) {
            if( !( VL < TMP1 && VU >= TMP1 ) ) M = 0;
         }
         if ( M == 1 ) {
            W( 1 ) = REAL( CTMP1 );
            if (WANTZ) Z( 1, 1 ) = CONE;
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
      ANRM = CLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK );
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if ( ISCALE == 1 ) {
         if ( LOWER ) {
            clascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            clascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
         if (ABSTOL > 0) ABSTLL = ABSTOL*SIGMA;
         if ( VALEIG ) {
            VLL = VL*SIGMA;
            VUU = VU*SIGMA;
         }
      }

      // Call CHBTRD to reduce Hermitian band matrix to tridiagonal form.

      INDD = 1;
      INDE = INDD + N;
      INDRWK = INDE + N;
      INDWRK = 1;
      chbtrd(JOBZ, UPLO, N, KD, AB, LDAB, RWORK( INDD ), RWORK( INDE ), Q, LDQ, WORK( INDWRK ), IINFO );

      // If all eigenvalues are desired and ABSTOL is less than or equal
      // to zero, then call SSTERF or CSTEQR.  If this fails for some
      // eigenvalue, then try SSTEBZ.

      TEST = false;
      if (INDEIG) {
         if (IL == 1 && IU == N) {
            TEST = true;
         }
      }
      if ((ALLEIG || TEST) && (ABSTOL <= ZERO)) {
         scopy(N, RWORK( INDD ), 1, W, 1 );
         INDEE = INDRWK + 2*N;
         if ( !WANTZ ) {
            scopy(N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 );
            ssterf(N, W, RWORK( INDEE ), INFO );
         } else {
            clacpy('A', N, N, Q, LDQ, Z, LDZ );
            scopy(N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 );
            csteqr(JOBZ, N, W, RWORK( INDEE ), Z, LDZ, RWORK( INDRWK ), INFO );
            if ( INFO == 0 ) {
               for (I = 1; I <= N; I++) { // 10
                  IFAIL( I ) = 0;
               } // 10
            }
         }
         if ( INFO == 0 ) {
            M = N;
            GO TO 30;
         }
         INFO = 0;
      }

      // Otherwise, call SSTEBZ and, if eigenvectors are desired, CSTEIN.

      if ( WANTZ ) {
         ORDER = 'B';
      } else {
         ORDER = 'E';
      }
      INDIBL = 1;
      INDISP = INDIBL + N;
      INDIWK = INDISP + N;
      sstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, RWORK( INDD ), RWORK( INDE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ), IWORK( INDIWK ), INFO );

      if ( WANTZ ) {
         cstein(N, RWORK( INDD ), RWORK( INDE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, RWORK( INDRWK ), IWORK( INDIWK ), IFAIL, INFO );

         // Apply unitary matrix used in reduction to tridiagonal
         // form to eigenvectors returned by CSTEIN.

         for (J = 1; J <= M; J++) { // 20
            ccopy(N, Z( 1, J ), 1, WORK( 1 ), 1 );
            cgemv('N', N, N, CONE, Q, LDQ, WORK, 1, CZERO, Z( 1, J ), 1 );
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
               W( I ) = W( J );
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 );
               W( J ) = TMP1;
               IWORK( INDIBL+J-1 ) = ITMP1;
               cswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
               if ( INFO != 0 ) {
                  ITMP1 = IFAIL( I );
                  IFAIL( I ) = IFAIL( J );
                  IFAIL( J ) = ITMP1;
               }
            }
         } // 50
      }

      return;

      // End of CHBEVX

      }
