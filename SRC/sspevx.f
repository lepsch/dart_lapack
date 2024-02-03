      void sspevx(JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDZ, M, N;
      REAL               ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      REAL               AP( * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, TEST, VALEIG, WANTZ;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDISP, INDIWO, INDTAU, INDWRK, ISCALE, ITMP1, J, JJ, NSPLIT;
      REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANSP;
      // EXTERNAL LSAME, SLAMCH, SLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SOPGTR, SOPMTR, SSCAL, SSPTRD, SSTEBZ, SSTEIN, SSTEQR, SSTERF, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' );
      ALLEIG = LSAME( RANGE, 'A' );
      VALEIG = LSAME( RANGE, 'V' );
      INDEIG = LSAME( RANGE, 'I' );

      INFO = 0;
      if ( !( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -2;
      } else if ( !( LSAME( UPLO, 'L' ) || LSAME( UPLO, 'U' ) ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else {
         if ( VALEIG ) {
            if (N > 0 && VU <= VL) INFO = -7;
         } else if ( INDEIG ) {
            if ( IL < 1 || IL > max( 1, N ) ) {
               INFO = -8;
            } else if ( IU < min( N, IL ) || IU > N ) {
               INFO = -9;
            }
         }
      }
      if ( INFO == 0 ) {
         if( LDZ < 1 || ( WANTZ && LDZ < N ) ) INFO = -14;
      }

      if ( INFO != 0 ) {
         xerbla('SSPEVX', -INFO );
         return;
      }

      // Quick return if possible

      M = 0;
      if (N == 0) return;

      if ( N == 1 ) {
         if ( ALLEIG || INDEIG ) {
            M = 1;
            W( 1 ) = AP( 1 );
         } else {
            if ( VL < AP( 1 ) && VU >= AP( 1 ) ) {
               M = 1;
               W( 1 ) = AP( 1 );
            }
         }
         if (WANTZ) Z( 1, 1 ) = ONE;
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
      ANRM = SLANSP( 'M', UPLO, N, AP, WORK );
      if ( ANRM > ZERO && ANRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / ANRM;
      } else if ( ANRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / ANRM;
      }
      if ( ISCALE == 1 ) {
         sscal(( N*( N+1 ) ) / 2, SIGMA, AP, 1 );
         if (ABSTOL > 0) ABSTLL = ABSTOL*SIGMA;
         if ( VALEIG ) {
            VLL = VL*SIGMA;
            VUU = VU*SIGMA;
         }
      }

      // Call SSPTRD to reduce symmetric packed matrix to tridiagonal form.

      INDTAU = 1;
      INDE = INDTAU + N;
      INDD = INDE + N;
      INDWRK = INDD + N;
      ssptrd(UPLO, N, AP, WORK( INDD ), WORK( INDE ), WORK( INDTAU ), IINFO );

      // If all eigenvalues are desired and ABSTOL is less than or equal
      // to zero, then call SSTERF or SOPGTR and SSTEQR.  If this fails
      // for some eigenvalue, then try SSTEBZ.

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
            sopgtr(UPLO, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
            scopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            ssteqr(JOBZ, N, W, WORK( INDEE ), Z, LDZ, WORK( INDWRK ), INFO );
            if ( INFO == 0 ) {
               for (I = 1; I <= N; I++) { // 10
                  IFAIL( I ) = 0;
               } // 10
            }
         }
         if ( INFO == 0 ) {
            M = N;
            GO TO 20;
         }
         INFO = 0;
      }

      // Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN.

      if ( WANTZ ) {
         ORDER = 'B';
      } else {
         ORDER = 'E';
      }
      INDISP = 1 + N;
      INDIWO = INDISP + N;
      sstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, WORK( INDD ), WORK( INDE ), M, NSPLIT, W, IWORK( 1 ), IWORK( INDISP ), WORK( INDWRK ), IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         sstein(N, WORK( INDD ), WORK( INDE ), M, W, IWORK( 1 ), IWORK( INDISP ), Z, LDZ, WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO );

         // Apply orthogonal matrix used in reduction to tridiagonal
         // form to eigenvectors returned by SSTEIN.

         sopmtr('L', UPLO, 'N', N, M, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      } // 20
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
         for (J = 1; J <= M - 1; J++) { // 40
            I = 0;
            TMP1 = W( J );
            for (JJ = J + 1; JJ <= M; JJ++) { // 30
               if ( W( JJ ) < TMP1 ) {
                  I = JJ;
                  TMP1 = W( JJ );
               }
            } // 30

            if ( I != 0 ) {
               ITMP1 = IWORK( 1 + I-1 );
               W( I ) = W( J );
               IWORK( 1 + I-1 ) = IWORK( 1 + J-1 );
               W( J ) = TMP1;
               IWORK( 1 + J-1 ) = ITMP1;
               sswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
               if ( INFO != 0 ) {
                  ITMP1 = IFAIL( I );
                  IFAIL( I ) = IFAIL( J );
                  IFAIL( J ) = ITMP1;
               }
            }
         } // 40
      }

      return;

      // End of SSPEVX

      }
