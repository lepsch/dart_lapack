      SUBROUTINE SSTEVX( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO );

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE;
      int                IL, INFO, IU, LDZ, M, N;
      REAL               ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      REAL               D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, TEST, VALEIG, WANTZ;
      String             ORDER;
      int                I, IMAX, INDISP, INDIWO, INDWRK, ISCALE, ITMP1, J, JJ, NSPLIT;
      REAL               BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, TNRM, VLL, VUU;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANST;
      // EXTERNAL LSAME, SLAMCH, SLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SSCAL, SSTEBZ, SSTEIN, SSTEQR, SSTERF, SSWAP, XERBLA
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
      } else if ( N < 0 ) {
         INFO = -3;
      } else {
         if ( VALEIG ) {
            if (N > 0 && VU <= VL) INFO = -7;
         } else if ( INDEIG ) {
            if ( IL < 1 || IL > MAX( 1, N ) ) {
               INFO = -8;
            } else if ( IU < MIN( N, IL ) || IU > N ) {
               INFO = -9;
            }
         }
      }
      if ( INFO == 0 ) {
         if( LDZ < 1 || ( WANTZ && LDZ < N ) ) INFO = -14;
      }

      if ( INFO != 0 ) {
         xerbla('SSTEVX', -INFO );
         return;
      }

      // Quick return if possible

      M = 0;
      if (N == 0) RETURN;

      if ( N == 1 ) {
         if ( ALLEIG || INDEIG ) {
            M = 1;
            W( 1 ) = D( 1 );
         } else {
            if ( VL < D( 1 ) && VU >= D( 1 ) ) {
               M = 1;
               W( 1 ) = D( 1 );
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
      RMIN = SQRT( SMLNUM );
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) );

      // Scale matrix to allowable range, if necessary.

      ISCALE = 0;
      if ( VALEIG ) {
         VLL = VL;
         VUU = VU;
      } else {
         VLL = ZERO;
         VUU = ZERO;
      }
      TNRM = SLANST( 'M', N, D, E );
      if ( TNRM > ZERO && TNRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / TNRM;
      } else if ( TNRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / TNRM;
      }
      if ( ISCALE == 1 ) {
         sscal(N, SIGMA, D, 1 );
         sscal(N-1, SIGMA, E( 1 ), 1 );
         if ( VALEIG ) {
            VLL = VL*SIGMA;
            VUU = VU*SIGMA;
         }
      }

      // If all eigenvalues are desired and ABSTOL is less than zero, then
      // call SSTERF or SSTEQR.  If this fails for some eigenvalue, then
      // try SSTEBZ.

      TEST = false;
      if ( INDEIG ) {
         if ( IL == 1 && IU == N ) {
            TEST = true;
         }
      }
      if ( ( ALLEIG || TEST ) && ( ABSTOL <= ZERO ) ) {
         scopy(N, D, 1, W, 1 );
         scopy(N-1, E( 1 ), 1, WORK( 1 ), 1 );
         INDWRK = N + 1;
         if ( !WANTZ ) {
            ssterf(N, W, WORK, INFO );
         } else {
            ssteqr('I', N, W, WORK, Z, LDZ, WORK( INDWRK ), INFO );
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
      INDWRK = 1;
      INDISP = 1 + N;
      INDIWO = INDISP + N;
      sstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTOL, D, E, M, NSPLIT, W, IWORK( 1 ), IWORK( INDISP ), WORK( INDWRK ), IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         sstein(N, D, E, M, W, IWORK( 1 ), IWORK( INDISP ), Z, LDZ, WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO );
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

      // End of SSTEVX

      }
