      void dstevr(JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE;
      int                IL, INFO, IU, LDZ, LIWORK, LWORK, M, N;
      double             ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * );
      double             D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, TEST, LQUERY, VALEIG, WANTZ, TRYRAC;
      String             ORDER;
      int                I, IEEEOK, IMAX, INDIBL, INDIFL, INDISP, INDIWO, ISCALE, ITMP1, J, JJ, LIWMIN, LWMIN, NSPLIT;
      double             BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, TNRM, VLL, VUU;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- double             DLAMCH, DLANST;
      // EXTERNAL lsame, ILAENV, DLAMCH, DLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSCAL, DSTEBZ, DSTEMR, DSTEIN, DSTERF, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..


      // Test the input parameters.

      IEEEOK = ILAENV( 10, 'DSTEVR', 'N', 1, 2, 3, 4 );

      WANTZ = lsame( JOBZ, 'V' );
      ALLEIG = lsame( RANGE, 'A' );
      VALEIG = lsame( RANGE, 'V' );
      INDEIG = lsame( RANGE, 'I' );

      LQUERY = ( ( LWORK == -1 ) || ( LIWORK == -1 ) );
      LWMIN = max( 1, 20*N );
      LIWMIN = max( 1, 10*N );


      INFO = 0;
      if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
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
         if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
            INFO = -14;
         }
      }

      if ( INFO == 0 ) {
         WORK[1] = LWMIN;
         IWORK[1] = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -17;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -19;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DSTEVR', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      M = 0;
      if (N == 0) return;

      if ( N == 1 ) {
         if ( ALLEIG || INDEIG ) {
            M = 1;
            W[1] = D( 1 );
         } else {
            if ( VL < D( 1 ) && VU >= D( 1 ) ) {
               M = 1;
               W[1] = D( 1 );
            }
         }
         if (WANTZ) Z( 1, 1 ) = ONE;
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
      if ( VALEIG ) {
         VLL = VL;
         VUU = VU;
      }

      TNRM = DLANST( 'M', N, D, E );
      if ( TNRM > ZERO && TNRM < RMIN ) {
         ISCALE = 1;
         SIGMA = RMIN / TNRM;
      } else if ( TNRM > RMAX ) {
         ISCALE = 1;
         SIGMA = RMAX / TNRM;
      }
      if ( ISCALE == 1 ) {
         dscal(N, SIGMA, D, 1 );
         dscal(N-1, SIGMA, E( 1 ), 1 );
         if ( VALEIG ) {
            VLL = VL*SIGMA;
            VUU = VU*SIGMA;
         }
      }

      // Initialize indices into workspaces.  Note: These indices are used only
      // if DSTERF or DSTEMR fail.

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
      INDIWO = INDISP + N;

      // If all eigenvalues are desired, then
      // call DSTERF or DSTEMR.  If this fails for some eigenvalue, then
      // try DSTEBZ.


      TEST = false;
      if ( INDEIG ) {
         if ( IL == 1 && IU == N ) {
            TEST = true;
         }
      }
      if ( ( ALLEIG || TEST ) && IEEEOK == 1 ) {
         dcopy(N-1, E( 1 ), 1, WORK( 1 ), 1 );
         if ( !WANTZ ) {
            dcopy(N, D, 1, W, 1 );
            dsterf(N, W, WORK, INFO );
         } else {
            dcopy(N, D, 1, WORK( N+1 ), 1 );
            if (ABSTOL <= TWO*N*EPS) {
               TRYRAC = true;
            } else {
               TRYRAC = false;
            }
            dstemr(JOBZ, 'A', N, WORK( N+1 ), WORK, VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK( 2*N+1 ), LWORK-2*N, IWORK, LIWORK, INFO );

         }
         if ( INFO == 0 ) {
            M = N;
            GO TO 10;
         }
         INFO = 0;
      }

      // Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN.

      if ( WANTZ ) {
         ORDER = 'B';
      } else {
         ORDER = 'E';
      }
       dstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTOL, D, E, M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), WORK, IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         dstein(N, D, E, M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, WORK, IWORK( INDIWO ), IWORK( INDIFL ), INFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      } // 10
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
         for (J = 1; J <= M - 1; J++) { // 30
            I = 0;
            TMP1 = W( J );
            for (JJ = J + 1; JJ <= M; JJ++) { // 20
               if ( W( JJ ) < TMP1 ) {
                  I = JJ;
                  TMP1 = W( JJ );
               }
            } // 20

            if ( I != 0 ) {
               ITMP1 = IWORK( I );
               W[I] = W( J );
               IWORK[I] = IWORK( J );
               W[J] = TMP1;
               IWORK[J] = ITMP1;
               dswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
            }
         } // 30
      }

       // Causes problems with tests 19 & 20:
       // IF (wantz && INDEIG ) Z( 1,1) = Z(1,1) / 1.002 + .002


      WORK[1] = LWMIN;
      IWORK[1] = LIWMIN;
      return;
      }
