      SUBROUTINE DSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDZ, M, N;
      double             ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      double             AP( * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, TEST, VALEIG, WANTZ;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDISP, INDIWO, INDTAU, INDWRK, ISCALE, ITMP1, J, JJ, NSPLIT;
      double             ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANSP;
      // EXTERNAL LSAME, DLAMCH, DLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DOPGTR, DOPMTR, DSCAL, DSPTRD, DSTEBZ, DSTEIN, DSTEQR, DSTERF, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )

      INFO = 0
      if ( .NOT.( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -2
      } else if ( .NOT.( LSAME( UPLO, 'L' ) || LSAME( UPLO, 'U' ) ) ) {
         INFO = -3
      } else if ( N < 0 ) {
         INFO = -4
      } else {
         if ( VALEIG ) {
            if (N.GT.0 && VU.LE.VL) INFO = -7;
         } else if ( INDEIG ) {
            if ( IL < 1 || IL.GT.MAX( 1, N ) ) {
               INFO = -8
            } else if ( IU < MIN( N, IL ) || IU.GT.N ) {
               INFO = -9
            }
         }
      }
      if ( INFO == 0 ) {
         IF( LDZ < 1 || ( WANTZ && LDZ < N ) ) INFO = -14
      }

      if ( INFO != 0 ) {
         xerbla('DSPEVX', -INFO );
         RETURN
      }

      // Quick return if possible

      M = 0
      if (N == 0) RETURN;

      if ( N == 1 ) {
         if ( ALLEIG || INDEIG ) {
            M = 1
            W( 1 ) = AP( 1 )
         } else {
            if ( VL < AP( 1 ) && VU.GE.AP( 1 ) ) {
               M = 1
               W( 1 ) = AP( 1 )
            }
         }
         if (WANTZ) Z( 1, 1 ) = ONE;
         RETURN
      }

      // Get machine constants.

      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )

      // Scale matrix to allowable range, if necessary.

      ISCALE = 0
      ABSTLL = ABSTOL
      if ( VALEIG ) {
         VLL = VL
         VUU = VU
      } else {
         VLL = ZERO
         VUU = ZERO
      }
      ANRM = DLANSP( 'M', UPLO, N, AP, WORK )
      if ( ANRM.GT.ZERO && ANRM < RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE == 1 ) {
         dscal(( N*( N+1 ) ) / 2, SIGMA, AP, 1 );
         if (ABSTOL.GT.0) ABSTLL = ABSTOL*SIGMA;
         if ( VALEIG ) {
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         }
      }

      // Call DSPTRD to reduce symmetric packed matrix to tridiagonal form.

      INDTAU = 1
      INDE = INDTAU + N
      INDD = INDE + N
      INDWRK = INDD + N
      dsptrd(UPLO, N, AP, WORK( INDD ), WORK( INDE ), WORK( INDTAU ), IINFO );

      // If all eigenvalues are desired and ABSTOL is less than or equal
      // to zero, then call DSTERF or DOPGTR and SSTEQR.  If this fails
      // for some eigenvalue, then try DSTEBZ.

      TEST = false;
      if (INDEIG) {
         if (IL == 1 && IU == N) {
            TEST = true;
         }
      }
      if ((ALLEIG || TEST) && (ABSTOL.LE.ZERO)) {
         dcopy(N, WORK( INDD ), 1, W, 1 );
         INDEE = INDWRK + 2*N
         if ( .NOT.WANTZ ) {
            dcopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            dsterf(N, W, WORK( INDEE ), INFO );
         } else {
            dopgtr(UPLO, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
            dcopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            dsteqr(JOBZ, N, W, WORK( INDEE ), Z, LDZ, WORK( INDWRK ), INFO );
            if ( INFO == 0 ) {
               for (I = 1; I <= N; I++) { // 10
                  IFAIL( I ) = 0
               } // 10
            }
         }
         if ( INFO == 0 ) {
            M = N
            GO TO 20
         }
         INFO = 0
      }

      // Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN.

      if ( WANTZ ) {
         ORDER = 'B'
      } else {
         ORDER = 'E'
      }
      INDISP = 1 + N
      INDIWO = INDISP + N
      dstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, WORK( INDD ), WORK( INDE ), M, NSPLIT, W, IWORK( 1 ), IWORK( INDISP ), WORK( INDWRK ), IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         dstein(N, WORK( INDD ), WORK( INDE ), M, W, IWORK( 1 ), IWORK( INDISP ), Z, LDZ, WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO );

         // Apply orthogonal matrix used in reduction to tridiagonal
         // form to eigenvectors returned by DSTEIN.

         dopmtr('L', UPLO, 'N', N, M, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      } // 20
      if ( ISCALE == 1 ) {
         if ( INFO == 0 ) {
            IMAX = M
         } else {
            IMAX = INFO - 1
         }
         dscal(IMAX, ONE / SIGMA, W, 1 );
      }

      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.

      if ( WANTZ ) {
         for (J = 1; J <= M - 1; J++) { // 40
            I = 0
            TMP1 = W( J )
            for (JJ = J + 1; JJ <= M; JJ++) { // 30
               if ( W( JJ ) < TMP1 ) {
                  I = JJ
                  TMP1 = W( JJ )
               }
            } // 30

            if ( I != 0 ) {
               ITMP1 = IWORK( 1 + I-1 )
               W( I ) = W( J )
               IWORK( 1 + I-1 ) = IWORK( 1 + J-1 )
               W( J ) = TMP1
               IWORK( 1 + J-1 ) = ITMP1
               dswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
               if ( INFO != 0 ) {
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               }
            }
         } // 40
      }

      RETURN

      // End of DSPEVX

      }
