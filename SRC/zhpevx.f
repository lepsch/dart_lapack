      SUBROUTINE ZHPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO )

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
      double             RWORK( * ), W( * );
      COMPLEX*16         AP( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D0, 0.0D0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, TEST, VALEIG, WANTZ;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDISP, INDIWK, INDRWK, INDTAU, INDWRK, ISCALE, ITMP1, J, JJ, NSPLIT;
      double             ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANHP;
      // EXTERNAL LSAME, DLAMCH, ZLANHP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSCAL, DSTEBZ, DSTERF, XERBLA, ZDSCAL, ZHPTRD, ZSTEIN, ZSTEQR, ZSWAP, ZUPGTR, ZUPMTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )

      INFO = 0
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) {
         INFO = -2
      } else if ( .NOT.( LSAME( UPLO, 'L' ) .OR. LSAME( UPLO, 'U' ) ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else {
         if ( VALEIG ) {
            if (N.GT.0 && VU.LE.VL) INFO = -7;
         } else if ( INDEIG ) {
            if ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) {
               INFO = -8
            } else if ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) {
               INFO = -9
            }
         }
      }
      if ( INFO == 0 ) {
         IF( LDZ.LT.1 .OR. ( WANTZ && LDZ.LT.N ) ) INFO = -14
      }

      if ( INFO != 0 ) {
         xerbla('ZHPEVX', -INFO );
         RETURN
      }

      // Quick return if possible

      M = 0
      if (N == 0) RETURN;

      if ( N == 1 ) {
         if ( ALLEIG .OR. INDEIG ) {
            M = 1
            W( 1 ) = DBLE( AP( 1 ) )
         } else {
            if ( VL.LT.DBLE( AP( 1 ) ) && VU.GE.DBLE( AP( 1 ) ) ) {
               M = 1
               W( 1 ) = DBLE( AP( 1 ) )
            }
         }
         if (WANTZ) Z( 1, 1 ) = CONE;
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
      ANRM = ZLANHP( 'M', UPLO, N, AP, RWORK )
      if ( ANRM.GT.ZERO && ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE == 1 ) {
         zdscal(( N*( N+1 ) ) / 2, SIGMA, AP, 1 );
         if (ABSTOL.GT.0) ABSTLL = ABSTOL*SIGMA;
         if ( VALEIG ) {
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         }
      }

      // Call ZHPTRD to reduce Hermitian packed matrix to tridiagonal form.

      INDD = 1
      INDE = INDD + N
      INDRWK = INDE + N
      INDTAU = 1
      INDWRK = INDTAU + N
      zhptrd(UPLO, N, AP, RWORK( INDD ), RWORK( INDE ), WORK( INDTAU ), IINFO );

      // If all eigenvalues are desired and ABSTOL is less than or equal
      // to zero, then call DSTERF or ZUPGTR and ZSTEQR.  If this fails
      // for some eigenvalue, then try DSTEBZ.

      TEST = false;
      if (INDEIG) {
         if (IL == 1 && IU == N) {
            TEST = true;
         }
      }
      if ((ALLEIG .OR. TEST) && (ABSTOL.LE.ZERO)) {
         dcopy(N, RWORK( INDD ), 1, W, 1 );
         INDEE = INDRWK + 2*N
         if ( .NOT.WANTZ ) {
            dcopy(N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 );
            dsterf(N, W, RWORK( INDEE ), INFO );
         } else {
            zupgtr(UPLO, N, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
            dcopy(N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 );
            zsteqr(JOBZ, N, W, RWORK( INDEE ), Z, LDZ, RWORK( INDRWK ), INFO );
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

      // Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.

      if ( WANTZ ) {
         ORDER = 'B'
      } else {
         ORDER = 'E'
      }
      INDISP = 1 + N
      INDIWK = INDISP + N
      dstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, RWORK( INDD ), RWORK( INDE ), M, NSPLIT, W, IWORK( 1 ), IWORK( INDISP ), RWORK( INDRWK ), IWORK( INDIWK ), INFO );

      if ( WANTZ ) {
         zstein(N, RWORK( INDD ), RWORK( INDE ), M, W, IWORK( 1 ), IWORK( INDISP ), Z, LDZ, RWORK( INDRWK ), IWORK( INDIWK ), IFAIL, INFO );

         // Apply unitary matrix used in reduction to tridiagonal
         // form to eigenvectors returned by ZSTEIN.

         INDWRK = INDTAU + N
         zupmtr('L', UPLO, 'N', N, M, AP, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), IINFO );
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
               if ( W( JJ ).LT.TMP1 ) {
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
               zswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
               if ( INFO != 0 ) {
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               }
            }
         } // 40
      }

      RETURN

      // End of ZHPEVX

      }
