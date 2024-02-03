      SUBROUTINE ZHEEVR_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDA, LDZ, LIWORK, LRWORK, LWORK, M, N;
      double             ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * );
      double             RWORK( * ), W( * );
      COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LOWER, LQUERY, TEST, VALEIG, WANTZ, TRYRAC;
      String             ORDER;
      int                I, IEEEOK, IINFO, IMAX, INDIBL, INDIFL, INDISP, INDIWO, INDRD, INDRDD, INDRE, INDREE, INDRWK, INDTAU, INDWK, INDWKN, ISCALE, ITMP1, J, JJ, LIWMIN, LLWORK, LLRWORK, LLWRKN, LRWMIN, LWMIN, NSPLIT, LHTRD, LWTRD, KD, IB, INDHOUS;
      double             ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV, ILAENV2STAGE;
      double             DLAMCH, ZLANSY;
      // EXTERNAL LSAME, DLAMCH, ZLANSY, ILAENV, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSCAL, DSTEBZ, DSTERF, XERBLA, ZDSCAL, ZHETRD_2STAGE, ZSTEMR, ZSTEIN, ZSWAP, ZUNMTR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      IEEEOK = ILAENV( 10, 'ZHEEVR', 'N', 1, 2, 3, 4 )

      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )

      LQUERY = ( ( LWORK.EQ.-1 ) .OR. ( LRWORK.EQ.-1 ) .OR. ( LIWORK.EQ.-1 ) )

      KD     = ILAENV2STAGE( 1, 'ZHETRD_2STAGE', JOBZ, N, -1, -1, -1 )
      IB     = ILAENV2STAGE( 2, 'ZHETRD_2STAGE', JOBZ, N, KD, -1, -1 )
      LHTRD  = ILAENV2STAGE( 3, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 )
      LWTRD  = ILAENV2STAGE( 4, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 )

      if ( N.LE.1 ) {
         LWMIN  = 1
         LRWMIN = 1
         LIWMIN = 1
      } else {
         LWMIN  = N + LHTRD + LWTRD
         LRWMIN = 24*N
         LIWMIN = 10*N
      }

      INFO = 0
      if ( .NOT.( LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) {
         INFO = -2
      } else if ( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      } else {
         if ( VALEIG ) {
            IF( N.GT.0 .AND. VU.LE.VL ) INFO = -8
         } else if ( INDEIG ) {
            if ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) {
               INFO = -9
            } else if ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) {
               INFO = -10
            }
         }
      }
      if ( INFO.EQ.0 ) {
         if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
            INFO = -15
         }
      }

      if ( INFO.EQ.0 ) {
         WORK( 1 )  = LWMIN
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -18
         } else if ( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) {
            INFO = -20
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
            INFO = -22
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('ZHEEVR_2STAGE', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      M = 0
      if ( N.EQ.0 ) {
         WORK( 1 ) = 1
         RETURN
      }

      if ( N.EQ.1 ) {
         WORK( 1 ) = 1
         if ( ALLEIG .OR. INDEIG ) {
            M = 1
            W( 1 ) = DBLE( A( 1, 1 ) )
         } else {
            if ( VL.LT.DBLE( A( 1, 1 ) ) .AND. VU.GE.DBLE( A( 1, 1 ) ) ) {
               M = 1
               W( 1 ) = DBLE( A( 1, 1 ) )
            }
         }
         if ( WANTZ ) {
            Z( 1, 1 ) = ONE
            ISUPPZ( 1 ) = 1
            ISUPPZ( 2 ) = 1
         }
         RETURN
      }

      // Get machine constants.

      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS    = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = SQRT( SMLNUM )
      RMAX   = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )

      // Scale matrix to allowable range, if necessary.

      ISCALE = 0
      ABSTLL = ABSTOL
      if (VALEIG) {
         VLL = VL
         VUU = VU
      }
      ANRM = ZLANSY( 'M', UPLO, N, A, LDA, RWORK )
      if ( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE.EQ.1 ) {
         if ( LOWER ) {
            DO 10 J = 1, N
               zdscal(N-J+1, SIGMA, A( J, J ), 1 );
   10       CONTINUE
         } else {
            DO 20 J = 1, N
               zdscal(J, SIGMA, A( 1, J ), 1 );
   20       CONTINUE
         }
         IF( ABSTOL.GT.0 ) ABSTLL = ABSTOL*SIGMA
         if ( VALEIG ) {
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         }
      }

      // Initialize indices into workspaces.  Note: The IWORK indices are
      // used only if DSTERF or ZSTEMR fail.

      // WORK(INDTAU:INDTAU+N-1) stores the complex scalar factors of the
      // elementary reflectors used in ZHETRD.
      INDTAU = 1
      // INDWK is the starting offset of the remaining complex workspace,
      // and LLWORK is the remaining complex workspace size.
      INDHOUS = INDTAU + N
      INDWK   = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWK + 1

      // RWORK(INDRD:INDRD+N-1) stores the real tridiagonal's diagonal
      // entries.
      INDRD = 1
      // RWORK(INDRE:INDRE+N-1) stores the off-diagonal entries of the
      // tridiagonal matrix from ZHETRD.
      INDRE = INDRD + N
      // RWORK(INDRDD:INDRDD+N-1) is a copy of the diagonal entries over
      // -written by ZSTEMR (the DSTERF path copies the diagonal to W).
      INDRDD = INDRE + N
      // RWORK(INDREE:INDREE+N-1) is a copy of the off-diagonal entries over
      // -written while computing the eigenvalues in DSTERF and ZSTEMR.
      INDREE = INDRDD + N
      // INDRWK is the starting offset of the left-over real workspace, and
      // LLRWORK is the remaining workspace size.
      INDRWK = INDREE + N
      LLRWORK = LRWORK - INDRWK + 1

      // IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and
      // stores the block indices of each of the M<=N eigenvalues.
      INDIBL = 1
      // IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and
      // stores the starting and finishing indices of each block.
      INDISP = INDIBL + N
      // IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
      // that corresponding to eigenvectors that fail to converge in
      // ZSTEIN.  This information is discarded; if any fail, the driver
      // returns INFO > 0.
      INDIFL = INDISP + N
      // INDIWO is the offset of the remaining integer workspace.
      INDIWO = INDIFL + N


      // Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.

      zhetrd_2stage(JOBZ, UPLO, N, A, LDA, RWORK( INDRD ), RWORK( INDRE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWK ), LLWORK, IINFO );

      // If all eigenvalues are desired
      // then call DSTERF or ZSTEMR and ZUNMTR.

      TEST = .FALSE.
      if ( INDEIG ) {
         if ( IL.EQ.1 .AND. IU.EQ.N ) {
            TEST = .TRUE.
         }
      }
      if ( ( ALLEIG.OR.TEST ) .AND. ( IEEEOK.EQ.1 ) ) {
         if ( .NOT.WANTZ ) {
            dcopy(N, RWORK( INDRD ), 1, W, 1 );
            dcopy(N-1, RWORK( INDRE ), 1, RWORK( INDREE ), 1 );
            dsterf(N, W, RWORK( INDREE ), INFO );
         } else {
            dcopy(N-1, RWORK( INDRE ), 1, RWORK( INDREE ), 1 );
            dcopy(N, RWORK( INDRD ), 1, RWORK( INDRDD ), 1 );

            if (ABSTOL .LE. TWO*N*EPS) {
               TRYRAC = .TRUE.
            } else {
               TRYRAC = .FALSE.
            }
            zstemr(JOBZ, 'A', N, RWORK( INDRDD ), RWORK( INDREE ), VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, RWORK( INDRWK ), LLRWORK, IWORK, LIWORK, INFO );

            // Apply unitary matrix used in reduction to tridiagonal
            // form to eigenvectors returned by ZSTEMR.

            if ( WANTZ .AND. INFO.EQ.0 ) {
               INDWKN = INDWK
               LLWRKN = LWORK - INDWKN + 1
               zunmtr('L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO );
            }
         }


         if ( INFO.EQ.0 ) {
            M = N
            GO TO 30
         }
         INFO = 0
      }

      // Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.
      // Also call DSTEBZ and ZSTEIN if ZSTEMR fails.

      if ( WANTZ ) {
         ORDER = 'B'
      } else {
         ORDER = 'E'
      }
       dstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, RWORK( INDRD ), RWORK( INDRE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ), IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         zstein(N, RWORK( INDRD ), RWORK( INDRE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, RWORK( INDRWK ), IWORK( INDIWO ), IWORK( INDIFL ), INFO );

         // Apply unitary matrix used in reduction to tridiagonal
         // form to eigenvectors returned by ZSTEIN.

         INDWKN = INDWK
         LLWRKN = LWORK - INDWKN + 1
         zunmtr('L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

   30 CONTINUE
      if ( ISCALE.EQ.1 ) {
         if ( INFO.EQ.0 ) {
            IMAX = M
         } else {
            IMAX = INFO - 1
         }
         dscal(IMAX, ONE / SIGMA, W, 1 );
      }

      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.

      if ( WANTZ ) {
         DO 50 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 40 JJ = J + 1, M
               if ( W( JJ ).LT.TMP1 ) {
                  I = JJ
                  TMP1 = W( JJ )
               }
   40       CONTINUE

            if ( I.NE.0 ) {
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               zswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
            }
   50    CONTINUE
      }

      // Set WORK(1) to optimal workspace size.

      WORK( 1 )  = LWMIN
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of ZHEEVR_2STAGE

      }
