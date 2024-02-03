      SUBROUTINE DSYEVR_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N;
      double             ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * );
      double             A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * );
      // ..

* =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LOWER, LQUERY, VALEIG, WANTZ, TRYRAC;
      String             ORDER;
      int                I, IEEEOK, IINFO, IMAX, INDD, INDDD, INDE, INDEE, INDIBL, INDIFL, INDISP, INDIWO, INDTAU, INDWK, INDWKN, ISCALE, J, JJ, LIWMIN, LLWORK, LLWRKN, LWMIN, NSPLIT, LHTRD, LWTRD, KD, IB, INDHOUS;
      double             ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV, ILAENV2STAGE;
      double             DLAMCH, DLANSY;
      // EXTERNAL LSAME, DLAMCH, DLANSY, ILAENV, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DORMTR, DSCAL, DSTEBZ, DSTEMR, DSTEIN, DSTERF, DSWAP, DSYTRD_2STAGE, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      IEEEOK = ILAENV( 10, 'DSYEVR', 'N', 1, 2, 3, 4 )

      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )

      LQUERY = ( ( LWORK.EQ.-1 ) .OR. ( LIWORK.EQ.-1 ) )

      KD     = ILAENV2STAGE( 1, 'DSYTRD_2STAGE', JOBZ, N, -1, -1, -1 )
      IB     = ILAENV2STAGE( 2, 'DSYTRD_2STAGE', JOBZ, N, KD, -1, -1 )
      LHTRD  = ILAENV2STAGE( 3, 'DSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )
      LWTRD  = ILAENV2STAGE( 4, 'DSYTRD_2STAGE', JOBZ, N, KD, IB, -1 )

      if ( N.LE.1 ) {
         LWMIN  = 1
         LIWMIN = 1
      } else {
         LWMIN  = MAX( 26*N, 5*N + LHTRD + LWTRD )
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
         } else if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -18
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
            INFO = -20
         }
      }

      if ( INFO.EQ.0 ) {
          // NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
          // NB = MAX( NB, ILAENV( 1, 'DORMTR', UPLO, N, -1, -1, -1 ) )
          // LWKOPT = MAX( ( NB+1 )*N, LWMIN )
         WORK( 1 )  = LWMIN
         IWORK( 1 ) = LIWMIN
      }

      if ( INFO.NE.0 ) {
         xerbla('DSYEVR_2STAGE', -INFO );
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
            W( 1 ) = A( 1, 1 )
         } else {
            if ( VL.LT.A( 1, 1 ) .AND. VU.GE.A( 1, 1 ) ) {
               M = 1
               W( 1 ) = A( 1, 1 )
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
      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
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
               dscal(N-J+1, SIGMA, A( J, J ), 1 );
   10       CONTINUE
         } else {
            DO 20 J = 1, N
               dscal(J, SIGMA, A( 1, J ), 1 );
   20       CONTINUE
         }
         IF( ABSTOL.GT.0 ) ABSTLL = ABSTOL*SIGMA
         if ( VALEIG ) {
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         }
      }

      // Initialize indices into workspaces.  Note: The IWORK indices are
      // used only if DSTERF or DSTEMR fail.

      // WORK(INDTAU:INDTAU+N-1) stores the scalar factors of the
      // elementary reflectors used in DSYTRD.
      INDTAU = 1
      // WORK(INDD:INDD+N-1) stores the tridiagonal's diagonal entries.
      INDD = INDTAU + N
      // WORK(INDE:INDE+N-1) stores the off-diagonal entries of the
      // tridiagonal matrix from DSYTRD.
      INDE = INDD + N
      // WORK(INDDD:INDDD+N-1) is a copy of the diagonal entries over
      // -written by DSTEMR (the DSTERF path copies the diagonal to W).
      INDDD = INDE + N
      // WORK(INDEE:INDEE+N-1) is a copy of the off-diagonal entries over
      // -written while computing the eigenvalues in DSTERF and DSTEMR.
      INDEE = INDDD + N
      // INDHOUS is the starting offset Householder storage of stage 2
      INDHOUS = INDEE + N
      // INDWK is the starting offset of the left-over workspace, and
      // LLWORK is the remaining workspace size.
      INDWK  = INDHOUS + LHTRD
      LLWORK = LWORK - INDWK + 1


      // IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and
      // stores the block indices of each of the M<=N eigenvalues.
      INDIBL = 1
      // IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and
      // stores the starting and finishing indices of each block.
      INDISP = INDIBL + N
      // IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
      // that corresponding to eigenvectors that fail to converge in
      // DSTEIN.  This information is discarded; if any fail, the driver
      // returns INFO > 0.
      INDIFL = INDISP + N
      // INDIWO is the offset of the remaining integer workspace.
      INDIWO = INDIFL + N


      // Call DSYTRD_2STAGE to reduce symmetric matrix to tridiagonal form.


      dsytrd_2stage(JOBZ, UPLO, N, A, LDA, WORK( INDD ), WORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWK ), LLWORK, IINFO );

      // If all eigenvalues are desired
      // then call DSTERF or DSTEMR and DORMTR.

      if ( ( ALLEIG .OR. ( INDEIG .AND. IL.EQ.1 .AND. IU.EQ.N ) ) .AND. IEEEOK.EQ.1 ) {
         if ( .NOT.WANTZ ) {
            dcopy(N, WORK( INDD ), 1, W, 1 );
            dcopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            dsterf(N, W, WORK( INDEE ), INFO );
         } else {
            dcopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            dcopy(N, WORK( INDD ), 1, WORK( INDDD ), 1 );

            if (ABSTOL .LE. TWO*N*EPS) {
               TRYRAC = .TRUE.
            } else {
               TRYRAC = .FALSE.
            }
            dstemr(JOBZ, 'A', N, WORK( INDDD ), WORK( INDEE ), VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK( INDWK ), LWORK, IWORK, LIWORK, INFO );



         // Apply orthogonal matrix used in reduction to tridiagonal
         // form to eigenvectors returned by DSTEMR.

            if ( WANTZ .AND. INFO.EQ.0 ) {
               INDWKN = INDE
               LLWRKN = LWORK - INDWKN + 1
               dormtr('L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO );
            }
         }


         if ( INFO.EQ.0 ) {
            // Everything worked.  Skip DSTEBZ/DSTEIN.  IWORK(:) are
            // undefined.
            M = N
            GO TO 30
         }
         INFO = 0
      }

      // Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN.
      // Also call DSTEBZ and DSTEIN if DSTEMR fails.

      if ( WANTZ ) {
         ORDER = 'B'
      } else {
         ORDER = 'E'
      }
       dstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, WORK( INDD ), WORK( INDE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWK ), IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         dstein(N, WORK( INDD ), WORK( INDE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, WORK( INDWK ), IWORK( INDIWO ), IWORK( INDIFL ), INFO );

         // Apply orthogonal matrix used in reduction to tridiagonal
         // form to eigenvectors returned by DSTEIN.

         INDWKN = INDE
         LLWRKN = LWORK - INDWKN + 1
         dormtr('L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

*  Jump here if DSTEMR/DSTEIN succeeded.
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
      // eigenvectors.  Note: We do not sort the IFAIL portion of IWORK.
      // It may not be initialized (if DSTEMR/DSTEIN succeeded), and we do
      // not return this detailed information to the user.

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
               W( I ) = W( J )
               W( J ) = TMP1
               dswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
            }
   50    CONTINUE
      }

      // Set WORK(1) to optimal workspace size.

      WORK( 1 )  = LWMIN
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of DSYEVR_2STAGE

      }
