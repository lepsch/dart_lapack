      SUBROUTINE SSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDA, LDZ, LWORK, M, N;
      REAL               ABSTOL, VL, VU
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      REAL               A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

* =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LOWER, LQUERY, TEST, VALEIG, WANTZ;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, INDISP, INDIWO, INDTAU, INDWKN, INDWRK, ISCALE, ITMP1, J, JJ, LLWORK, LLWRKN, LWKMIN, LWKOPT, NB, NSPLIT;
      REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANSY, SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SLAMCH, SLANSY, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLACPY, SORGTR, SORMTR, SSCAL, SSTEBZ, SSTEIN, SSTEQR, SSTERF, SSWAP, SSYTRD, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
      LQUERY = ( LWORK.EQ.-1 )

      INFO = 0
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
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
            if (N.GT.0 .AND. VU.LE.VL) INFO = -8;
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
         if ( N.LE.1 ) {
            LWKMIN = 1
            LWKOPT = 1
         } else {
            LWKMIN = 8*N
            NB = ILAENV( 1, 'SSYTRD', UPLO, N, -1, -1, -1 )
            NB = MAX( NB, ILAENV( 1, 'SORMTR', UPLO, N, -1, -1, -1 ) )
            LWKOPT = MAX( LWKMIN, ( NB + 3 )*N )
         }
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

         if (LWORK.LT.LWKMIN .AND. .NOT.LQUERY) INFO = -17;
      }

      if ( INFO.NE.0 ) {
         xerbla('SSYEVX', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      M = 0
      if ( N.EQ.0 ) {
         RETURN
      }

      if ( N.EQ.1 ) {
         if ( ALLEIG .OR. INDEIG ) {
            M = 1
            W( 1 ) = A( 1, 1 )
         } else {
            if ( VL.LT.A( 1, 1 ) .AND. VU.GE.A( 1, 1 ) ) {
               M = 1
               W( 1 ) = A( 1, 1 )
            }
         }
         if (WANTZ) Z( 1, 1 ) = ONE;
         RETURN
      }

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
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
      }
      ANRM = SLANSY( 'M', UPLO, N, A, LDA, WORK )
      if ( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE.EQ.1 ) {
         if ( LOWER ) {
            for (J = 1; J <= N; J++) { // 10
               sscal(N-J+1, SIGMA, A( J, J ), 1 );
            } // 10
         } else {
            for (J = 1; J <= N; J++) { // 20
               sscal(J, SIGMA, A( 1, J ), 1 );
            } // 20
         }
         if (ABSTOL.GT.0) ABSTLL = ABSTOL*SIGMA;
         if ( VALEIG ) {
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         }
      }

      // Call SSYTRD to reduce symmetric matrix to tridiagonal form.

      INDTAU = 1
      INDE = INDTAU + N
      INDD = INDE + N
      INDWRK = INDD + N
      LLWORK = LWORK - INDWRK + 1
      ssytrd(UPLO, N, A, LDA, WORK( INDD ), WORK( INDE ), WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );

      // If all eigenvalues are desired and ABSTOL is less than or equal to
      // zero, then call SSTERF or SORGTR and SSTEQR.  If this fails for
      // some eigenvalue, then try SSTEBZ.

      TEST = .FALSE.
      if ( INDEIG ) {
         if ( IL.EQ.1 .AND. IU.EQ.N ) {
            TEST = .TRUE.
         }
      }
      if ( ( ALLEIG .OR. TEST ) .AND. ( ABSTOL.LE.ZERO ) ) {
         scopy(N, WORK( INDD ), 1, W, 1 );
         INDEE = INDWRK + 2*N
         if ( .NOT.WANTZ ) {
            scopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            ssterf(N, W, WORK( INDEE ), INFO );
         } else {
            slacpy('A', N, N, A, LDA, Z, LDZ );
            sorgtr(UPLO, N, Z, LDZ, WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO );
            scopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            ssteqr(JOBZ, N, W, WORK( INDEE ), Z, LDZ, WORK( INDWRK ), INFO );
            if ( INFO.EQ.0 ) {
               for (I = 1; I <= N; I++) { // 30
                  IFAIL( I ) = 0
               } // 30
            }
         }
         if ( INFO.EQ.0 ) {
            M = N
            GO TO 40
         }
         INFO = 0
      }

      // Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN.

      if ( WANTZ ) {
         ORDER = 'B'
      } else {
         ORDER = 'E'
      }
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWO = INDISP + N
      sstebz(RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, WORK( INDD ), WORK( INDE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWRK ), IWORK( INDIWO ), INFO );

      if ( WANTZ ) {
         sstein(N, WORK( INDD ), WORK( INDE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO );

         // Apply orthogonal matrix used in reduction to tridiagonal
         // form to eigenvectors returned by SSTEIN.

         INDWKN = INDE
         LLWRKN = LWORK - INDWKN + 1
         sormtr('L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWKN ), LLWRKN, IINFO );
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      } // 40
      if ( ISCALE.EQ.1 ) {
         if ( INFO.EQ.0 ) {
            IMAX = M
         } else {
            IMAX = INFO - 1
         }
         sscal(IMAX, ONE / SIGMA, W, 1 );
      }

      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.

      if ( WANTZ ) {
         for (J = 1; J <= M - 1; J++) { // 60
            I = 0
            TMP1 = W( J )
            for (JJ = J + 1; JJ <= M; JJ++) { // 50
               if ( W( JJ ).LT.TMP1 ) {
                  I = JJ
                  TMP1 = W( JJ )
               }
            } // 50

            if ( I.NE.0 ) {
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               sswap(N, Z( 1, I ), 1, Z( 1, J ), 1 );
               if ( INFO.NE.0 ) {
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               }
            }
         } // 60
      }

      // Set WORK(1) to optimal workspace size.

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

      RETURN

      // End of SSYEVX

      }
