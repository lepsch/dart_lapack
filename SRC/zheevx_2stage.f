      SUBROUTINE ZHEEVX_2STAGE( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, LDA, LDZ, LWORK, M, N;
      double             ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      double             RWORK( * ), W( * );
      COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LOWER, LQUERY, TEST, VALEIG, WANTZ;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, INDISP, INDIWK, INDRWK, INDTAU, INDWRK, ISCALE, ITMP1, J, JJ, LLWORK, NSPLIT, LWMIN, LHTRD, LWTRD, KD, IB, INDHOUS;
      double             ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      double             DLAMCH, ZLANHE;
      // EXTERNAL LSAME, DLAMCH, ZLANHE, ILAENV2STAGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DSCAL, DSTEBZ, DSTERF, XERBLA, ZDSCAL, ZLACPY, ZSTEIN, ZSTEQR, ZSWAP, ZUNGTR, ZUNMTR, ZHETRD_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, SQRT
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
         if ( N.LE.1 ) {
            LWMIN = 1
            WORK( 1 ) = LWMIN
         } else {
            KD    = ILAENV2STAGE( 1, 'ZHETRD_2STAGE', JOBZ, N, -1, -1, -1 )             IB    = ILAENV2STAGE( 2, 'ZHETRD_2STAGE', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'ZHETRD_2STAGE', JOBZ, N, KD, IB, -1 )
            LWMIN = N + LHTRD + LWTRD
            WORK( 1 )  = LWMIN
         }

         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) INFO = -17
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZHEEVX_2STAGE', -INFO )
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
         W( 1 ) = DBLE( A( 1, 1 ) )
         } else if ( VALEIG ) {
            if ( VL.LT.DBLE( A( 1, 1 ) ) .AND. VU.GE.DBLE( A( 1, 1 ) ) ) {
               M = 1
               W( 1 ) = DBLE( A( 1, 1 ) )
            }
         }
         IF( WANTZ ) Z( 1, 1 ) = CONE
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
      if ( VALEIG ) {
         VLL = VL
         VUU = VU
      }
      ANRM = ZLANHE( 'M', UPLO, N, A, LDA, RWORK )
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
               CALL ZDSCAL( N-J+1, SIGMA, A( J, J ), 1 )
   10       CONTINUE
         } else {
            DO 20 J = 1, N
               CALL ZDSCAL( J, SIGMA, A( 1, J ), 1 )
   20       CONTINUE
         }
         IF( ABSTOL.GT.0 ) ABSTLL = ABSTOL*SIGMA
         if ( VALEIG ) {
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         }
      }

      // Call ZHETRD_2STAGE to reduce Hermitian matrix to tridiagonal form.

      INDD    = 1
      INDE    = INDD + N
      INDRWK  = INDE + N
      INDTAU  = 1
      INDHOUS = INDTAU + N
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1

      CALL ZHETRD_2STAGE( JOBZ, UPLO, N, A, LDA, RWORK( INDD ), RWORK( INDE ), WORK( INDTAU ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO )

      // If all eigenvalues are desired and ABSTOL is less than or equal to
      // zero, then call DSTERF or ZUNGTR and ZSTEQR.  If this fails for
      // some eigenvalue, then try DSTEBZ.

      TEST = .FALSE.
      if ( INDEIG ) {
         if ( IL.EQ.1 .AND. IU.EQ.N ) {
            TEST = .TRUE.
         }
      }
      if ( ( ALLEIG .OR. TEST ) .AND. ( ABSTOL.LE.ZERO ) ) {
         CALL DCOPY( N, RWORK( INDD ), 1, W, 1 )
         INDEE = INDRWK + 2*N
         if ( .NOT.WANTZ ) {
            CALL DCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
            CALL DSTERF( N, W, RWORK( INDEE ), INFO )
         } else {
            CALL ZLACPY( 'A', N, N, A, LDA, Z, LDZ )
            CALL ZUNGTR( UPLO, N, Z, LDZ, WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO )
            CALL DCOPY( N-1, RWORK( INDE ), 1, RWORK( INDEE ), 1 )
            CALL ZSTEQR( JOBZ, N, W, RWORK( INDEE ), Z, LDZ, RWORK( INDRWK ), INFO )
            if ( INFO.EQ.0 ) {
               DO 30 I = 1, N
                  IFAIL( I ) = 0
   30          CONTINUE
            }
         }
         if ( INFO.EQ.0 ) {
            M = N
            GO TO 40
         }
         INFO = 0
      }

      // Otherwise, call DSTEBZ and, if eigenvectors are desired, ZSTEIN.

      if ( WANTZ ) {
         ORDER = 'B'
      } else {
         ORDER = 'E'
      }
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWK = INDISP + N
      CALL DSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL, RWORK( INDD ), RWORK( INDE ), M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), RWORK( INDRWK ), IWORK( INDIWK ), INFO )

      if ( WANTZ ) {
         CALL ZSTEIN( N, RWORK( INDD ), RWORK( INDE ), M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, RWORK( INDRWK ), IWORK( INDIWK ), IFAIL, INFO )

         // Apply unitary matrix used in reduction to tridiagonal
         // form to eigenvectors returned by ZSTEIN.

         CALL ZUNMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z, LDZ, WORK( INDWRK ), LLWORK, IINFO )
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

   40 CONTINUE
      if ( ISCALE.EQ.1 ) {
         if ( INFO.EQ.0 ) {
            IMAX = M
         } else {
            IMAX = INFO - 1
         }
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      }

      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.

      if ( WANTZ ) {
         DO 60 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 50 JJ = J + 1, M
               if ( W( JJ ).LT.TMP1 ) {
                  I = JJ
                  TMP1 = W( JJ )
               }
   50       CONTINUE

            if ( I.NE.0 ) {
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
               if ( INFO.NE.0 ) {
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               }
            }
   60    CONTINUE
      }

      // Set WORK(1) to optimal complex workspace size.

      WORK( 1 ) = LWMIN

      RETURN

      // End of ZHEEVX_2STAGE

      }
