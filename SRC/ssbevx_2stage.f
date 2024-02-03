      SUBROUTINE SSBEVX_2STAGE( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO )

      IMPLICIT NONE

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, IU, KD, LDAB, LDQ, LDZ, M, N, LWORK;
      REAL               ABSTOL, VL, VU
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      REAL               AB( LDAB, * ), Q( LDQ, * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LOWER, TEST, VALEIG, WANTZ, LQUERY;
      String             ORDER;
      int                I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL, INDISP, INDIWO, INDWRK, ISCALE, ITMP1, J, JJ, LLWORK, LWMIN, LHTRD, LWTRD, IB, INDHOUS, NSPLIT;
      REAL               ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, VLL, VUU
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      REAL               SLAMCH, SLANSB, SROUNDUP_LWORK
      // EXTERNAL LSAME, SLAMCH, SLANSB, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMV, SLACPY, SLASCL, SSCAL, SSTEBZ, SSTEIN, SSTEQR, SSTERF, SSWAP, XERBLA, SSYTRD_SB2ST
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
      LOWER = LSAME( UPLO, 'L' )
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
      } else if ( KD.LT.0 ) {
         INFO = -5
      } else if ( LDAB.LT.KD+1 ) {
         INFO = -7
      } else if ( WANTZ .AND. LDQ.LT.MAX( 1, N ) ) {
         INFO = -9
      } else {
         if ( VALEIG ) {
            IF( N.GT.0 .AND. VU.LE.VL ) INFO = -11
         } else if ( INDEIG ) {
            if ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) {
               INFO = -12
            } else if ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) {
               INFO = -13
            }
         }
      }
      if ( INFO.EQ.0 ) {
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) INFO = -18
      }

      if ( INFO.EQ.0 ) {
         if ( N.LE.1 ) {
            LWMIN = 1
            WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         } else {
            IB    = ILAENV2STAGE( 2, 'SSYTRD_SB2ST', JOBZ, N, KD, -1, -1 )             LHTRD = ILAENV2STAGE( 3, 'SSYTRD_SB2ST', JOBZ, N, KD, IB, -1 )             LWTRD = ILAENV2STAGE( 4, 'SSYTRD_SB2ST', JOBZ, N, KD, IB, -1 )
            LWMIN = 2*N + LHTRD + LWTRD
            WORK( 1 )  = SROUNDUP_LWORK(LWMIN)
         }

         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) INFO = -20
      }

      if ( INFO.NE.0 ) {
         xerbla('SSBEVX_2STAGE ', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      M = 0
      IF( N.EQ.0 ) RETURN

      if ( N.EQ.1 ) {
         M = 1
         if ( LOWER ) {
            TMP1 = AB( 1, 1 )
         } else {
            TMP1 = AB( KD+1, 1 )
         }
         if ( VALEIG ) {
            IF( .NOT.( VL.LT.TMP1 .AND. VU.GE.TMP1 ) ) M = 0
         }
         if ( M.EQ.1 ) {
            W( 1 ) = TMP1
            IF( WANTZ ) Z( 1, 1 ) = ONE
         }
         RETURN
      }

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS    = SLAMCH( 'Precision' )
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
      } else {
         VLL = ZERO
         VUU = ZERO
      }
      ANRM = SLANSB( 'M', UPLO, N, KD, AB, LDAB, WORK )
      if ( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / ANRM
      } else if ( ANRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / ANRM
      }
      if ( ISCALE.EQ.1 ) {
         if ( LOWER ) {
            slascl('B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         } else {
            slascl('Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO );
         }
         IF( ABSTOL.GT.0 ) ABSTLL = ABSTOL*SIGMA
         if ( VALEIG ) {
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         }
      }

      // Call SSYTRD_SB2ST to reduce symmetric band matrix to tridiagonal form.

      INDD    = 1
      INDE    = INDD + N
      INDHOUS = INDE + N
      INDWRK  = INDHOUS + LHTRD
      LLWORK  = LWORK - INDWRK + 1

      ssytrd_sb2st("N", JOBZ, UPLO, N, KD, AB, LDAB, WORK( INDD ), WORK( INDE ), WORK( INDHOUS ), LHTRD, WORK( INDWRK ), LLWORK, IINFO );

      // If all eigenvalues are desired and ABSTOL is less than or equal
      // to zero, then call SSTERF or SSTEQR.  If this fails for some
      // eigenvalue, then try SSTEBZ.

      TEST = .FALSE.
      if (INDEIG) {
         if (IL.EQ.1 .AND. IU.EQ.N) {
            TEST = .TRUE.
         }
      }
      if ((ALLEIG .OR. TEST) .AND. (ABSTOL.LE.ZERO)) {
         scopy(N, WORK( INDD ), 1, W, 1 );
         INDEE = INDWRK + 2*N
         if ( .NOT.WANTZ ) {
            scopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            ssterf(N, W, WORK( INDEE ), INFO );
         } else {
            slacpy('A', N, N, Q, LDQ, Z, LDZ );
            scopy(N-1, WORK( INDE ), 1, WORK( INDEE ), 1 );
            ssteqr(JOBZ, N, W, WORK( INDEE ), Z, LDZ, WORK( INDWRK ), INFO );
            if ( INFO.EQ.0 ) {
               for (I = 1; I <= N; I++) { // 10
                  IFAIL( I ) = 0
               } // 10
            }
         }
         if ( INFO.EQ.0 ) {
            M = N
            GO TO 30
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

         for (J = 1; J <= M; J++) { // 20
            scopy(N, Z( 1, J ), 1, WORK( 1 ), 1 );
            sgemv('N', N, N, ONE, Q, LDQ, WORK, 1, ZERO, Z( 1, J ), 1 );
         } // 20
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

      } // 30
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
         DO 50 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 40 JJ = J + 1, M
               if ( W( JJ ).LT.TMP1 ) {
                  I = JJ
                  TMP1 = W( JJ )
               }
            } // 40

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
         } // 50
      }

      // Set WORK(1) to optimal workspace size.

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)

      RETURN

      // End of SSBEVX_2STAGE

      }
