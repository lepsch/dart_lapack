      SUBROUTINE SSTEVR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE;
      int                IL, INFO, IU, LDZ, LIWORK, LWORK, M, N;
      REAL               ABSTOL, VL, VU
      // ..
      // .. Array Arguments ..
      int                ISUPPZ( * ), IWORK( * );
      REAL               D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, TEST, LQUERY, VALEIG, WANTZ, TRYRAC;
      String             ORDER;
      int                I, IEEEOK, IMAX, INDIBL, INDIFL, INDISP, INDIWO, ISCALE, J, JJ, LIWMIN, LWMIN, NSPLIT;
      REAL               BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, TNRM, VLL, VUU;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANST, SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SLAMCH, SLANST, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SSCAL, SSTEBZ, SSTEMR, SSTEIN, SSTERF, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..


      // Test the input parameters.

      IEEEOK = ILAENV( 10, 'SSTEVR', 'N', 1, 2, 3, 4 )

      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )

      LQUERY = ( ( LWORK.EQ.-1 ) .OR. ( LIWORK.EQ.-1 ) )
      LWMIN = MAX( 1, 20*N )
      LIWMIN = MAX(1, 10*N )


      INFO = 0
      if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else {
         if ( VALEIG ) {
            IF( N.GT.0 .AND. VU.LE.VL ) INFO = -7
         } else if ( INDEIG ) {
            if ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) {
               INFO = -8
            } else if ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) {
               INFO = -9
            }
         }
      }
      if ( INFO.EQ.0 ) {
         if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
            INFO = -14
         }
      }

      if ( INFO.EQ.0 ) {
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -17
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
            INFO = -19
         }
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SSTEVR', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      M = 0
      IF( N.EQ.0 ) RETURN

      if ( N.EQ.1 ) {
         if ( ALLEIG .OR. INDEIG ) {
            M = 1
            W( 1 ) = D( 1 )
         } else {
            if ( VL.LT.D( 1 ) .AND. VU.GE.D( 1 ) ) {
               M = 1
               W( 1 ) = D( 1 )
            }
         }
         IF( WANTZ ) Z( 1, 1 ) = ONE
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
      if ( VALEIG ) {
         VLL = VL
         VUU = VU
      }

      TNRM = SLANST( 'M', N, D, E )
      if ( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) {
         ISCALE = 1
         SIGMA = RMIN / TNRM
      } else if ( TNRM.GT.RMAX ) {
         ISCALE = 1
         SIGMA = RMAX / TNRM
      }
      if ( ISCALE.EQ.1 ) {
         CALL SSCAL( N, SIGMA, D, 1 )
         CALL SSCAL( N-1, SIGMA, E( 1 ), 1 )
         if ( VALEIG ) {
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         }
      }

      // Initialize indices into workspaces.  Note: These indices are used only
      // if SSTERF or SSTEMR fail.

      // IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in SSTEBZ and
      // stores the block indices of each of the M<=N eigenvalues.
      INDIBL = 1
      // IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in SSTEBZ and
      // stores the starting and finishing indices of each block.
      INDISP = INDIBL + N
      // IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
      // that corresponding to eigenvectors that fail to converge in
      // SSTEIN.  This information is discarded; if any fail, the driver
      // returns INFO > 0.
      INDIFL = INDISP + N
      // INDIWO is the offset of the remaining integer workspace.
      INDIWO = INDISP + N

      // If all eigenvalues are desired, then
      // call SSTERF or SSTEMR.  If this fails for some eigenvalue, then
      // try SSTEBZ.


      TEST = .FALSE.
      if ( INDEIG ) {
         if ( IL.EQ.1 .AND. IU.EQ.N ) {
            TEST = .TRUE.
         }
      }
      if ( ( ALLEIG .OR. TEST ) .AND. IEEEOK.EQ.1 ) {
         CALL SCOPY( N-1, E( 1 ), 1, WORK( 1 ), 1 )
         if ( .NOT.WANTZ ) {
            CALL SCOPY( N, D, 1, W, 1 )
            CALL SSTERF( N, W, WORK, INFO )
         } else {
            CALL SCOPY( N, D, 1, WORK( N+1 ), 1 )
            if (ABSTOL .LE. TWO*N*EPS) {
               TRYRAC = .TRUE.
            } else {
               TRYRAC = .FALSE.
            }
            CALL SSTEMR( JOBZ, 'A', N, WORK( N+1 ), WORK, VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK( 2*N+1 ), LWORK-2*N, IWORK, LIWORK, INFO )

         }
         if ( INFO.EQ.0 ) {
            M = N
            GO TO 10
         }
         INFO = 0
      }

      // Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN.

      if ( WANTZ ) {
         ORDER = 'B'
      } else {
         ORDER = 'E'
      }
       CALL SSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTOL, D, E, M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), WORK, IWORK( INDIWO ), INFO )

      if ( WANTZ ) {
         CALL SSTEIN( N, D, E, M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, WORK, IWORK( INDIWO ), IWORK( INDIFL ), INFO )
      }

      // If matrix was scaled, then rescale eigenvalues appropriately.

   10 CONTINUE
      if ( ISCALE.EQ.1 ) {
         if ( INFO.EQ.0 ) {
            IMAX = M
         } else {
            IMAX = INFO - 1
         }
         CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
      }

      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.

      if ( WANTZ ) {
         DO 30 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 20 JJ = J + 1, M
               if ( W( JJ ).LT.TMP1 ) {
                  I = JJ
                  TMP1 = W( JJ )
               }
   20       CONTINUE

            if ( I.NE.0 ) {
               W( I ) = W( J )
               W( J ) = TMP1
               CALL SSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
            }
   30    CONTINUE
      }

       // Causes problems with tests 19 & 20:
       // IF (wantz .and. INDEIG ) Z( 1,1) = Z(1,1) / 1.002 + .002


      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN
      RETURN

      // End of SSTEVR

      }
