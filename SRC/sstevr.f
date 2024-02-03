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
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, TEST, LQUERY, VALEIG, WANTZ, TRYRAC;
      String             ORDER;
      int                I, IEEEOK, IMAX, INDIBL, INDIFL, INDISP, INDIWO, ISCALE, J, JJ, LIWMIN, LWMIN, NSPLIT       REAL               BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM, TMP1, TNRM, VLL, VUU;
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
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE
         IF( VALEIG ) THEN
            IF( N.GT.0 .AND. VU.LE.VL ) INFO = -7
         ELSE IF( INDEIG ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) THEN
               INFO = -8
            ELSE IF( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) THEN
               INFO = -9
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
            INFO = -14
         END IF
      END IF

      IF( INFO.EQ.0 ) THEN
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         IWORK( 1 ) = LIWMIN

         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -17
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -19
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SSTEVR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      M = 0
      IF( N.EQ.0 ) RETURN

      IF( N.EQ.1 ) THEN
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = D( 1 )
         ELSE
            IF( VL.LT.D( 1 ) .AND. VU.GE.D( 1 ) ) THEN
               M = 1
               W( 1 ) = D( 1 )
            END IF
         END IF
         IF( WANTZ ) Z( 1, 1 ) = ONE
         RETURN
      END IF

      // Get machine constants.

      SAFMIN = SLAMCH( 'Safe minimum' )
      EPS = SLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )


      // Scale matrix to allowable range, if necessary.

      ISCALE = 0
      IF( VALEIG ) THEN
         VLL = VL
         VUU = VU
      END IF

      TNRM = SLANST( 'M', N, D, E )
      IF( TNRM.GT.ZERO .AND. TNRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / TNRM
      ELSE IF( TNRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / TNRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         CALL SSCAL( N, SIGMA, D, 1 )
         CALL SSCAL( N-1, SIGMA, E( 1 ), 1 )
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF

      // Initialize indices into workspaces.  Note: These indices are used only
      // if SSTERF or SSTEMR fail.

      // IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in SSTEBZ and
      // stores the block indices of each of the M<=N eigenvalues.
      INDIBL = 1
      // IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in SSTEBZ and
      // stores the starting and finishing indices of each block.
      INDISP = INDIBL + N
      // IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
     t // hat corresponding to eigenvectors that fail to converge in
      // SSTEIN.  This information is discarded; if any fail, the driver
      // returns INFO > 0.
      INDIFL = INDISP + N
      // INDIWO is the offset of the remaining integer workspace.
      INDIWO = INDISP + N

      // If all eigenvalues are desired, then
      // call SSTERF or SSTEMR.  If this fails for some eigenvalue, then
     t // ry SSTEBZ.


      TEST = .FALSE.
      IF( INDEIG ) THEN
         IF( IL.EQ.1 .AND. IU.EQ.N ) THEN
            TEST = .TRUE.
         END IF
      END IF
      IF( ( ALLEIG .OR. TEST ) .AND. IEEEOK.EQ.1 ) THEN
         CALL SCOPY( N-1, E( 1 ), 1, WORK( 1 ), 1 )
         IF( .NOT.WANTZ ) THEN
            CALL SCOPY( N, D, 1, W, 1 )
            CALL SSTERF( N, W, WORK, INFO )
         ELSE
            CALL SCOPY( N, D, 1, WORK( N+1 ), 1 )
            IF (ABSTOL .LE. TWO*N*EPS) THEN
               TRYRAC = .TRUE.
            ELSE
               TRYRAC = .FALSE.
            END IF
            CALL SSTEMR( JOBZ, 'A', N, WORK( N+1 ), WORK, VL, VU, IL, IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC, WORK( 2*N+1 ), LWORK-2*N, IWORK, LIWORK, INFO )

         END IF
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 10
         END IF
         INFO = 0
      END IF

      // Otherwise, call SSTEBZ and, if eigenvectors are desired, SSTEIN.

      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
       CALL SSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTOL, D, E, M, NSPLIT, W, IWORK( INDIBL ), IWORK( INDISP ), WORK, IWORK( INDIWO ), INFO )

      IF( WANTZ ) THEN
         CALL SSTEIN( N, D, E, M, W, IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ, WORK, IWORK( INDIWO ), IWORK( INDIFL ), INFO )
      END IF

      // If matrix was scaled, then rescale eigenvalues appropriately.

   10 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = M
         ELSE
            IMAX = INFO - 1
         END IF
         CALL SSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF

      // If eigenvalues are not in order, then sort them, along with
      // eigenvectors.

      IF( WANTZ ) THEN
         DO 30 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 20 JJ = J + 1, M
               IF( W( JJ ).LT.TMP1 ) THEN
                  I = JJ
                  TMP1 = W( JJ )
               END IF
   20       CONTINUE

            IF( I.NE.0 ) THEN
               W( I ) = W( J )
               W( J ) = TMP1
               CALL SSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
            END IF
   30    CONTINUE
      END IF

       // Causes problems with tests 19 & 20:
       // IF (wantz .and. INDEIG ) Z( 1,1) = Z(1,1) / 1.002 + .002


      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN
      RETURN

      // End of SSTEVR

      END
