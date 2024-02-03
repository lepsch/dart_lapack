      SUBROUTINE DGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILIMIT, ILV, ILVL, ILVR, LQUERY
      CHARACTER          CHTEMP
      INTEGER            ICOLS, IHI, IINFO, IJOBVL, IJOBVR, ILEFT, ILO, IN, IRIGHT, IROWS, ITAU, IWORK, JC, JR, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3       DOUBLE PRECISION   ABSAI, ABSAR, ABSB, ANRM, ANRM1, ANRM2, BNRM, BNRM1, BNRM2, EPS, ONEPLS, SAFMAX, SAFMIN, SALFAI, SALFAR, SBETA, SCALE, TEMP
*     ..
*     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQRF, DGGBAK, DGGBAL, DGGHRD, DHGEQZ, DLACPY, DLASCL, DLASET, DORGQR, DORMQR, DTGEVC, XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, MAX
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( LSAME( JOBVL, 'N' ) ) THEN
         IJOBVL = 1
         ILVL = .FALSE.
      ELSE IF( LSAME( JOBVL, 'V' ) ) THEN
         IJOBVL = 2
         ILVL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVL = .FALSE.
      END IF
*
      IF( LSAME( JOBVR, 'N' ) ) THEN
         IJOBVR = 1
         ILVR = .FALSE.
      ELSE IF( LSAME( JOBVR, 'V' ) ) THEN
         IJOBVR = 2
         ILVR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVR = .FALSE.
      END IF
      ILV = ILVL .OR. ILVR
*
*     Test the input arguments
*
      LWKMIN = MAX( 8*N, 1 )
      LWKOPT = LWKMIN
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      INFO = 0
      IF( IJOBVL.LE.0 ) THEN
         INFO = -1
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( ILVL .AND. LDVL.LT.N ) ) THEN
         INFO = -12
      ELSE IF( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -16
      END IF
*
      IF( INFO.EQ.0 ) THEN
         NB1 = ILAENV( 1, 'DGEQRF', ' ', N, N, -1, -1 )
         NB2 = ILAENV( 1, 'DORMQR', ' ', N, N, N, -1 )
         NB3 = ILAENV( 1, 'DORGQR', ' ', N, N, N, -1 )
         NB = MAX( NB1, NB2, NB3 )
         LOPT = 2*N + MAX( 6*N, N*( NB+1 ) )
         WORK( 1 ) = LOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEGV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'E' )*DLAMCH( 'B' )
      SAFMIN = DLAMCH( 'S' )
      SAFMIN = SAFMIN + SAFMIN
      SAFMAX = ONE / SAFMIN
      ONEPLS = ONE + ( 4*EPS )
*
*     Scale A
*
      ANRM = DLANGE( 'M', N, N, A, LDA, WORK )
      ANRM1 = ANRM
      ANRM2 = ONE
      IF( ANRM.LT.ONE ) THEN
         IF( SAFMAX*ANRM.LT.ONE ) THEN
            ANRM1 = SAFMIN
            ANRM2 = SAFMAX*ANRM
         END IF
      END IF
*
      IF( ANRM.GT.ZERO ) THEN
         CALL DLASCL( 'G', -1, -1, ANRM, ONE, N, N, A, LDA, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 10
            RETURN
         END IF
      END IF
*
*     Scale B
*
      BNRM = DLANGE( 'M', N, N, B, LDB, WORK )
      BNRM1 = BNRM
      BNRM2 = ONE
      IF( BNRM.LT.ONE ) THEN
         IF( SAFMAX*BNRM.LT.ONE ) THEN
            BNRM1 = SAFMIN
            BNRM2 = SAFMAX*BNRM
         END IF
      END IF
*
      IF( BNRM.GT.ZERO ) THEN
         CALL DLASCL( 'G', -1, -1, BNRM, ONE, N, N, B, LDB, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 10
            RETURN
         END IF
      END IF
*
*     Permute the matrix to make it more nearly triangular
*     Workspace layout:  (8*N words -- "work" requires 6*N words)
*        left_permutation, right_permutation, work...
*
      ILEFT = 1
      IRIGHT = N + 1
      IWORK = IRIGHT + N
      CALL DGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWORK ), IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 1
         GO TO 120
      END IF
*
*     Reduce B to triangular form, and initialize VL and/or VR
*     Workspace layout:  ("work..." must have at least N words)
*        left_permutation, right_permutation, tau, work...
*
      IROWS = IHI + 1 - ILO
      IF( ILV ) THEN
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      END IF
      ITAU = IWORK
      IWORK = ITAU + IROWS
      CALL DGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )       IF( IINFO.GE.0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 2
         GO TO 120
      END IF
*
      CALL DORMQR( 'L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWORK ), LWORK+1-IWORK, IINFO )
      IF( IINFO.GE.0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      IF( IINFO.NE.0 ) THEN
         INFO = N + 3
         GO TO 120
      END IF
*
      IF( ILVL ) THEN
         CALL DLASET( 'Full', N, N, ZERO, ONE, VL, LDVL )
         CALL DLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL )          CALL DORGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWORK ), LWORK+1-IWORK, IINFO )
         IF( IINFO.GE.0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 4
            GO TO 120
         END IF
      END IF
*
      IF( ILVR ) CALL DLASET( 'Full', N, N, ZERO, ONE, VR, LDVR )
*
*     Reduce to generalized Hessenberg form
*
      IF( ILV ) THEN
*
*        Eigenvectors requested -- work on whole matrix.
*
         CALL DGGHRD( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, IINFO )
      ELSE
         CALL DGGHRD( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IINFO )
      END IF
      IF( IINFO.NE.0 ) THEN
         INFO = N + 5
         GO TO 120
      END IF
*
*     Perform QZ algorithm
*     Workspace layout:  ("work..." must have at least 1 word)
*        left_permutation, right_permutation, work...
*
      IWORK = ITAU
      IF( ILV ) THEN
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      END IF
      CALL DHGEQZ( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK( IWORK ), LWORK+1-IWORK, IINFO )
      IF( IINFO.GE.0 ) LWKOPT = MAX( LWKOPT, INT( WORK( IWORK ) )+IWORK-1 )
      IF( IINFO.NE.0 ) THEN
         IF( IINFO.GT.0 .AND. IINFO.LE.N ) THEN
            INFO = IINFO
         ELSE IF( IINFO.GT.N .AND. IINFO.LE.2*N ) THEN
            INFO = IINFO - N
         ELSE
            INFO = N + 6
         END IF
         GO TO 120
      END IF
*
      IF( ILV ) THEN
*
*        Compute Eigenvectors  (DTGEVC requires 6*N words of workspace)
*
         IF( ILVL ) THEN
            IF( ILVR ) THEN
               CHTEMP = 'B'
            ELSE
               CHTEMP = 'L'
            END IF
         ELSE
            CHTEMP = 'R'
         END IF
*
         CALL DTGEVC( CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWORK ), IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = N + 7
            GO TO 120
         END IF
*
*        Undo balancing on VL and VR, rescale
*
         IF( ILVL ) THEN
            CALL DGGBAK( 'P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VL, LDVL, IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = N + 8
               GO TO 120
            END IF
            DO 50 JC = 1, N
               IF( ALPHAI( JC ).LT.ZERO ) GO TO 50
               TEMP = ZERO
               IF( ALPHAI( JC ).EQ.ZERO ) THEN
                  DO 10 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) ) )
   10             CONTINUE
               ELSE
                  DO 20 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) )+ ABS( VL( JR, JC+1 ) ) )
   20             CONTINUE
               END IF
               IF( TEMP.LT.SAFMIN ) GO TO 50
               TEMP = ONE / TEMP
               IF( ALPHAI( JC ).EQ.ZERO ) THEN
                  DO 30 JR = 1, N
                     VL( JR, JC ) = VL( JR, JC )*TEMP
   30             CONTINUE
               ELSE
                  DO 40 JR = 1, N
                     VL( JR, JC ) = VL( JR, JC )*TEMP
                     VL( JR, JC+1 ) = VL( JR, JC+1 )*TEMP
   40             CONTINUE
               END IF
   50       CONTINUE
         END IF
         IF( ILVR ) THEN
            CALL DGGBAK( 'P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VR, LDVR, IINFO )
            IF( IINFO.NE.0 ) THEN
               INFO = N + 9
               GO TO 120
            END IF
            DO 100 JC = 1, N
               IF( ALPHAI( JC ).LT.ZERO ) GO TO 100
               TEMP = ZERO
               IF( ALPHAI( JC ).EQ.ZERO ) THEN
                  DO 60 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) ) )
   60             CONTINUE
               ELSE
                  DO 70 JR = 1, N
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) )+ ABS( VR( JR, JC+1 ) ) )
   70             CONTINUE
               END IF
               IF( TEMP.LT.SAFMIN ) GO TO 100
               TEMP = ONE / TEMP
               IF( ALPHAI( JC ).EQ.ZERO ) THEN
                  DO 80 JR = 1, N
                     VR( JR, JC ) = VR( JR, JC )*TEMP
   80             CONTINUE
               ELSE
                  DO 90 JR = 1, N
                     VR( JR, JC ) = VR( JR, JC )*TEMP
                     VR( JR, JC+1 ) = VR( JR, JC+1 )*TEMP
   90             CONTINUE
               END IF
  100       CONTINUE
         END IF
*
*        End of eigenvector calculation
*
      END IF
*
*     Undo scaling in alpha, beta
*
*     Note: this does not give the alpha and beta for the unscaled
*     problem.
*
*     Un-scaling is limited to avoid underflow in alpha and beta
*     if they are significant.
*
      DO 110 JC = 1, N
         ABSAR = ABS( ALPHAR( JC ) )
         ABSAI = ABS( ALPHAI( JC ) )
         ABSB = ABS( BETA( JC ) )
         SALFAR = ANRM*ALPHAR( JC )
         SALFAI = ANRM*ALPHAI( JC )
         SBETA = BNRM*BETA( JC )
         ILIMIT = .FALSE.
         SCALE = ONE
*
*        Check for significant underflow in ALPHAI
*
         IF( ABS( SALFAI ).LT.SAFMIN .AND. ABSAI.GE. MAX( SAFMIN, EPS*ABSAR, EPS*ABSB ) ) THEN
            ILIMIT = .TRUE.
            SCALE = ( ONEPLS*SAFMIN / ANRM1 ) / MAX( ONEPLS*SAFMIN, ANRM2*ABSAI )
*
         ELSE IF( SALFAI.EQ.ZERO ) THEN
*
*           If insignificant underflow in ALPHAI, then make the
*           conjugate eigenvalue real.
*
            IF( ALPHAI( JC ).LT.ZERO .AND. JC.GT.1 ) THEN
               ALPHAI( JC-1 ) = ZERO
            ELSE IF( ALPHAI( JC ).GT.ZERO .AND. JC.LT.N ) THEN
               ALPHAI( JC+1 ) = ZERO
            END IF
         END IF
*
*        Check for significant underflow in ALPHAR
*
         IF( ABS( SALFAR ).LT.SAFMIN .AND. ABSAR.GE. MAX( SAFMIN, EPS*ABSAI, EPS*ABSB ) ) THEN
            ILIMIT = .TRUE.
            SCALE = MAX( SCALE, ( ONEPLS*SAFMIN / ANRM1 ) / MAX( ONEPLS*SAFMIN, ANRM2*ABSAR ) )
         END IF
*
*        Check for significant underflow in BETA
*
         IF( ABS( SBETA ).LT.SAFMIN .AND. ABSB.GE. MAX( SAFMIN, EPS*ABSAR, EPS*ABSAI ) ) THEN
            ILIMIT = .TRUE.
            SCALE = MAX( SCALE, ( ONEPLS*SAFMIN / BNRM1 ) / MAX( ONEPLS*SAFMIN, BNRM2*ABSB ) )
         END IF
*
*        Check for possible overflow when limiting scaling
*
         IF( ILIMIT ) THEN
            TEMP = ( SCALE*SAFMIN )*MAX( ABS( SALFAR ), ABS( SALFAI ), ABS( SBETA ) )             IF( TEMP.GT.ONE ) SCALE = SCALE / TEMP             IF( SCALE.LT.ONE ) ILIMIT = .FALSE.
         END IF
*
*        Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary.
*
         IF( ILIMIT ) THEN
            SALFAR = ( SCALE*ALPHAR( JC ) )*ANRM
            SALFAI = ( SCALE*ALPHAI( JC ) )*ANRM
            SBETA = ( SCALE*BETA( JC ) )*BNRM
         END IF
         ALPHAR( JC ) = SALFAR
         ALPHAI( JC ) = SALFAI
         BETA( JC ) = SBETA
  110 CONTINUE
*
  120 CONTINUE
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of DGEGV
*
      END
