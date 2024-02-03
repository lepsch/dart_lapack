      SUBROUTINE CGGEV3( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E0, 0.0E0 ), CONE = ( 1.0E0, 0.0E0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY;
      String             CHTEMP;
      int                ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO, IN, IRIGHT, IROWS, IRWRK, ITAU, IWRK, JC, JR, LWKOPT, LWKMIN;
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SMLNUM, TEMP
      COMPLEX            X
      // ..
      // .. Local Arrays ..
      bool               LDUMMA( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQRF, CGGBAK, CGGBAL, CGGHD3, CLAQZ0, CLACPY, CLASCL, CLASET, CTGEVC, CUNGQR, CUNMQR, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, SLAMCH, SROUNDUP_LWORK
      // EXTERNAL LSAME, CLANGE, SLAMCH, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, REAL, SQRT
      // ..
      // .. Statement Functions ..
      REAL               ABS1
      // ..
      // .. Statement Function definitions ..
      ABS1( X ) = ABS( REAL( X ) ) + ABS( AIMAG( X ) )
      // ..
      // .. Executable Statements ..

      // Decode the input arguments

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

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      LWKMIN = MAX( 1, 2*N )
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
         INFO = -11
      ELSE IF( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) THEN
         INFO = -13
      ELSE IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -15
      END IF

      // Compute workspace

      IF( INFO.EQ.0 ) THEN
         CALL CGEQRF( N, N, B, LDB, WORK, WORK, -1, IERR )
         LWKOPT = MAX( LWKMIN, N+INT( WORK( 1 ) ) )
         CALL CUNMQR( 'L', 'C', N, N, N, B, LDB, WORK, A, LDA, WORK, -1, IERR )
         LWKOPT = MAX( LWKOPT, N+INT( WORK( 1 ) ) )
         IF( ILVL ) THEN
            CALL CUNGQR( N, N, N, VL, LDVL, WORK, WORK, -1, IERR )
            LWKOPT = MAX( LWKOPT, N+INT( WORK( 1 ) ) )
         END IF
         IF( ILV ) THEN
            CALL CGGHD3( JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK, -1, IERR )
            LWKOPT = MAX( LWKOPT, N+INT( WORK( 1 ) ) )
            CALL CLAQZ0( 'S', JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK, -1, RWORK, 0, IERR )
            LWKOPT = MAX( LWKOPT, N+INT( WORK( 1 ) ) )
         ELSE
            CALL CGGHD3( 'N', 'N', N, 1, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK, -1, IERR )
            LWKOPT = MAX( LWKOPT, N+INT( WORK( 1 ) ) )
            CALL CLAQZ0( 'E', JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK, -1, RWORK, 0, IERR )
            LWKOPT = MAX( LWKOPT, N+INT( WORK( 1 ) ) )
         END IF
         IF( N.EQ.0 ) THEN
            WORK( 1 ) = 1
         ELSE
            WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGGEV3 ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Get machine constants

      EPS = SLAMCH( 'E' )*SLAMCH( 'B' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
      IF( ILASCL ) CALL CLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = CLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      ELSE IF( BNRM.GT.BIGNUM ) THEN
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      END IF
      IF( ILBSCL ) CALL CLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )

      // Permute the matrices A, B to isolate eigenvalues if possible

      ILEFT = 1
      IRIGHT = N + 1
      IRWRK = IRIGHT + N
      CALL CGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), RWORK( IRWRK ), IERR )

      // Reduce B to triangular form (QR decomposition of B)

      IROWS = IHI + 1 - ILO
      IF( ILV ) THEN
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      END IF
      ITAU = 1
      IWRK = ITAU + IROWS
      CALL CGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )

      // Apply the orthogonal transformation to matrix A

      CALL CUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR )

      // Initialize VL

      IF( ILVL ) THEN
         CALL CLASET( 'Full', N, N, CZERO, CONE, VL, LDVL )
         IF( IROWS.GT.1 ) THEN
            CALL CLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL )
         END IF
         CALL CUNGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF

      // Initialize VR

      IF( ILVR ) CALL CLASET( 'Full', N, N, CZERO, CONE, VR, LDVR )

      // Reduce to generalized Hessenberg form

      IF( ILV ) THEN

         // Eigenvectors requested -- work on whole matrix.

         CALL CGGHD3( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, IERR )
      ELSE
         CALL CGGHD3( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF

      // Perform QZ algorithm (Compute eigenvalues, and optionally, the
      // Schur form and Schur vectors)

      IWRK = ITAU
      IF( ILV ) THEN
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      END IF
      CALL CLAQZ0( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, RWORK( IRWRK ), 0, IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 70
      END IF

      // Compute Eigenvectors

      IF( ILV ) THEN
         IF( ILVL ) THEN
            IF( ILVR ) THEN
               CHTEMP = 'B'
            ELSE
               CHTEMP = 'L'
            END IF
         ELSE
            CHTEMP = 'R'
         END IF

         CALL CTGEVC( CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWRK ), RWORK( IRWRK ), IERR )
         IF( IERR.NE.0 ) THEN
            INFO = N + 2
            GO TO 70
         END IF

         // Undo balancing on VL and VR and normalization

         IF( ILVL ) THEN
            CALL CGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VL, LDVL, IERR )
            DO 30 JC = 1, N
               TEMP = ZERO
               DO 10 JR = 1, N
                  TEMP = MAX( TEMP, ABS1( VL( JR, JC ) ) )
   10          CONTINUE
               IF( TEMP.LT.SMLNUM ) GO TO 30
               TEMP = ONE / TEMP
               DO 20 JR = 1, N
                  VL( JR, JC ) = VL( JR, JC )*TEMP
   20          CONTINUE
   30       CONTINUE
         END IF
         IF( ILVR ) THEN
            CALL CGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VR, LDVR, IERR )
            DO 60 JC = 1, N
               TEMP = ZERO
               DO 40 JR = 1, N
                  TEMP = MAX( TEMP, ABS1( VR( JR, JC ) ) )
   40          CONTINUE
               IF( TEMP.LT.SMLNUM ) GO TO 60
               TEMP = ONE / TEMP
               DO 50 JR = 1, N
                  VR( JR, JC ) = VR( JR, JC )*TEMP
   50          CONTINUE
   60       CONTINUE
         END IF
      END IF

      // Undo scaling if necessary

   70 CONTINUE

      IF( ILASCL ) CALL CLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )

      IF( ILBSCL ) CALL CLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      RETURN

      // End of CGGEV3

      }
