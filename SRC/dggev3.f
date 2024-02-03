      SUBROUTINE DGGEV3( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY;
      String             CHTEMP;
      int                ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO, IN, IRIGHT, IROWS, ITAU, IWRK, JC, JR, LWKOPT, LWKMIN;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SMLNUM, TEMP;
      // ..
      // .. Local Arrays ..
      bool               LDUMMA( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEQRF, DGGBAK, DGGBAL, DGGHD3, DLAQZ0, DLACPY, DLASCL, DLASET, DORGQR, DORMQR, DTGEVC, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, DLAMCH, DLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
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
      LWKMIN = MAX( 1, 8*N )
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

      // Compute workspace

      IF( INFO.EQ.0 ) THEN
         CALL DGEQRF( N, N, B, LDB, WORK, WORK, -1, IERR )
         LWKOPT = MAX( LWKMIN, 3*N+INT( WORK( 1 ) ) )
         CALL DORMQR( 'L', 'T', N, N, N, B, LDB, WORK, A, LDA, WORK, -1, IERR )
         LWKOPT = MAX( LWKOPT, 3*N+INT( WORK( 1 ) ) )
         IF( ILVL ) THEN
            CALL DORGQR( N, N, N, VL, LDVL, WORK, WORK, -1, IERR )
            LWKOPT = MAX( LWKOPT, 3*N+INT( WORK( 1 ) ) )
         END IF
         IF( ILV ) THEN
            CALL DGGHD3( JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK, -1, IERR )
            LWKOPT = MAX( LWKOPT, 3*N+INT( WORK ( 1 ) ) )
            CALL DLAQZ0( 'S', JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, -1, 0, IERR )
            LWKOPT = MAX( LWKOPT, 2*N+INT( WORK( 1 ) ) )
         ELSE
            CALL DGGHD3( 'N', 'N', N, 1, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK, -1, IERR )
            LWKOPT = MAX( LWKOPT, 3*N+INT( WORK( 1 ) ) )
            CALL DLAQZ0( 'E', JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, -1, 0, IERR )
            LWKOPT = MAX( LWKOPT, 2*N+INT( WORK( 1 ) ) )
         END IF
         IF( N.EQ.0 ) THEN
            WORK( 1 ) = 1
         ELSE
            WORK( 1 ) = LWKOPT
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGEV3 ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Get machine constants

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', N, N, A, LDA, WORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
      IF( ILASCL ) CALL DLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = DLANGE( 'M', N, N, B, LDB, WORK )
      ILBSCL = .FALSE.
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      ELSE IF( BNRM.GT.BIGNUM ) THEN
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      END IF
      IF( ILBSCL ) CALL DLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )

      // Permute the matrices A, B to isolate eigenvalues if possible

      ILEFT = 1
      IRIGHT = N + 1
      IWRK = IRIGHT + N
      CALL DGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWRK ), IERR )

      // Reduce B to triangular form (QR decomposition of B)

      IROWS = IHI + 1 - ILO
      IF( ILV ) THEN
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      END IF
      ITAU = IWRK
      IWRK = ITAU + IROWS
      CALL DGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )

      // Apply the orthogonal transformation to matrix A

      CALL DORMQR( 'L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR )

      // Initialize VL

      IF( ILVL ) THEN
         CALL DLASET( 'Full', N, N, ZERO, ONE, VL, LDVL )
         IF( IROWS.GT.1 ) THEN
            CALL DLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL )
         END IF
         CALL DORGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF

      // Initialize VR

      IF( ILVR ) CALL DLASET( 'Full', N, N, ZERO, ONE, VR, LDVR )

      // Reduce to generalized Hessenberg form

      IF( ILV ) THEN

         // Eigenvectors requested -- work on whole matrix.

         CALL DGGHD3( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, IERR )
      ELSE
         CALL DGGHD3( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF

      // Perform QZ algorithm (Compute eigenvalues, and optionally, the
      // Schur forms and Schur vectors)

      IWRK = ITAU
      IF( ILV ) THEN
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      END IF
      CALL DLAQZ0( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, 0, IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 110
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
         CALL DTGEVC( CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWRK ), IERR )
         IF( IERR.NE.0 ) THEN
            INFO = N + 2
            GO TO 110
         END IF

         // Undo balancing on VL and VR and normalization

         IF( ILVL ) THEN
            CALL DGGBAK( 'P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VL, LDVL, IERR )
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
               IF( TEMP.LT.SMLNUM ) GO TO 50
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
            CALL DGGBAK( 'P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VR, LDVR, IERR )
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
               IF( TEMP.LT.SMLNUM ) GO TO 100
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

         // End of eigenvector calculation

      END IF

      // Undo scaling if necessary

  110 CONTINUE

      IF( ILASCL ) THEN
         CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR )
         CALL DLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR )
      END IF

      IF( ILBSCL ) THEN
         CALL DLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
      END IF

      WORK( 1 ) = LWKOPT
      RETURN

      // End of DGGEV3

      END
