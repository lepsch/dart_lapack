      SUBROUTINE CGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
      implicit none

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDVL, LDVR, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL   RWORK( * )
      COMPLEX         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL   ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, SCALEA, WANTVL, WANTVR;
      String             SIDE;
      int                HSWORK, I, IBAL, IERR, IHI, ILO, IRWORK, ITAU, IWRK, K, LWORK_TREVC, MAXWRK, MINWRK, NOUT;
      REAL   ANRM, BIGNUM, CSCALE, EPS, SCL, SMLNUM
      COMPLEX         TMP
      // ..
      // .. Local Arrays ..
      bool               SELECT( 1 );
      REAL   DUM( 1 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CSSCAL, CGEBAK, CGEBAL, CGEHRD, CHSEQR, CLACPY, CLASCL, CSCAL, CTREVC3, CUNGHR
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX, ILAENV;
      REAL               SLAMCH, SCNRM2, CLANGE, SROUNDUP_LWORK
      // EXTERNAL LSAME, ISAMAX, ILAENV, SLAMCH, SCNRM2, CLANGE, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL, CMPLX, CONJG, AIMAG, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVL = LSAME( JOBVL, 'V' )
      WANTVR = LSAME( JOBVR, 'V' )
      IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      END IF

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // CWorkspace refers to complex workspace, and RWorkspace to real
        // workspace. NB refers to the optimal block size for the
        // immediately following subroutine, as returned by ILAENV.
        // HSWORK refers to the workspace preferred by CHSEQR, as
        // calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
       t // he worst case.)

      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         } else {
            MAXWRK = N + N*ILAENV( 1, 'CGEHRD', ' ', N, 1, N, 0 )
            MINWRK = 2*N
            IF( WANTVL ) THEN
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'CUNGHR', ' ', N, 1, N, -1 ) )                CALL CTREVC3( 'L', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, RWORK, -1, IERR )
               LWORK_TREVC = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
               CALL CHSEQR( 'S', 'V', N, 1, N, A, LDA, W, VL, LDVL, WORK, -1, INFO )
            ELSE IF( WANTVR ) THEN
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'CUNGHR', ' ', N, 1, N, -1 ) )                CALL CTREVC3( 'R', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, RWORK, -1, IERR )
               LWORK_TREVC = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
               CALL CHSEQR( 'S', 'V', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO )
            } else {
               CALL CHSEQR( 'E', 'N', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO )
            END IF
            HSWORK = INT( WORK(1) )
            MAXWRK = MAX( MAXWRK, HSWORK, MINWRK )
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)

         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGEEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA ) CALL CLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )

      // Balance the matrix
      // (CWorkspace: none)
      // (RWorkspace: need N)

      IBAL = 1
      CALL CGEBAL( 'B', N, A, LDA, ILO, IHI, RWORK( IBAL ), IERR )

      // Reduce to upper Hessenberg form
      // (CWorkspace: need 2*N, prefer N+N*NB)
      // (RWorkspace: none)

      ITAU = 1
      IWRK = ITAU + N
      CALL CGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )

      IF( WANTVL ) THEN

         // Want left eigenvectors
         // Copy Householder vectors to VL

         SIDE = 'L'
         CALL CLACPY( 'L', N, N, A, LDA, VL, LDVL )

         // Generate unitary matrix in VL
         // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
         // (RWorkspace: none)

         CALL CUNGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )

         // Perform QR iteration, accumulating Schur vectors in VL
         // (CWorkspace: need 1, prefer HSWORK (see comments) )
         // (RWorkspace: none)

         IWRK = ITAU
         CALL CHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, W, VL, LDVL, WORK( IWRK ), LWORK-IWRK+1, INFO )

         IF( WANTVR ) THEN

            // Want left and right eigenvectors
            // Copy Schur vectors to VR

            SIDE = 'B'
            CALL CLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         END IF

      ELSE IF( WANTVR ) THEN

         // Want right eigenvectors
         // Copy Householder vectors to VR

         SIDE = 'R'
         CALL CLACPY( 'L', N, N, A, LDA, VR, LDVR )

         // Generate unitary matrix in VR
         // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
         // (RWorkspace: none)

         CALL CUNGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )

         // Perform QR iteration, accumulating Schur vectors in VR
         // (CWorkspace: need 1, prefer HSWORK (see comments) )
         // (RWorkspace: none)

         IWRK = ITAU
         CALL CHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, W, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO )

      } else {

         // Compute eigenvalues only
         // (CWorkspace: need 1, prefer HSWORK (see comments) )
         // (RWorkspace: none)

         IWRK = ITAU
         CALL CHSEQR( 'E', 'N', N, ILO, IHI, A, LDA, W, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO )
      END IF

      // If INFO .NE. 0 from CHSEQR, then quit

      IF( INFO.NE.0 ) GO TO 50

      IF( WANTVL .OR. WANTVR ) THEN

         // Compute left and/or right eigenvectors
         // (CWorkspace: need 2*N, prefer N + 2*N*NB)
         // (RWorkspace: need 2*N)

         IRWORK = IBAL + N
         CALL CTREVC3( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK( IWRK ), LWORK-IWRK+1, RWORK( IRWORK ), N, IERR )
      END IF

      IF( WANTVL ) THEN

         // Undo balancing of left eigenvectors
         // (CWorkspace: none)
         // (RWorkspace: need N)

         CALL CGEBAK( 'B', 'L', N, ILO, IHI, RWORK( IBAL ), N, VL, LDVL, IERR )

         // Normalize left eigenvectors and make largest component real

         DO 20 I = 1, N
            SCL = ONE / SCNRM2( N, VL( 1, I ), 1 )
            CALL CSSCAL( N, SCL, VL( 1, I ), 1 )
            DO 10 K = 1, N
               RWORK( IRWORK+K-1 ) = REAL( VL( K, I ) )**2 + AIMAG( VL( K, I ) )**2
   10       CONTINUE
            K = ISAMAX( N, RWORK( IRWORK ), 1 )
            TMP = CONJG( VL( K, I ) ) / SQRT( RWORK( IRWORK+K-1 ) )
            CALL CSCAL( N, TMP, VL( 1, I ), 1 )
            VL( K, I ) = CMPLX( REAL( VL( K, I ) ), ZERO )
   20    CONTINUE
      END IF

      IF( WANTVR ) THEN

         // Undo balancing of right eigenvectors
         // (CWorkspace: none)
         // (RWorkspace: need N)

         CALL CGEBAK( 'B', 'R', N, ILO, IHI, RWORK( IBAL ), N, VR, LDVR, IERR )

         // Normalize right eigenvectors and make largest component real

         DO 40 I = 1, N
            SCL = ONE / SCNRM2( N, VR( 1, I ), 1 )
            CALL CSSCAL( N, SCL, VR( 1, I ), 1 )
            DO 30 K = 1, N
               RWORK( IRWORK+K-1 ) = REAL( VR( K, I ) )**2 + AIMAG( VR( K, I ) )**2
   30       CONTINUE
            K = ISAMAX( N, RWORK( IRWORK ), 1 )
            TMP = CONJG( VR( K, I ) ) / SQRT( RWORK( IRWORK+K-1 ) )
            CALL CSCAL( N, TMP, VR( 1, I ), 1 )
            VR( K, I ) = CMPLX( REAL( VR( K, I ) ), ZERO )
   40    CONTINUE
      END IF

      // Undo scaling if necessary

   50 CONTINUE
      IF( SCALEA ) THEN
         CALL CLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, W( INFO+1 ), MAX( N-INFO, 1 ), IERR )
         IF( INFO.GT.0 ) THEN
            CALL CLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, W, N, IERR )
         END IF
      END IF

      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
      RETURN

      // End of CGEEV

      }
