      SUBROUTINE ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, INFO )
      implicit none
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             BALANC, JOBVL, JOBVR, SENSE;
      int                IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N;
      double             ABNRM;
      // ..
      // .. Array Arguments ..
      double             RCONDE( * ), RCONDV( * ), RWORK( * ), SCALE( * )       COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), WORK( * );
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      // ..
      // .. Local Scalars ..
      bool               LQUERY, SCALEA, WANTVL, WANTVR, WNTSNB, WNTSNE, WNTSNN, WNTSNV;
      String             JOB, SIDE;
      int                HSWORK, I, ICOND, IERR, ITAU, IWRK, K, LWORK_TREVC, MAXWRK, MINWRK, NOUT;
      double             ANRM, BIGNUM, CSCALE, EPS, SCL, SMLNUM;
      COMPLEX*16         TMP
      // ..
      // .. Local Arrays ..
      bool               SELECT( 1 );
      double             DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASCL, XERBLA, ZDSCAL, ZGEBAK, ZGEBAL, ZGEHRD, ZHSEQR, ZLACPY, ZLASCL, ZSCAL, ZTREVC3, ZTRSNA, ZUNGHR
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX, ILAENV;
      double             DLAMCH, DZNRM2, ZLANGE;
      // EXTERNAL LSAME, IDAMAX, ILAENV, DLAMCH, DZNRM2, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, CONJG, AIMAG, MAX, SQRT
      // ..
      // .. Executable Statements ..
*
      // Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVL = LSAME( JOBVL, 'V' )
      WANTVR = LSAME( JOBVR, 'V' )
      WNTSNN = LSAME( SENSE, 'N' )
      WNTSNE = LSAME( SENSE, 'E' )
      WNTSNV = LSAME( SENSE, 'V' )
      WNTSNB = LSAME( SENSE, 'B' )
      IF( .NOT.( LSAME( BALANC, 'N' ) .OR. LSAME( BALANC, 'S' ) .OR. LSAME( BALANC, 'P' ) .OR. LSAME( BALANC, 'B' ) ) ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT.( WNTSNN .OR. WNTSNE .OR. WNTSNB .OR. WNTSNV ) .OR. ( ( WNTSNE .OR. WNTSNB ) .AND. .NOT.( WANTVL .AND. WANTVR ) ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
         INFO = -10
      ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
         INFO = -12
      END IF
*
      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // CWorkspace refers to complex workspace, and RWorkspace to real
        // workspace. NB refers to the optimal block size for the
        // immediately following subroutine, as returned by ILAENV.
        // HSWORK refers to the workspace preferred by ZHSEQR, as
        // calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
       t // he worst case.)
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            MAXWRK = N + N*ILAENV( 1, 'ZGEHRD', ' ', N, 1, N, 0 )
*
            IF( WANTVL ) THEN
               CALL ZTREVC3( 'L', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, RWORK, -1, IERR )
               LWORK_TREVC = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, LWORK_TREVC )
               CALL ZHSEQR( 'S', 'V', N, 1, N, A, LDA, W, VL, LDVL, WORK, -1, INFO )
            ELSE IF( WANTVR ) THEN
               CALL ZTREVC3( 'R', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, RWORK, -1, IERR )
               LWORK_TREVC = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, LWORK_TREVC )
               CALL ZHSEQR( 'S', 'V', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO )
            ELSE
               IF( WNTSNN ) THEN
                  CALL ZHSEQR( 'E', 'N', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO )
               ELSE
                  CALL ZHSEQR( 'S', 'N', N, 1, N, A, LDA, W, VR, LDVR, WORK, -1, INFO )
               END IF
            END IF
            HSWORK = INT( WORK(1) )
*
            IF( ( .NOT.WANTVL ) .AND. ( .NOT.WANTVR ) ) THEN
               MINWRK = 2*N
               IF( .NOT.( WNTSNN .OR. WNTSNE ) ) MINWRK = MAX( MINWRK, N*N + 2*N )
               MAXWRK = MAX( MAXWRK, HSWORK )
               IF( .NOT.( WNTSNN .OR. WNTSNE ) ) MAXWRK = MAX( MAXWRK, N*N + 2*N )
            ELSE
               MINWRK = 2*N
               IF( .NOT.( WNTSNN .OR. WNTSNE ) ) MINWRK = MAX( MINWRK, N*N + 2*N )
               MAXWRK = MAX( MAXWRK, HSWORK )
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'ZUNGHR', ' ', N, 1, N, -1 ) )                IF( .NOT.( WNTSNN .OR. WNTSNE ) ) MAXWRK = MAX( MAXWRK, N*N + 2*N )
               MAXWRK = MAX( MAXWRK, 2*N )
            END IF
            MAXWRK = MAX( MAXWRK, MINWRK )
         END IF
         WORK( 1 ) = MAXWRK
*
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -20
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEEVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      // Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
      // Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ICOND = 0
      ANRM = ZLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA ) CALL ZLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )
*
      // Balance the matrix and compute ABNRM
*
      CALL ZGEBAL( BALANC, N, A, LDA, ILO, IHI, SCALE, IERR )
      ABNRM = ZLANGE( '1', N, N, A, LDA, DUM )
      IF( SCALEA ) THEN
         DUM( 1 ) = ABNRM
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, 1, 1, DUM, 1, IERR )
         ABNRM = DUM( 1 )
      END IF
*
      // Reduce to upper Hessenberg form
      // (CWorkspace: need 2*N, prefer N+N*NB)
      // (RWorkspace: none)
*
      ITAU = 1
      IWRK = ITAU + N
      CALL ZGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )
*
      IF( WANTVL ) THEN
*
         // Want left eigenvectors
         // Copy Householder vectors to VL
*
         SIDE = 'L'
         CALL ZLACPY( 'L', N, N, A, LDA, VL, LDVL )
*
         // Generate unitary matrix in VL
         // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
         // (RWorkspace: none)
*
         CALL ZUNGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )
*
         // Perform QR iteration, accumulating Schur vectors in VL
         // (CWorkspace: need 1, prefer HSWORK (see comments) )
         // (RWorkspace: none)
*
         IWRK = ITAU
         CALL ZHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, W, VL, LDVL, WORK( IWRK ), LWORK-IWRK+1, INFO )
*
         IF( WANTVR ) THEN
*
            // Want left and right eigenvectors
            // Copy Schur vectors to VR
*
            SIDE = 'B'
            CALL ZLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         END IF
*
      ELSE IF( WANTVR ) THEN
*
         // Want right eigenvectors
         // Copy Householder vectors to VR
*
         SIDE = 'R'
         CALL ZLACPY( 'L', N, N, A, LDA, VR, LDVR )
*
         // Generate unitary matrix in VR
         // (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
         // (RWorkspace: none)
*
         CALL ZUNGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )
*
         // Perform QR iteration, accumulating Schur vectors in VR
         // (CWorkspace: need 1, prefer HSWORK (see comments) )
         // (RWorkspace: none)
*
         IWRK = ITAU
         CALL ZHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, W, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO )
*
      ELSE
*
         // Compute eigenvalues only
         // If condition numbers desired, compute Schur form
*
         IF( WNTSNN ) THEN
            JOB = 'E'
         ELSE
            JOB = 'S'
         END IF
*
         // (CWorkspace: need 1, prefer HSWORK (see comments) )
         // (RWorkspace: none)
*
         IWRK = ITAU
         CALL ZHSEQR( JOB, 'N', N, ILO, IHI, A, LDA, W, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO )
      END IF
*
      // If INFO .NE. 0 from ZHSEQR, then quit
*
      IF( INFO.NE.0 ) GO TO 50
*
      IF( WANTVL .OR. WANTVR ) THEN
*
         // Compute left and/or right eigenvectors
         // (CWorkspace: need 2*N, prefer N + 2*N*NB)
         // (RWorkspace: need N)
*
         CALL ZTREVC3( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK( IWRK ), LWORK-IWRK+1, RWORK, N, IERR )
      END IF
*
      // Compute condition numbers if desired
      // (CWorkspace: need N*N+2*N unless SENSE = 'E')
      // (RWorkspace: need 2*N unless SENSE = 'E')
*
      IF( .NOT.WNTSNN ) THEN
         CALL ZTRSNA( SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, RCONDE, RCONDV, N, NOUT, WORK( IWRK ), N, RWORK, ICOND )
      END IF
*
      IF( WANTVL ) THEN
*
         // Undo balancing of left eigenvectors
*
         CALL ZGEBAK( BALANC, 'L', N, ILO, IHI, SCALE, N, VL, LDVL, IERR )
*
         // Normalize left eigenvectors and make largest component real
*
         DO 20 I = 1, N
            SCL = ONE / DZNRM2( N, VL( 1, I ), 1 )
            CALL ZDSCAL( N, SCL, VL( 1, I ), 1 )
            DO 10 K = 1, N
               RWORK( K ) = DBLE( VL( K, I ) )**2 + AIMAG( VL( K, I ) )**2
   10       CONTINUE
            K = IDAMAX( N, RWORK, 1 )
            TMP = CONJG( VL( K, I ) ) / SQRT( RWORK( K ) )
            CALL ZSCAL( N, TMP, VL( 1, I ), 1 )
            VL( K, I ) = DCMPLX( DBLE( VL( K, I ) ), ZERO )
   20    CONTINUE
      END IF
*
      IF( WANTVR ) THEN
*
         // Undo balancing of right eigenvectors
*
         CALL ZGEBAK( BALANC, 'R', N, ILO, IHI, SCALE, N, VR, LDVR, IERR )
*
         // Normalize right eigenvectors and make largest component real
*
         DO 40 I = 1, N
            SCL = ONE / DZNRM2( N, VR( 1, I ), 1 )
            CALL ZDSCAL( N, SCL, VR( 1, I ), 1 )
            DO 30 K = 1, N
               RWORK( K ) = DBLE( VR( K, I ) )**2 + AIMAG( VR( K, I ) )**2
   30       CONTINUE
            K = IDAMAX( N, RWORK, 1 )
            TMP = CONJG( VR( K, I ) ) / SQRT( RWORK( K ) )
            CALL ZSCAL( N, TMP, VR( 1, I ), 1 )
            VR( K, I ) = DCMPLX( DBLE( VR( K, I ) ), ZERO )
   40    CONTINUE
      END IF
*
      // Undo scaling if necessary
*
   50 CONTINUE
      IF( SCALEA ) THEN
         CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, W( INFO+1 ), MAX( N-INFO, 1 ), IERR )
         IF( INFO.EQ.0 ) THEN
            IF( ( WNTSNV .OR. WNTSNB ) .AND. ICOND.EQ.0 ) CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, RCONDV, N, IERR )
         ELSE
            CALL ZLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, W, N, IERR )
         END IF
      END IF
*
      WORK( 1 ) = MAXWRK
      RETURN
*
      // End of ZGEEVX
*
      END
