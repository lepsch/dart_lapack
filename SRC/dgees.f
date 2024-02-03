      SUBROUTINE DGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI, VS, LDVS, WORK, LWORK, BWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVS, SORT;
      int                INFO, LDA, LDVS, LWORK, N, SDIM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      double             A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ), WR( * );
      // ..
      // .. Function Arguments ..
      bool               SELECT;
      // EXTERNAL SELECT
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      // ..
      // .. Local Scalars ..
      bool               CURSL, LASTSL, LQUERY, LST2SL, SCALEA, WANTST, WANTVS       int                HSWORK, I, I1, I2, IBAL, ICOND, IERR, IEVAL, IHI, ILO, INXT, IP, ITAU, IWRK, MAXWRK, MINWRK;;
      double             ANRM, BIGNUM, CSCALE, EPS, S, SEP, SMLNUM;
      // ..
      // .. Local Arrays ..
      int                IDUM( 1 );
      double             DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLACPY, DLASCL, DORGHR, DSWAP, DTRSEN, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, DLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVS = LSAME( JOBVS, 'V' )
      WANTST = LSAME( SORT, 'S' )
      IF( ( .NOT.WANTVS ) .AND. ( .NOT.LSAME( JOBVS, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVS.LT.1 .OR. ( WANTVS .AND. LDVS.LT.N ) ) THEN
         INFO = -11
      END IF

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.
        // HSWORK refers to the workspace preferred by DHSEQR, as
        // calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
       t // he worst case.)

      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            MAXWRK = 2*N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 )
            MINWRK = 3*N

            CALL DHSEQR( 'S', JOBVS, N, 1, N, A, LDA, WR, WI, VS, LDVS, WORK, -1, IEVAL )
            HSWORK = INT( WORK( 1 ) )

            IF( .NOT.WANTVS ) THEN
               MAXWRK = MAX( MAXWRK, N + HSWORK )
            ELSE
               MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1, 'DORGHR', ' ', N, 1, N, -1 ) )
               MAXWRK = MAX( MAXWRK, N + HSWORK )
            END IF
         END IF
         WORK( 1 ) = MAXWRK

         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -13
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEES ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) THEN
         SDIM = 0
         RETURN
      END IF

      // Get machine constants

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = DLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      END IF
      IF( SCALEA ) CALL DLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )

      // Permute the matrix to make it more nearly triangular
      // (Workspace: need N)

      IBAL = 1
      CALL DGEBAL( 'P', N, A, LDA, ILO, IHI, WORK( IBAL ), IERR )

      // Reduce to upper Hessenberg form
      // (Workspace: need 3*N, prefer 2*N+N*NB)

      ITAU = N + IBAL
      IWRK = N + ITAU
      CALL DGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )

      IF( WANTVS ) THEN

         // Copy Householder vectors to VS

         CALL DLACPY( 'L', N, N, A, LDA, VS, LDVS )

         // Generate orthogonal matrix in VS
         // (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)

         CALL DORGHR( N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )
      END IF

      SDIM = 0

      // Perform QR iteration, accumulating Schur vectors in VS if desired
      // (Workspace: need N+1, prefer N+HSWORK (see comments) )

      IWRK = ITAU
      CALL DHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, WR, WI, VS, LDVS, WORK( IWRK ), LWORK-IWRK+1, IEVAL )       IF( IEVAL.GT.0 ) INFO = IEVAL

      // Sort eigenvalues if desired

      IF( WANTST .AND. INFO.EQ.0 ) THEN
         IF( SCALEA ) THEN
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, WR, N, IERR )
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, WI, N, IERR )
         END IF
         DO 10 I = 1, N
            BWORK( I ) = SELECT( WR( I ), WI( I ) )
   10    CONTINUE

         // Reorder eigenvalues and transform Schur vectors
         // (Workspace: none needed)

         CALL DTRSEN( 'N', JOBVS, BWORK, N, A, LDA, VS, LDVS, WR, WI, SDIM, S, SEP, WORK( IWRK ), LWORK-IWRK+1, IDUM, 1, ICOND )
         IF( ICOND.GT.0 ) INFO = N + ICOND
      END IF

      IF( WANTVS ) THEN

         // Undo balancing
         // (Workspace: need N)

         CALL DGEBAK( 'P', 'R', N, ILO, IHI, WORK( IBAL ), N, VS, LDVS, IERR )
      END IF

      IF( SCALEA ) THEN

         // Undo scaling for the Schur form of A

         CALL DLASCL( 'H', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR )
         CALL DCOPY( N, A, LDA+1, WR, 1 )
         IF( CSCALE.EQ.SMLNUM ) THEN

            // If scaling back towards underflow, adjust WI if an
            // offdiagonal element of a 2-by-2 block in the Schur form
            // underflows.

            IF( IEVAL.GT.0 ) THEN
               I1 = IEVAL + 1
               I2 = IHI - 1
               CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, MAX( ILO-1, 1 ), IERR )
            ELSE IF( WANTST ) THEN
               I1 = 1
               I2 = N - 1
            ELSE
               I1 = ILO
               I2 = IHI - 1
            END IF
            INXT = I1 - 1
            DO 20 I = I1, I2
               IF( I.LT.INXT ) GO TO 20
               IF( WI( I ).EQ.ZERO ) THEN
                  INXT = I + 1
               ELSE
                  IF( A( I+1, I ).EQ.ZERO ) THEN
                     WI( I ) = ZERO
                     WI( I+1 ) = ZERO
                  ELSE IF( A( I+1, I ).NE.ZERO .AND. A( I, I+1 ).EQ. ZERO ) THEN
                     WI( I ) = ZERO
                     WI( I+1 ) = ZERO
                     IF( I.GT.1 ) CALL DSWAP( I-1, A( 1, I ), 1, A( 1, I+1 ), 1 )                      IF( N.GT.I+1 ) CALL DSWAP( N-I-1, A( I, I+2 ), LDA, A( I+1, I+2 ), LDA )
                     IF( WANTVS ) THEN
                        CALL DSWAP( N, VS( 1, I ), 1, VS( 1, I+1 ), 1 )
                     END IF
                     A( I, I+1 ) = A( I+1, I )
                     A( I+1, I ) = ZERO
                  END IF
                  INXT = I + 2
               END IF
   20       CONTINUE
         END IF

         // Undo scaling for the imaginary part of the eigenvalues

         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-IEVAL, 1, WI( IEVAL+1 ), MAX( N-IEVAL, 1 ), IERR )
      END IF

      IF( WANTST .AND. INFO.EQ.0 ) THEN

         // Check if reordering successful

         LASTSL = .TRUE.
         LST2SL = .TRUE.
         SDIM = 0
         IP = 0
         DO 30 I = 1, N
            CURSL = SELECT( WR( I ), WI( I ) )
            IF( WI( I ).EQ.ZERO ) THEN
               IF( CURSL ) SDIM = SDIM + 1
               IP = 0
               IF( CURSL .AND. .NOT.LASTSL ) INFO = N + 2
            ELSE
               IF( IP.EQ.1 ) THEN

                  // Last eigenvalue of conjugate pair

                  CURSL = CURSL .OR. LASTSL
                  LASTSL = CURSL
                  IF( CURSL ) SDIM = SDIM + 2
                  IP = -1
                  IF( CURSL .AND. .NOT.LST2SL ) INFO = N + 2
               ELSE

                  // First eigenvalue of conjugate pair

                  IP = 1
               END IF
            END IF
            LST2SL = LASTSL
            LASTSL = CURSL
   30    CONTINUE
      END IF

      WORK( 1 ) = MAXWRK
      RETURN

      // End of DGEES

      END
