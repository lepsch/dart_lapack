      SUBROUTINE SGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI, VS, LDVS, WORK, LWORK, BWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVS, SORT;
      int                INFO, LDA, LDVS, LWORK, N, SDIM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      REAL               A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ), WR( * )
      // ..
      // .. Function Arguments ..
      bool               SELECT;
      // EXTERNAL SELECT
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               CURSL, LASTSL, LQUERY, LST2SL, SCALEA, WANTST, WANTVS;
      int                HSWORK, I, I1, I2, IBAL, ICOND, IERR, IEVAL, IHI, ILO, INXT, IP, ITAU, IWRK, MAXWRK, MINWRK;
      REAL               ANRM, BIGNUM, CSCALE, EPS, S, SEP, SMLNUM
      // ..
      // .. Local Arrays ..
      int                IDUM( 1 );
      REAL               DUM( 1 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEBAK, SGEBAL, SGEHRD, SHSEQR, SLACPY, SLASCL, SORGHR, SSWAP, STRSEN, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANGE, SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SLAMCH, SLANGE, SROUNDUP_LWORK
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
      if ( ( .NOT.WANTVS ) .AND. ( .NOT.LSAME( JOBVS, 'N' ) ) ) {
         INFO = -1
      } else if ( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDVS.LT.1 .OR. ( WANTVS .AND. LDVS.LT.N ) ) {
         INFO = -11
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.
        // HSWORK refers to the workspace preferred by SHSEQR, as
        // calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
       t // he worst case.)

      if ( INFO.EQ.0 ) {
         if ( N.EQ.0 ) {
            MINWRK = 1
            MAXWRK = 1
         } else {
            MAXWRK = 2*N + N*ILAENV( 1, 'SGEHRD', ' ', N, 1, N, 0 )
            MINWRK = 3*N

            CALL SHSEQR( 'S', JOBVS, N, 1, N, A, LDA, WR, WI, VS, LDVS, WORK, -1, IEVAL )
            HSWORK = INT( WORK( 1 ) )

            if ( .NOT.WANTVS ) {
               MAXWRK = MAX( MAXWRK, N + HSWORK )
            } else {
               MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1, 'SORGHR', ' ', N, 1, N, -1 ) )
               MAXWRK = MAX( MAXWRK, N + HSWORK )
            }
         }
         WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)

         if ( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) {
            INFO = -13
         }
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SGEES ', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         SDIM = 0
         RETURN
      }

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = SLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      } else if ( ANRM.GT.BIGNUM ) {
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      }
      IF( SCALEA ) CALL SLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )

      // Permute the matrix to make it more nearly triangular
      // (Workspace: need N)

      IBAL = 1
      CALL SGEBAL( 'P', N, A, LDA, ILO, IHI, WORK( IBAL ), IERR )

      // Reduce to upper Hessenberg form
      // (Workspace: need 3*N, prefer 2*N+N*NB)

      ITAU = N + IBAL
      IWRK = N + ITAU
      CALL SGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )

      if ( WANTVS ) {

         // Copy Householder vectors to VS

         CALL SLACPY( 'L', N, N, A, LDA, VS, LDVS )

         // Generate orthogonal matrix in VS
         // (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)

         CALL SORGHR( N, ILO, IHI, VS, LDVS, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )
      }

      SDIM = 0

      // Perform QR iteration, accumulating Schur vectors in VS if desired
      // (Workspace: need N+1, prefer N+HSWORK (see comments) )

      IWRK = ITAU
      CALL SHSEQR( 'S', JOBVS, N, ILO, IHI, A, LDA, WR, WI, VS, LDVS, WORK( IWRK ), LWORK-IWRK+1, IEVAL )       IF( IEVAL.GT.0 ) INFO = IEVAL

      // Sort eigenvalues if desired

      if ( WANTST .AND. INFO.EQ.0 ) {
         if ( SCALEA ) {
            CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, WR, N, IERR )
            CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, WI, N, IERR )
         }
         DO 10 I = 1, N
            BWORK( I ) = SELECT( WR( I ), WI( I ) )
   10    CONTINUE

         // Reorder eigenvalues and transform Schur vectors
         // (Workspace: none needed)

         CALL STRSEN( 'N', JOBVS, BWORK, N, A, LDA, VS, LDVS, WR, WI, SDIM, S, SEP, WORK( IWRK ), LWORK-IWRK+1, IDUM, 1, ICOND )
         IF( ICOND.GT.0 ) INFO = N + ICOND
      }

      if ( WANTVS ) {

         // Undo balancing
         // (Workspace: need N)

         CALL SGEBAK( 'P', 'R', N, ILO, IHI, WORK( IBAL ), N, VS, LDVS, IERR )
      }

      if ( SCALEA ) {

         // Undo scaling for the Schur form of A

         CALL SLASCL( 'H', 0, 0, CSCALE, ANRM, N, N, A, LDA, IERR )
         CALL SCOPY( N, A, LDA+1, WR, 1 )
         if ( CSCALE.EQ.SMLNUM ) {

            // If scaling back towards underflow, adjust WI if an
            // offdiagonal element of a 2-by-2 block in the Schur form
            // underflows.

            if ( IEVAL.GT.0 ) {
               I1 = IEVAL + 1
               I2 = IHI - 1
               CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, MAX( ILO-1, 1 ), IERR )
            } else if ( WANTST ) {
               I1 = 1
               I2 = N - 1
            } else {
               I1 = ILO
               I2 = IHI - 1
            }
            INXT = I1 - 1
            DO 20 I = I1, I2
               IF( I.LT.INXT ) GO TO 20
               if ( WI( I ).EQ.ZERO ) {
                  INXT = I + 1
               } else {
                  if ( A( I+1, I ).EQ.ZERO ) {
                     WI( I ) = ZERO
                     WI( I+1 ) = ZERO
                  } else if ( A( I+1, I ).NE.ZERO .AND. A( I, I+1 ).EQ. ZERO ) {
                     WI( I ) = ZERO
                     WI( I+1 ) = ZERO
                     IF( I.GT.1 ) CALL SSWAP( I-1, A( 1, I ), 1, A( 1, I+1 ), 1 )                      IF( N.GT.I+1 ) CALL SSWAP( N-I-1, A( I, I+2 ), LDA, A( I+1, I+2 ), LDA )
                     if ( WANTVS ) {
                        CALL SSWAP( N, VS( 1, I ), 1, VS( 1, I+1 ), 1 )
                     }
                     A( I, I+1 ) = A( I+1, I )
                     A( I+1, I ) = ZERO
                  }
                  INXT = I + 2
               }
   20       CONTINUE
         }

         // Undo scaling for the imaginary part of the eigenvalues

         CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, N-IEVAL, 1, WI( IEVAL+1 ), MAX( N-IEVAL, 1 ), IERR )
      }

      if ( WANTST .AND. INFO.EQ.0 ) {

         // Check if reordering successful

         LASTSL = .TRUE.
         LST2SL = .TRUE.
         SDIM = 0
         IP = 0
         DO 30 I = 1, N
            CURSL = SELECT( WR( I ), WI( I ) )
            if ( WI( I ).EQ.ZERO ) {
               IF( CURSL ) SDIM = SDIM + 1
               IP = 0
               IF( CURSL .AND. .NOT.LASTSL ) INFO = N + 2
            } else {
               if ( IP.EQ.1 ) {

                  // Last eigenvalue of conjugate pair

                  CURSL = CURSL .OR. LASTSL
                  LASTSL = CURSL
                  IF( CURSL ) SDIM = SDIM + 2
                  IP = -1
                  IF( CURSL .AND. .NOT.LST2SL ) INFO = N + 2
               } else {

                  // First eigenvalue of conjugate pair

                  IP = 1
               }
            }
            LST2SL = LASTSL
            LASTSL = CURSL
   30    CONTINUE
      }

      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
      RETURN

      // End of SGEES

      }
