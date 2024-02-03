      SUBROUTINE SGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
      implicit none

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVL, JOBVR;
      int                INFO, LDA, LDVL, LDVR, LWORK, N;
      // ..
      // .. Array Arguments ..
      REAL   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), WORK( * ), WR( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL   ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, SCALEA, WANTVL, WANTVR;
      String             SIDE;
      int                HSWORK, I, IBAL, IERR, IHI, ILO, ITAU, IWRK, K, LWORK_TREVC, MAXWRK, MINWRK, NOUT;
      REAL   ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM, SN;
      // ..
      // .. Local Arrays ..
      bool               SELECT( 1 );
      REAL   DUM( 1 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEBAK, SGEBAL, SGEHRD, SHSEQR, SLACPY, SLARTG, SLASCL, SORGHR, SROT, SSCAL, STREVC3, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ISAMAX, ILAENV;
      REAL               SLAMCH, SLANGE, SLAPY2, SNRM2, SROUNDUP_LWORK
      // EXTERNAL LSAME, ISAMAX, ILAENV, SLAMCH, SLANGE, SLAPY2, SNRM2, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      WANTVL = LSAME( JOBVL, 'V' )
      WANTVR = LSAME( JOBVR, 'V' )
      if ( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) {
         INFO = -1
      } else if ( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) {
         INFO = -9
      } else if ( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) {
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
        // the worst case.)

      if ( INFO.EQ.0 ) {
         if ( N.EQ.0 ) {
            MINWRK = 1
            MAXWRK = 1
         } else {
            MAXWRK = 2*N + N*ILAENV( 1, 'SGEHRD', ' ', N, 1, N, 0 )
            if ( WANTVL ) {
               MINWRK = 4*N
               MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1, 'SORGHR', ' ', N, 1, N, -1 ) )                CALL SHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VL, LDVL, WORK, -1, INFO )
               HSWORK = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
               CALL STREVC3( 'L', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, IERR )
               LWORK_TREVC = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
               MAXWRK = MAX( MAXWRK, 4*N )
            } else if ( WANTVR ) {
               MINWRK = 4*N
               MAXWRK = MAX( MAXWRK, 2*N + ( N - 1 )*ILAENV( 1, 'SORGHR', ' ', N, 1, N, -1 ) )                CALL SHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO )
               HSWORK = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
               CALL STREVC3( 'R', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, IERR )
               LWORK_TREVC = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
               MAXWRK = MAX( MAXWRK, 4*N )
            } else {
               MINWRK = 3*N
               CALL SHSEQR( 'E', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO )
               HSWORK = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + 1, N + HSWORK )
            }
            MAXWRK = MAX( MAXWRK, MINWRK )
         }
         WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)

         if ( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) {
            INFO = -13
         }
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SGEEV ', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

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

      // Balance the matrix
      // (Workspace: need N)

      IBAL = 1
      CALL SGEBAL( 'B', N, A, LDA, ILO, IHI, WORK( IBAL ), IERR )

      // Reduce to upper Hessenberg form
      // (Workspace: need 3*N, prefer 2*N+N*NB)

      ITAU = IBAL + N
      IWRK = ITAU + N
      CALL SGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )

      if ( WANTVL ) {

         // Want left eigenvectors
         // Copy Householder vectors to VL

         SIDE = 'L'
         CALL SLACPY( 'L', N, N, A, LDA, VL, LDVL )

         // Generate orthogonal matrix in VL
         // (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)

         CALL SORGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )

         // Perform QR iteration, accumulating Schur vectors in VL
         // (Workspace: need N+1, prefer N+HSWORK (see comments) )

         IWRK = ITAU
         CALL SHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL, WORK( IWRK ), LWORK-IWRK+1, INFO )

         if ( WANTVR ) {

            // Want left and right eigenvectors
            // Copy Schur vectors to VR

            SIDE = 'B'
            CALL SLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         }

      } else if ( WANTVR ) {

         // Want right eigenvectors
         // Copy Householder vectors to VR

         SIDE = 'R'
         CALL SLACPY( 'L', N, N, A, LDA, VR, LDVR )

         // Generate orthogonal matrix in VR
         // (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)

         CALL SORGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )

         // Perform QR iteration, accumulating Schur vectors in VR
         // (Workspace: need N+1, prefer N+HSWORK (see comments) )

         IWRK = ITAU
         CALL SHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO )

      } else {

         // Compute eigenvalues only
         // (Workspace: need N+1, prefer N+HSWORK (see comments) )

         IWRK = ITAU
         CALL SHSEQR( 'E', 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO )
      }

      // If INFO .NE. 0 from SHSEQR, then quit

      IF( INFO.NE.0 ) GO TO 50

      if ( WANTVL .OR. WANTVR ) {

         // Compute left and/or right eigenvectors
         // (Workspace: need 4*N, prefer N + N + 2*N*NB)

         CALL STREVC3( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK( IWRK ), LWORK-IWRK+1, IERR )
      }

      if ( WANTVL ) {

         // Undo balancing of left eigenvectors
         // (Workspace: need N)

         CALL SGEBAK( 'B', 'L', N, ILO, IHI, WORK( IBAL ), N, VL, LDVL, IERR )

         // Normalize left eigenvectors and make largest component real

         DO 20 I = 1, N
            if ( WI( I ).EQ.ZERO ) {
               SCL = ONE / SNRM2( N, VL( 1, I ), 1 )
               CALL SSCAL( N, SCL, VL( 1, I ), 1 )
            } else if ( WI( I ).GT.ZERO ) {
               SCL = ONE / SLAPY2( SNRM2( N, VL( 1, I ), 1 ), SNRM2( N, VL( 1, I+1 ), 1 ) )
               CALL SSCAL( N, SCL, VL( 1, I ), 1 )
               CALL SSCAL( N, SCL, VL( 1, I+1 ), 1 )
               DO 10 K = 1, N
                  WORK( IWRK+K-1 ) = VL( K, I )**2 + VL( K, I+1 )**2
   10          CONTINUE
               K = ISAMAX( N, WORK( IWRK ), 1 )
               CALL SLARTG( VL( K, I ), VL( K, I+1 ), CS, SN, R )
               CALL SROT( N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN )
               VL( K, I+1 ) = ZERO
            }
   20    CONTINUE
      }

      if ( WANTVR ) {

         // Undo balancing of right eigenvectors
         // (Workspace: need N)

         CALL SGEBAK( 'B', 'R', N, ILO, IHI, WORK( IBAL ), N, VR, LDVR, IERR )

         // Normalize right eigenvectors and make largest component real

         DO 40 I = 1, N
            if ( WI( I ).EQ.ZERO ) {
               SCL = ONE / SNRM2( N, VR( 1, I ), 1 )
               CALL SSCAL( N, SCL, VR( 1, I ), 1 )
            } else if ( WI( I ).GT.ZERO ) {
               SCL = ONE / SLAPY2( SNRM2( N, VR( 1, I ), 1 ), SNRM2( N, VR( 1, I+1 ), 1 ) )
               CALL SSCAL( N, SCL, VR( 1, I ), 1 )
               CALL SSCAL( N, SCL, VR( 1, I+1 ), 1 )
               DO 30 K = 1, N
                  WORK( IWRK+K-1 ) = VR( K, I )**2 + VR( K, I+1 )**2
   30          CONTINUE
               K = ISAMAX( N, WORK( IWRK ), 1 )
               CALL SLARTG( VR( K, I ), VR( K, I+1 ), CS, SN, R )
               CALL SROT( N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN )
               VR( K, I+1 ) = ZERO
            }
   40    CONTINUE
      }

      // Undo scaling if necessary

   50 CONTINUE
      if ( SCALEA ) {
         CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WR( INFO+1 ), MAX( N-INFO, 1 ), IERR )          CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WI( INFO+1 ), MAX( N-INFO, 1 ), IERR )
         if ( INFO.GT.0 ) {
            CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N, IERR )             CALL SLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N, IERR )
         }
      }

      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
      RETURN

      // End of SGEEV

      }
