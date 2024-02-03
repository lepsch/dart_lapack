      SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
      implicit none

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             BALANC, JOBVL, JOBVR, SENSE;
      int                IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N;
      double             ABNRM;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             A( LDA, * ), RCONDE( * ), RCONDV( * ), SCALE( * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), WORK( * ), WR( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, SCALEA, WANTVL, WANTVR, WNTSNB, WNTSNE, WNTSNN, WNTSNV;
      String             JOB, SIDE;
      int                HSWORK, I, ICOND, IERR, ITAU, IWRK, K, LWORK_TREVC, MAXWRK, MINWRK, NOUT;
      double             ANRM, BIGNUM, CS, CSCALE, EPS, R, SCL, SMLNUM, SN;
      // ..
      // .. Local Arrays ..
      bool               SELECT( 1 );
      double             DUM( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLACPY, DLARTG, DLASCL, DORGHR, DROT, DSCAL, DTREVC3, DTRSNA, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                IDAMAX, ILAENV;
      double             DLAMCH, DLANGE, DLAPY2, DNRM2;
      // EXTERNAL LSAME, IDAMAX, ILAENV, DLAMCH, DLANGE, DLAPY2, DNRM2
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
      WNTSNN = LSAME( SENSE, 'N' )
      WNTSNE = LSAME( SENSE, 'E' )
      WNTSNV = LSAME( SENSE, 'V' )
      WNTSNB = LSAME( SENSE, 'B' )
      if ( .NOT.( LSAME( BALANC, 'N' ) .OR. LSAME( BALANC, 'S' ) .OR. LSAME( BALANC, 'P' ) .OR. LSAME( BALANC, 'B' ) ) ) {
         INFO = -1
      } else if ( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( JOBVL, 'N' ) ) ) {
         INFO = -2
      } else if ( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( JOBVR, 'N' ) ) ) {
         INFO = -3
      } else if ( .NOT.( WNTSNN .OR. WNTSNE .OR. WNTSNB .OR. WNTSNV ) .OR. ( ( WNTSNE .OR. WNTSNB ) .AND. .NOT.( WANTVL .AND. WANTVR ) ) ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) {
         INFO = -11
      } else if ( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) {
         INFO = -13
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.
        // HSWORK refers to the workspace preferred by DHSEQR, as
        // calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
        // the worst case.)

      if ( INFO.EQ.0 ) {
         if ( N.EQ.0 ) {
            MINWRK = 1
            MAXWRK = 1
         } else {
            MAXWRK = N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 )

            if ( WANTVL ) {
               CALL DTREVC3( 'L', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, IERR )
               LWORK_TREVC = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
               CALL DHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VL, LDVL, WORK, -1, INFO )
            } else if ( WANTVR ) {
               CALL DTREVC3( 'R', 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK, -1, IERR )
               LWORK_TREVC = INT( WORK(1) )
               MAXWRK = MAX( MAXWRK, N + LWORK_TREVC )
               CALL DHSEQR( 'S', 'V', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO )
            } else {
               if ( WNTSNN ) {
                  CALL DHSEQR( 'E', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO )
               } else {
                  CALL DHSEQR( 'S', 'N', N, 1, N, A, LDA, WR, WI, VR, LDVR, WORK, -1, INFO )
               }
            }
            HSWORK = INT( WORK(1) )

            if ( ( .NOT.WANTVL ) .AND. ( .NOT.WANTVR ) ) {
               MINWRK = 2*N
               IF( .NOT.WNTSNN ) MINWRK = MAX( MINWRK, N*N+6*N )
               MAXWRK = MAX( MAXWRK, HSWORK )
               IF( .NOT.WNTSNN ) MAXWRK = MAX( MAXWRK, N*N + 6*N )
            } else {
               MINWRK = 3*N
               IF( ( .NOT.WNTSNN ) .AND. ( .NOT.WNTSNE ) ) MINWRK = MAX( MINWRK, N*N + 6*N )
               MAXWRK = MAX( MAXWRK, HSWORK )
               MAXWRK = MAX( MAXWRK, N + ( N - 1 )*ILAENV( 1, 'DORGHR', ' ', N, 1, N, -1 ) )                IF( ( .NOT.WNTSNN ) .AND. ( .NOT.WNTSNE ) ) MAXWRK = MAX( MAXWRK, N*N + 6*N )
               MAXWRK = MAX( MAXWRK, 3*N )
            }
            MAXWRK = MAX( MAXWRK, MINWRK )
         }
         WORK( 1 ) = MAXWRK

         if ( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) {
            INFO = -21
         }
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DGEEVX', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Get machine constants

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ICOND = 0
      ANRM = DLANGE( 'M', N, N, A, LDA, DUM )
      SCALEA = .FALSE.
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         SCALEA = .TRUE.
         CSCALE = SMLNUM
      } else if ( ANRM.GT.BIGNUM ) {
         SCALEA = .TRUE.
         CSCALE = BIGNUM
      }
      IF( SCALEA ) CALL DLASCL( 'G', 0, 0, ANRM, CSCALE, N, N, A, LDA, IERR )

      // Balance the matrix and compute ABNRM

      CALL DGEBAL( BALANC, N, A, LDA, ILO, IHI, SCALE, IERR )
      ABNRM = DLANGE( '1', N, N, A, LDA, DUM )
      if ( SCALEA ) {
         DUM( 1 ) = ABNRM
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, 1, 1, DUM, 1, IERR )
         ABNRM = DUM( 1 )
      }

      // Reduce to upper Hessenberg form
      // (Workspace: need 2*N, prefer N+N*NB)

      ITAU = 1
      IWRK = ITAU + N
      CALL DGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )

      if ( WANTVL ) {

         // Want left eigenvectors
         // Copy Householder vectors to VL

         SIDE = 'L'
         CALL DLACPY( 'L', N, N, A, LDA, VL, LDVL )

         // Generate orthogonal matrix in VL
         // (Workspace: need 2*N-1, prefer N+(N-1)*NB)

         CALL DORGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )

         // Perform QR iteration, accumulating Schur vectors in VL
         // (Workspace: need 1, prefer HSWORK (see comments) )

         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL, WORK( IWRK ), LWORK-IWRK+1, INFO )

         if ( WANTVR ) {

            // Want left and right eigenvectors
            // Copy Schur vectors to VR

            SIDE = 'B'
            CALL DLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
         }

      } else if ( WANTVR ) {

         // Want right eigenvectors
         // Copy Householder vectors to VR

         SIDE = 'R'
         CALL DLACPY( 'L', N, N, A, LDA, VR, LDVR )

         // Generate orthogonal matrix in VR
         // (Workspace: need 2*N-1, prefer N+(N-1)*NB)

         CALL DORGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ), LWORK-IWRK+1, IERR )

         // Perform QR iteration, accumulating Schur vectors in VR
         // (Workspace: need 1, prefer HSWORK (see comments) )

         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO )

      } else {

         // Compute eigenvalues only
         // If condition numbers desired, compute Schur form

         if ( WNTSNN ) {
            JOB = 'E'
         } else {
            JOB = 'S'
         }

         // (Workspace: need 1, prefer HSWORK (see comments) )

         IWRK = ITAU
         CALL DHSEQR( JOB, 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR, WORK( IWRK ), LWORK-IWRK+1, INFO )
      }

      // If INFO .NE. 0 from DHSEQR, then quit

      IF( INFO.NE.0 ) GO TO 50

      if ( WANTVL .OR. WANTVR ) {

         // Compute left and/or right eigenvectors
         // (Workspace: need 3*N, prefer N + 2*N*NB)

         CALL DTREVC3( SIDE, 'B', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, N, NOUT, WORK( IWRK ), LWORK-IWRK+1, IERR )
      }

      // Compute condition numbers if desired
      // (Workspace: need N*N+6*N unless SENSE = 'E')

      if ( .NOT.WNTSNN ) {
         CALL DTRSNA( SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR, LDVR, RCONDE, RCONDV, N, NOUT, WORK( IWRK ), N, IWORK, ICOND )
      }

      if ( WANTVL ) {

         // Undo balancing of left eigenvectors

         CALL DGEBAK( BALANC, 'L', N, ILO, IHI, SCALE, N, VL, LDVL, IERR )

         // Normalize left eigenvectors and make largest component real

         DO 20 I = 1, N
            if ( WI( I ).EQ.ZERO ) {
               SCL = ONE / DNRM2( N, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
            } else if ( WI( I ).GT.ZERO ) {
               SCL = ONE / DLAPY2( DNRM2( N, VL( 1, I ), 1 ), DNRM2( N, VL( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I+1 ), 1 )
               DO 10 K = 1, N
                  WORK( K ) = VL( K, I )**2 + VL( K, I+1 )**2
   10          CONTINUE
               K = IDAMAX( N, WORK, 1 )
               CALL DLARTG( VL( K, I ), VL( K, I+1 ), CS, SN, R )
               CALL DROT( N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN )
               VL( K, I+1 ) = ZERO
            }
   20    CONTINUE
      }

      if ( WANTVR ) {

         // Undo balancing of right eigenvectors

         CALL DGEBAK( BALANC, 'R', N, ILO, IHI, SCALE, N, VR, LDVR, IERR )

         // Normalize right eigenvectors and make largest component real

         DO 40 I = 1, N
            if ( WI( I ).EQ.ZERO ) {
               SCL = ONE / DNRM2( N, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
            } else if ( WI( I ).GT.ZERO ) {
               SCL = ONE / DLAPY2( DNRM2( N, VR( 1, I ), 1 ), DNRM2( N, VR( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I+1 ), 1 )
               DO 30 K = 1, N
                  WORK( K ) = VR( K, I )**2 + VR( K, I+1 )**2
   30          CONTINUE
               K = IDAMAX( N, WORK, 1 )
               CALL DLARTG( VR( K, I ), VR( K, I+1 ), CS, SN, R )
               CALL DROT( N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN )
               VR( K, I+1 ) = ZERO
            }
   40    CONTINUE
      }

      // Undo scaling if necessary

   50 CONTINUE
      if ( SCALEA ) {
         CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WR( INFO+1 ), MAX( N-INFO, 1 ), IERR )          CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N-INFO, 1, WI( INFO+1 ), MAX( N-INFO, 1 ), IERR )
         if ( INFO.EQ.0 ) {
            IF( ( WNTSNV .OR. WNTSNB ) .AND. ICOND.EQ.0 ) CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, N, 1, RCONDV, N, IERR )
         } else {
            CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WR, N, IERR )             CALL DLASCL( 'G', 0, 0, CSCALE, ANRM, ILO-1, 1, WI, N, IERR )
         }
      }

      WORK( 1 ) = MAXWRK
      RETURN

      // End of DGEEVX

      }
