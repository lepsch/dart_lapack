      SUBROUTINE ZGGESX( JOBVSL, JOBVSR, SORT, SELCTG, SENSE, N, A, LDA, B, LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, RCONDE, RCONDV, WORK, LWORK, RWORK, IWORK, LIWORK, BWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVSL, JOBVSR, SENSE, SORT;
      int                INFO, LDA, LDB, LDVSL, LDVSR, LIWORK, LWORK, N, SDIM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      int                IWORK( * );
      double             RCONDE( 2 ), RCONDV( 2 ), RWORK( * );
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * )
      // ..
      // .. Function Arguments ..
      bool               SELCTG;
      // EXTERNAL SELCTG
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL, LQUERY, WANTSB, WANTSE, WANTSN, WANTST, WANTSV;
      int                I, ICOLS, IERR, IHI, IJOB, IJOBVL, IJOBVR, ILEFT, ILO, IRIGHT, IROWS, IRWRK, ITAU, IWRK, LIWMIN, LWRK, MAXWRK, MINWRK;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PL, PR, SMLNUM;
      // ..
      // .. Local Arrays ..
      double             DIF( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQRF, ZGGBAK, ZGGBAL, ZGGHRD, ZHGEQZ, ZLACPY, ZLASCL, ZLASET, ZTGSEN, ZUNGQR, ZUNMQR
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Decode the input arguments

      if ( LSAME( JOBVSL, 'N' ) ) {
         IJOBVL = 1
         ILVSL = .FALSE.
      } else if ( LSAME( JOBVSL, 'V' ) ) {
         IJOBVL = 2
         ILVSL = .TRUE.
      } else {
         IJOBVL = -1
         ILVSL = .FALSE.
      }

      if ( LSAME( JOBVSR, 'N' ) ) {
         IJOBVR = 1
         ILVSR = .FALSE.
      } else if ( LSAME( JOBVSR, 'V' ) ) {
         IJOBVR = 2
         ILVSR = .TRUE.
      } else {
         IJOBVR = -1
         ILVSR = .FALSE.
      }

      WANTST = LSAME( SORT, 'S' )
      WANTSN = LSAME( SENSE, 'N' )
      WANTSE = LSAME( SENSE, 'E' )
      WANTSV = LSAME( SENSE, 'V' )
      WANTSB = LSAME( SENSE, 'B' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
      if ( WANTSN ) {
         IJOB = 0
      } else if ( WANTSE ) {
         IJOB = 1
      } else if ( WANTSV ) {
         IJOB = 2
      } else if ( WANTSB ) {
         IJOB = 4
      }

      // Test the input arguments

      INFO = 0
      if ( IJOBVL.LE.0 ) {
         INFO = -1
      } else if ( IJOBVR.LE.0 ) {
         INFO = -2
      } else if ( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) {
         INFO = -3
      } else if ( .NOT.( WANTSN .OR. WANTSE .OR. WANTSV .OR. WANTSB ) .OR. ( .NOT.WANTST .AND. .NOT.WANTSN ) ) {
         INFO = -5
      } else if ( N.LT.0 ) {
         INFO = -6
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -10
      } else if ( LDVSL.LT.1 .OR. ( ILVSL .AND. LDVSL.LT.N ) ) {
         INFO = -15
      } else if ( LDVSR.LT.1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) {
         INFO = -17
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.)

      if ( INFO.EQ.0 ) {
         if ( N.GT.0) {
            MINWRK = 2*N
            MAXWRK = N*(1 + ILAENV( 1, 'ZGEQRF', ' ', N, 1, N, 0 ) )
            MAXWRK = MAX( MAXWRK, N*( 1 + ILAENV( 1, 'ZUNMQR', ' ', N, 1, N, -1 ) ) )
            if ( ILVSL ) {
               MAXWRK = MAX( MAXWRK, N*( 1 + ILAENV( 1, 'ZUNGQR', ' ', N, 1, N, -1 ) ) )
            }
            LWRK = MAXWRK
            IF( IJOB.GE.1 ) LWRK = MAX( LWRK, N*N/2 )
         } else {
            MINWRK = 1
            MAXWRK = 1
            LWRK   = 1
         }
         WORK( 1 ) = LWRK
         if ( WANTSN .OR. N.EQ.0 ) {
            LIWMIN = 1
         } else {
            LIWMIN = N + 2
         }
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) {
            INFO = -21
         } else if ( LIWORK.LT.LIWMIN  .AND. .NOT.LQUERY) {
            INFO = -24
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('ZGGESX', -INFO );
         RETURN
      } else if (LQUERY) {
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         SDIM = 0
         RETURN
      }

      // Get machine constants

      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = ZLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      } else if ( ANRM.GT.BIGNUM ) {
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      }
      IF( ILASCL ) CALL ZLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = ZLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      } else if ( BNRM.GT.BIGNUM ) {
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      }
      IF( ILBSCL ) CALL ZLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )

      // Permute the matrix to make it more nearly triangular
      // (Real Workspace: need 6*N)

      ILEFT = 1
      IRIGHT = N + 1
      IRWRK = IRIGHT + N
      zggbal('P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), RWORK( IRWRK ), IERR );

      // Reduce B to triangular form (QR decomposition of B)
      // (Complex Workspace: need N, prefer N*NB)

      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = 1
      IWRK = ITAU + IROWS
      zgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the unitary transformation to matrix A
      // (Complex Workspace: need N, prefer N*NB)

      zunmqr('L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VSL
      // (Complex Workspace: need N, prefer N*NB)

      if ( ILVSL ) {
         zlaset('Full', N, N, CZERO, CONE, VSL, LDVSL );
         if ( IROWS.GT.1 ) {
            zlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL );
         }
         zungqr(IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Initialize VSR

      IF( ILVSR ) CALL ZLASET( 'Full', N, N, CZERO, CONE, VSR, LDVSR )

      // Reduce to generalized Hessenberg form
      // (Workspace: none needed)

      zgghrd(JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, IERR );

      SDIM = 0

      // Perform QZ algorithm, computing Schur vectors if desired
      // (Complex Workspace: need N)
      // (Real Workspace:    need N)

      IWRK = ITAU
      zhgeqz('S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, RWORK( IRWRK ), IERR );
      if ( IERR.NE.0 ) {
         if ( IERR.GT.0 .AND. IERR.LE.N ) {
            INFO = IERR
         } else if ( IERR.GT.N .AND. IERR.LE.2*N ) {
            INFO = IERR - N
         } else {
            INFO = N + 1
         }
         GO TO 40
      }

      // Sort eigenvalues ALPHA/BETA and compute the reciprocal of
      // condition number(s)

      if ( WANTST ) {

         // Undo scaling on eigenvalues before SELCTGing

         IF( ILASCL ) CALL ZLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )          IF( ILBSCL ) CALL ZLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )

         // Select eigenvalues

         for (I = 1; I <= N; I++) { // 10
            BWORK( I ) = SELCTG( ALPHA( I ), BETA( I ) )
         } // 10

         // Reorder eigenvalues, transform Generalized Schur vectors, and
         // compute reciprocal condition numbers
         // (Complex Workspace: If IJOB >= 1, need MAX(1, 2*SDIM*(N-SDIM))
                             // otherwise, need 1 )

         ztgsen(IJOB, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PL, PR, DIF, WORK( IWRK ), LWORK-IWRK+1, IWORK, LIWORK, IERR );

         IF( IJOB.GE.1 ) MAXWRK = MAX( MAXWRK, 2*SDIM*( N-SDIM ) )
         if ( IERR.EQ.-21 ) {

             // not enough complex workspace

            INFO = -21
         } else {
            if ( IJOB.EQ.1 .OR. IJOB.EQ.4 ) {
               RCONDE( 1 ) = PL
               RCONDE( 2 ) = PR
            }
            if ( IJOB.EQ.2 .OR. IJOB.EQ.4 ) {
               RCONDV( 1 ) = DIF( 1 )
               RCONDV( 2 ) = DIF( 2 )
            }
            IF( IERR.EQ.1 ) INFO = N + 3
         }

      }

      // Apply permutation to VSL and VSR
      // (Workspace: none needed)

      IF( ILVSL ) CALL ZGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSL, LDVSL, IERR )

      IF( ILVSR ) CALL ZGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSR, LDVSR, IERR )

      // Undo scaling

      if ( ILASCL ) {
         zlascl('U', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR );
         zlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR );
      }

      if ( ILBSCL ) {
         zlascl('U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR );
         zlascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );
      }

      if ( WANTST ) {

         // Check if reordering is correct

         LASTSL = .TRUE.
         SDIM = 0
         for (I = 1; I <= N; I++) { // 30
            CURSL = SELCTG( ALPHA( I ), BETA( I ) )
            IF( CURSL ) SDIM = SDIM + 1             IF( CURSL .AND. .NOT.LASTSL ) INFO = N + 2
            LASTSL = CURSL
         } // 30

      }

      } // 40

      WORK( 1 ) = MAXWRK
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of ZGGESX

      }
