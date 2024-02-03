      SUBROUTINE ZGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, ILO, IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, IWORK, BWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             BALANC, JOBVL, JOBVR, SENSE;
      int                IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      double             ABNRM, BBNRM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      int                IWORK( * );
      double             LSCALE( * ), RCONDE( * ), RCONDV( * ), RSCALE( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY, NOSCL, WANTSB, WANTSE, WANTSN, WANTSV;
      String             CHTEMP;
      int                I, ICOLS, IERR, IJOBVL, IJOBVR, IN, IROWS, ITAU, IWRK, IWRK1, J, JC, JR, M, MAXWRK, MINWRK;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SMLNUM, TEMP;
      COMPLEX*16         X
      // ..
      // .. Local Arrays ..
      bool               LDUMMA( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASCL, XERBLA, ZGEQRF, ZGGBAK, ZGGBAL, ZGGHRD, ZHGEQZ, ZLACPY, ZLASCL, ZLASET, ZTGEVC, ZTGSNA, ZUNGQR, ZUNMQR
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, ILAENV, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, SQRT
      // ..
      // .. Statement Functions ..
      double             ABS1;
      // ..
      // .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
      // ..
      // .. Executable Statements ..

      // Decode the input arguments

      if ( LSAME( JOBVL, 'N' ) ) {
         IJOBVL = 1
         ILVL = .FALSE.
      } else if ( LSAME( JOBVL, 'V' ) ) {
         IJOBVL = 2
         ILVL = .TRUE.
      } else {
         IJOBVL = -1
         ILVL = .FALSE.
      }

      if ( LSAME( JOBVR, 'N' ) ) {
         IJOBVR = 1
         ILVR = .FALSE.
      } else if ( LSAME( JOBVR, 'V' ) ) {
         IJOBVR = 2
         ILVR = .TRUE.
      } else {
         IJOBVR = -1
         ILVR = .FALSE.
      }
      ILV = ILVL .OR. ILVR

      NOSCL  = LSAME( BALANC, 'N' ) .OR. LSAME( BALANC, 'P' )
      WANTSN = LSAME( SENSE, 'N' )
      WANTSE = LSAME( SENSE, 'E' )
      WANTSV = LSAME( SENSE, 'V' )
      WANTSB = LSAME( SENSE, 'B' )

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      if ( .NOT.( NOSCL .OR. LSAME( BALANC,'S' ) .OR. LSAME( BALANC, 'B' ) ) ) {
         INFO = -1
      } else if ( IJOBVL.LE.0 ) {
         INFO = -2
      } else if ( IJOBVR.LE.0 ) {
         INFO = -3
      } else if ( .NOT.( WANTSN .OR. WANTSE .OR. WANTSB .OR. WANTSV ) ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -9
      } else if ( LDVL.LT.1 .OR. ( ILVL .AND. LDVL.LT.N ) ) {
         INFO = -13
      } else if ( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) {
         INFO = -15
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV. The workspace is
        // computed assuming ILO = 1 and IHI = N, the worst case.)

      if ( INFO.EQ.0 ) {
         if ( N.EQ.0 ) {
            MINWRK = 1
            MAXWRK = 1
         } else {
            MINWRK = 2*N
            if ( WANTSE ) {
               MINWRK = 4*N
            } else if ( WANTSV .OR. WANTSB ) {
               MINWRK = 2*N*( N + 1)
            }
            MAXWRK = MINWRK
            MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'ZGEQRF', ' ', N, 1, N, 0 ) )             MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'ZUNMQR', ' ', N, 1, N, 0 ) )
            if ( ILVL ) {
               MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'ZUNGQR', ' ', N, 1, N, 0 ) )
            }
         }
         WORK( 1 ) = MAXWRK

         if ( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) {
            INFO = -25
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('ZGGEVX', -INFO );
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

      // Permute and/or balance the matrix pair (A,B)
      // (Real Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise)

      zggbal(BALANC, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, RWORK, IERR );

      // Compute ABNRM and BBNRM

      ABNRM = ZLANGE( '1', N, N, A, LDA, RWORK( 1 ) )
      if ( ILASCL ) {
         RWORK( 1 ) = ABNRM
         dlascl('G', 0, 0, ANRMTO, ANRM, 1, 1, RWORK( 1 ), 1, IERR );
         ABNRM = RWORK( 1 )
      }

      BBNRM = ZLANGE( '1', N, N, B, LDB, RWORK( 1 ) )
      if ( ILBSCL ) {
         RWORK( 1 ) = BBNRM
         dlascl('G', 0, 0, BNRMTO, BNRM, 1, 1, RWORK( 1 ), 1, IERR );
         BBNRM = RWORK( 1 )
      }

      // Reduce B to triangular form (QR decomposition of B)
      // (Complex Workspace: need N, prefer N*NB )

      IROWS = IHI + 1 - ILO
      if ( ILV .OR. .NOT.WANTSN ) {
         ICOLS = N + 1 - ILO
      } else {
         ICOLS = IROWS
      }
      ITAU = 1
      IWRK = ITAU + IROWS
      zgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the unitary transformation to A
      // (Complex Workspace: need N, prefer N*NB)

      zunmqr('L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VL and/or VR
      // (Workspace: need N, prefer N*NB)

      if ( ILVL ) {
         zlaset('Full', N, N, CZERO, CONE, VL, LDVL );
         if ( IROWS.GT.1 ) {
            zlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         }
         zungqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      IF( ILVR ) CALL ZLASET( 'Full', N, N, CZERO, CONE, VR, LDVR )

      // Reduce to generalized Hessenberg form
      // (Workspace: none needed)

      if ( ILV .OR. .NOT.WANTSN ) {

         // Eigenvectors requested -- work on whole matrix.

         zgghrd(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, IERR );
      } else {
         zgghrd('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR );
      }

      // Perform QZ algorithm (Compute eigenvalues, and optionally, the
      // Schur forms and Schur vectors)
      // (Complex Workspace: need N)
      // (Real Workspace: need N)

      IWRK = ITAU
      if ( ILV .OR. .NOT.WANTSN ) {
         CHTEMP = 'S'
      } else {
         CHTEMP = 'E'
      }

      zhgeqz(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, RWORK, IERR );
      if ( IERR.NE.0 ) {
         if ( IERR.GT.0 .AND. IERR.LE.N ) {
            INFO = IERR
         } else if ( IERR.GT.N .AND. IERR.LE.2*N ) {
            INFO = IERR - N
         } else {
            INFO = N + 1
         }
         GO TO 90
      }

      // Compute Eigenvectors and estimate condition numbers if desired
      // ZTGEVC: (Complex Workspace: need 2*N )
              // (Real Workspace:    need 2*N )
      // ZTGSNA: (Complex Workspace: need 2*N*N if SENSE='V' or 'B')
              // (Integer Workspace: need N+2 )

      if ( ILV .OR. .NOT.WANTSN ) {
         if ( ILV ) {
            if ( ILVL ) {
               if ( ILVR ) {
                  CHTEMP = 'B'
               } else {
                  CHTEMP = 'L'
               }
            } else {
               CHTEMP = 'R'
            }

            ztgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWRK ), RWORK, IERR );
            if ( IERR.NE.0 ) {
               INFO = N + 2
               GO TO 90
            }
         }

         if ( .NOT.WANTSN ) {

            // compute eigenvectors (ZTGEVC) and estimate condition
            // numbers (ZTGSNA). Note that the definition of the condition
            // number is not invariant under transformation (u,v) to
            // (Q*u, Z*v), where (u,v) are eigenvectors of the generalized
            // Schur form (S,T), Q and Z are orthogonal matrices. In order
            // to avoid using extra 2*N*N workspace, we have to
            // re-calculate eigenvectors and estimate the condition numbers
            // one at a time.

            for (I = 1; I <= N; I++) { // 20

               for (J = 1; J <= N; J++) { // 10
                  BWORK( J ) = .FALSE.
   10          CONTINUE
               BWORK( I ) = .TRUE.

               IWRK = N + 1
               IWRK1 = IWRK + N

               if ( WANTSE .OR. WANTSB ) {
                  ztgevc('B', 'S', BWORK, N, A, LDA, B, LDB, WORK( 1 ), N, WORK( IWRK ), N, 1, M, WORK( IWRK1 ), RWORK, IERR );
                  if ( IERR.NE.0 ) {
                     INFO = N + 2
                     GO TO 90
                  }
               }

               ztgsna(SENSE, 'S', BWORK, N, A, LDA, B, LDB, WORK( 1 ), N, WORK( IWRK ), N, RCONDE( I ), RCONDV( I ), 1, M, WORK( IWRK1 ), LWORK-IWRK1+1, IWORK, IERR );

   20       CONTINUE
         }
      }

      // Undo balancing on VL and VR and normalization
      // (Workspace: none needed)

      if ( ILVL ) {
         zggbak(BALANC, 'L', N, ILO, IHI, LSCALE, RSCALE, N, VL, LDVL, IERR );

         for (JC = 1; JC <= N; JC++) { // 50
            TEMP = ZERO
            for (JR = 1; JR <= N; JR++) { // 30
               TEMP = MAX( TEMP, ABS1( VL( JR, JC ) ) )
   30       CONTINUE
            IF( TEMP.LT.SMLNUM ) GO TO 50
            TEMP = ONE / TEMP
            for (JR = 1; JR <= N; JR++) { // 40
               VL( JR, JC ) = VL( JR, JC )*TEMP
   40       CONTINUE
   50    CONTINUE
      }

      if ( ILVR ) {
         zggbak(BALANC, 'R', N, ILO, IHI, LSCALE, RSCALE, N, VR, LDVR, IERR );
         for (JC = 1; JC <= N; JC++) { // 80
            TEMP = ZERO
            for (JR = 1; JR <= N; JR++) { // 60
               TEMP = MAX( TEMP, ABS1( VR( JR, JC ) ) )
   60       CONTINUE
            IF( TEMP.LT.SMLNUM ) GO TO 80
            TEMP = ONE / TEMP
            for (JR = 1; JR <= N; JR++) { // 70
               VR( JR, JC ) = VR( JR, JC )*TEMP
   70       CONTINUE
   80    CONTINUE
      }

      // Undo scaling if necessary

   90 CONTINUE

      IF( ILASCL ) CALL ZLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )

      IF( ILBSCL ) CALL ZLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )

      WORK( 1 ) = MAXWRK
      RETURN

      // End of ZGGEVX

      }
