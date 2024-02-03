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
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
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

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      LWKMIN = MAX( 1, 8*N )
      if ( IJOBVL.LE.0 ) {
         INFO = -1
      } else if ( IJOBVR.LE.0 ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDVL.LT.1 .OR. ( ILVL .AND. LDVL.LT.N ) ) {
         INFO = -12
      } else if ( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) {
         INFO = -14
      } else if ( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) {
         INFO = -16
      }

      // Compute workspace

      if ( INFO.EQ.0 ) {
         dgeqrf(N, N, B, LDB, WORK, WORK, -1, IERR );
         LWKOPT = MAX( LWKMIN, 3*N+INT( WORK( 1 ) ) )
         dormqr('L', 'T', N, N, N, B, LDB, WORK, A, LDA, WORK, -1, IERR );
         LWKOPT = MAX( LWKOPT, 3*N+INT( WORK( 1 ) ) )
         if ( ILVL ) {
            dorgqr(N, N, N, VL, LDVL, WORK, WORK, -1, IERR );
            LWKOPT = MAX( LWKOPT, 3*N+INT( WORK( 1 ) ) )
         }
         if ( ILV ) {
            dgghd3(JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK, -1, IERR );
            LWKOPT = MAX( LWKOPT, 3*N+INT( WORK ( 1 ) ) )
            dlaqz0('S', JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, -1, 0, IERR );
            LWKOPT = MAX( LWKOPT, 2*N+INT( WORK( 1 ) ) )
         } else {
            dgghd3('N', 'N', N, 1, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK, -1, IERR );
            LWKOPT = MAX( LWKOPT, 3*N+INT( WORK( 1 ) ) )
            dlaqz0('E', JOBVL, JOBVR, N, 1, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, -1, 0, IERR );
            LWKOPT = MAX( LWKOPT, 2*N+INT( WORK( 1 ) ) )
         }
         if ( N.EQ.0 ) {
            WORK( 1 ) = 1
         } else {
            WORK( 1 ) = LWKOPT
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('DGGEV3 ', -INFO );
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

      ANRM = DLANGE( 'M', N, N, A, LDA, WORK )
      ILASCL = .FALSE.
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      } else if ( ANRM.GT.BIGNUM ) {
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      }
      IF( ILASCL ) CALL DLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = DLANGE( 'M', N, N, B, LDB, WORK )
      ILBSCL = .FALSE.
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      } else if ( BNRM.GT.BIGNUM ) {
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      }
      IF( ILBSCL ) CALL DLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )

      // Permute the matrices A, B to isolate eigenvalues if possible

      ILEFT = 1
      IRIGHT = N + 1
      IWRK = IRIGHT + N
      dggbal('P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWRK ), IERR );

      // Reduce B to triangular form (QR decomposition of B)

      IROWS = IHI + 1 - ILO
      if ( ILV ) {
         ICOLS = N + 1 - ILO
      } else {
         ICOLS = IROWS
      }
      ITAU = IWRK
      IWRK = ITAU + IROWS
      dgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the orthogonal transformation to matrix A

      dormqr('L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VL

      if ( ILVL ) {
         dlaset('Full', N, N, ZERO, ONE, VL, LDVL );
         if ( IROWS.GT.1 ) {
            dlacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         }
         dorgqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Initialize VR

      IF( ILVR ) CALL DLASET( 'Full', N, N, ZERO, ONE, VR, LDVR )

      // Reduce to generalized Hessenberg form

      if ( ILV ) {

         // Eigenvectors requested -- work on whole matrix.

         dgghd3(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, IERR );
      } else {
         dgghd3('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      // Perform QZ algorithm (Compute eigenvalues, and optionally, the
      // Schur forms and Schur vectors)

      IWRK = ITAU
      if ( ILV ) {
         CHTEMP = 'S'
      } else {
         CHTEMP = 'E'
      }
      dlaqz0(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ), LWORK+1-IWRK, 0, IERR );
      if ( IERR.NE.0 ) {
         if ( IERR.GT.0 .AND. IERR.LE.N ) {
            INFO = IERR
         } else if ( IERR.GT.N .AND. IERR.LE.2*N ) {
            INFO = IERR - N
         } else {
            INFO = N + 1
         }
         GO TO 110
      }

      // Compute Eigenvectors

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
         dtgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK( IWRK ), IERR );
         if ( IERR.NE.0 ) {
            INFO = N + 2
            GO TO 110
         }

         // Undo balancing on VL and VR and normalization

         if ( ILVL ) {
            dggbak('P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VL, LDVL, IERR );
            for (JC = 1; JC <= N; JC++) { // 50
               IF( ALPHAI( JC ).LT.ZERO ) GO TO 50
               TEMP = ZERO
               if ( ALPHAI( JC ).EQ.ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 10
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) ) )
                  } // 10
               } else {
                  for (JR = 1; JR <= N; JR++) { // 20
                     TEMP = MAX( TEMP, ABS( VL( JR, JC ) )+ ABS( VL( JR, JC+1 ) ) )
                  } // 20
               }
               IF( TEMP.LT.SMLNUM ) GO TO 50
               TEMP = ONE / TEMP
               if ( ALPHAI( JC ).EQ.ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 30
                     VL( JR, JC ) = VL( JR, JC )*TEMP
                  } // 30
               } else {
                  for (JR = 1; JR <= N; JR++) { // 40
                     VL( JR, JC ) = VL( JR, JC )*TEMP
                     VL( JR, JC+1 ) = VL( JR, JC+1 )*TEMP
                  } // 40
               }
            } // 50
         }
         if ( ILVR ) {
            dggbak('P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VR, LDVR, IERR );
            for (JC = 1; JC <= N; JC++) { // 100
               IF( ALPHAI( JC ).LT.ZERO ) GO TO 100
               TEMP = ZERO
               if ( ALPHAI( JC ).EQ.ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 60
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) ) )
                  } // 60
               } else {
                  for (JR = 1; JR <= N; JR++) { // 70
                     TEMP = MAX( TEMP, ABS( VR( JR, JC ) )+ ABS( VR( JR, JC+1 ) ) )
                  } // 70
               }
               IF( TEMP.LT.SMLNUM ) GO TO 100
               TEMP = ONE / TEMP
               if ( ALPHAI( JC ).EQ.ZERO ) {
                  for (JR = 1; JR <= N; JR++) { // 80
                     VR( JR, JC ) = VR( JR, JC )*TEMP
                  } // 80
               } else {
                  for (JR = 1; JR <= N; JR++) { // 90
                     VR( JR, JC ) = VR( JR, JC )*TEMP
                     VR( JR, JC+1 ) = VR( JR, JC+1 )*TEMP
                  } // 90
               }
            } // 100
         }

         // End of eigenvector calculation

      }

      // Undo scaling if necessary

      } // 110

      if ( ILASCL ) {
         dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR );
         dlascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR );
      }

      if ( ILBSCL ) {
         dlascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );
      }

      WORK( 1 ) = LWKOPT
      RETURN

      // End of DGGEV3

      }
