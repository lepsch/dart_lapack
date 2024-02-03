      SUBROUTINE SGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO, IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, BWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             BALANC, JOBVL, JOBVR, SENSE;
      int                IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N;
      REAL               ABNRM, BBNRM
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      int                IWORK( * );
      REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), LSCALE( * ), RCONDE( * ), RCONDV( * ), RSCALE( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY, NOSCL, PAIR, WANTSB, WANTSE, WANTSN, WANTSV;
      String             CHTEMP;
      int                I, ICOLS, IERR, IJOBVL, IJOBVR, IN, IROWS, ITAU, IWRK, IWRK1, J, JC, JR, M, MAXWRK, MINWRK, MM;
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SMLNUM, TEMP
      // ..
      // .. Local Arrays ..
      bool               LDUMMA( 1 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQRF, SGGBAK, SGGBAL, SGGHRD, SHGEQZ, SLACPY, SLASCL, SLASET, SORGQR, SORMQR, STGEVC, STGSNA, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANGE, SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV, SLAMCH, SLANGE, SROUNDUP_LWORK
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

      NOSCL  = LSAME( BALANC, 'N' ) .OR. LSAME( BALANC, 'P' )
      WANTSN = LSAME( SENSE, 'N' )
      WANTSE = LSAME( SENSE, 'E' )
      WANTSV = LSAME( SENSE, 'V' )
      WANTSB = LSAME( SENSE, 'B' )

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      if ( .NOT.( NOSCL .OR. LSAME( BALANC, 'S' ) .OR. LSAME( BALANC, 'B' ) ) ) {
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
         INFO = -14
      } else if ( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) {
         INFO = -16
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
            if ( NOSCL .AND. .NOT.ILV ) {
               MINWRK = 2*N
            } else {
               MINWRK = 6*N
            }
            if ( WANTSE ) {
               MINWRK = 10*N
            } else if ( WANTSV .OR. WANTSB ) {
               MINWRK = 2*N*( N + 4 ) + 16
            }
            MAXWRK = MINWRK
            MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'SGEQRF', ' ', N, 1, N, 0 ) )             MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'SORMQR', ' ', N, 1, N, 0 ) )
            if ( ILVL ) {
               MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'SORGQR', ' ', N, 1, N, 0 ) )
            }
         }
         WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)

         if ( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) {
            INFO = -26
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('SGGEVX', -INFO );
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

      ANRM = SLANGE( 'M', N, N, A, LDA, WORK )
      ILASCL = .FALSE.
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      } else if ( ANRM.GT.BIGNUM ) {
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      }
      IF( ILASCL ) CALL SLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = SLANGE( 'M', N, N, B, LDB, WORK )
      ILBSCL = .FALSE.
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      } else if ( BNRM.GT.BIGNUM ) {
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      }
      IF( ILBSCL ) CALL SLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )

      // Permute and/or balance the matrix pair (A,B)
      // (Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise)

      sggbal(BALANC, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, IERR );

      // Compute ABNRM and BBNRM

      ABNRM = SLANGE( '1', N, N, A, LDA, WORK( 1 ) )
      if ( ILASCL ) {
         WORK( 1 ) = ABNRM
         slascl('G', 0, 0, ANRMTO, ANRM, 1, 1, WORK( 1 ), 1, IERR );
         ABNRM = WORK( 1 )
      }

      BBNRM = SLANGE( '1', N, N, B, LDB, WORK( 1 ) )
      if ( ILBSCL ) {
         WORK( 1 ) = BBNRM
         slascl('G', 0, 0, BNRMTO, BNRM, 1, 1, WORK( 1 ), 1, IERR );
         BBNRM = WORK( 1 )
      }

      // Reduce B to triangular form (QR decomposition of B)
      // (Workspace: need N, prefer N*NB )

      IROWS = IHI + 1 - ILO
      if ( ILV .OR. .NOT.WANTSN ) {
         ICOLS = N + 1 - ILO
      } else {
         ICOLS = IROWS
      }
      ITAU = 1
      IWRK = ITAU + IROWS
      sgeqrf(IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Apply the orthogonal transformation to A
      // (Workspace: need N, prefer N*NB)

      sormqr('L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR );

      // Initialize VL and/or VR
      // (Workspace: need N, prefer N*NB)

      if ( ILVL ) {
         slaset('Full', N, N, ZERO, ONE, VL, LDVL );
         if ( IROWS.GT.1 ) {
            slacpy('L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL );
         }
         sorgqr(IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR );
      }

      IF( ILVR ) CALL SLASET( 'Full', N, N, ZERO, ONE, VR, LDVR )

      // Reduce to generalized Hessenberg form
      // (Workspace: none needed)

      if ( ILV .OR. .NOT.WANTSN ) {

         // Eigenvectors requested -- work on whole matrix.

         sgghrd(JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, IERR );
      } else {
         sgghrd('N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR );
      }

      // Perform QZ algorithm (Compute eigenvalues, and optionally, the
      // Schur forms and Schur vectors)
      // (Workspace: need N)

      if ( ILV .OR. .NOT.WANTSN ) {
         CHTEMP = 'S'
      } else {
         CHTEMP = 'E'
      }

      shgeqz(CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, IERR );
      if ( IERR.NE.0 ) {
         if ( IERR.GT.0 .AND. IERR.LE.N ) {
            INFO = IERR
         } else if ( IERR.GT.N .AND. IERR.LE.2*N ) {
            INFO = IERR - N
         } else {
            INFO = N + 1
         }
         GO TO 130
      }

      // Compute Eigenvectors and estimate condition numbers if desired
      // (Workspace: STGEVC: need 6*N
                  // STGSNA: need 2*N*(N+2)+16 if SENSE = 'V' or 'B',
                          // need N otherwise )

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

            stgevc(CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK, IERR );
            if ( IERR.NE.0 ) {
               INFO = N + 2
               GO TO 130
            }
         }

         if ( .NOT.WANTSN ) {

            // compute eigenvectors (STGEVC) and estimate condition
            // numbers (STGSNA). Note that the definition of the condition
            // number is not invariant under transformation (u,v) to
            // (Q*u, Z*v), where (u,v) are eigenvectors of the generalized
            // Schur form (S,T), Q and Z are orthogonal matrices. In order
            // to avoid using extra 2*N*N workspace, we have to recalculate
            // eigenvectors and estimate one condition numbers at a time.

            PAIR = .FALSE.
            DO 20 I = 1, N

               if ( PAIR ) {
                  PAIR = .FALSE.
                  GO TO 20
               }
               MM = 1
               if ( I.LT.N ) {
                  if ( A( I+1, I ).NE.ZERO ) {
                     PAIR = .TRUE.
                     MM = 2
                  }
               }

               DO 10 J = 1, N
                  BWORK( J ) = .FALSE.
   10          CONTINUE
               if ( MM.EQ.1 ) {
                  BWORK( I ) = .TRUE.
               } else if ( MM.EQ.2 ) {
                  BWORK( I ) = .TRUE.
                  BWORK( I+1 ) = .TRUE.
               }

               IWRK = MM*N + 1
               IWRK1 = IWRK + MM*N

               // Compute a pair of left and right eigenvectors.
               // (compute workspace: need up to 4*N + 6*N)

               if ( WANTSE .OR. WANTSB ) {
                  stgevc('B', 'S', BWORK, N, A, LDA, B, LDB, WORK( 1 ), N, WORK( IWRK ), N, MM, M, WORK( IWRK1 ), IERR );
                  if ( IERR.NE.0 ) {
                     INFO = N + 2
                     GO TO 130
                  }
               }

               stgsna(SENSE, 'S', BWORK, N, A, LDA, B, LDB, WORK( 1 ), N, WORK( IWRK ), N, RCONDE( I ), RCONDV( I ), MM, M, WORK( IWRK1 ), LWORK-IWRK1+1, IWORK, IERR );

   20       CONTINUE
         }
      }

      // Undo balancing on VL and VR and normalization
      // (Workspace: none needed)

      if ( ILVL ) {
         sggbak(BALANC, 'L', N, ILO, IHI, LSCALE, RSCALE, N, VL, LDVL, IERR );

         DO 70 JC = 1, N
            IF( ALPHAI( JC ).LT.ZERO ) GO TO 70
            TEMP = ZERO
            if ( ALPHAI( JC ).EQ.ZERO ) {
               DO 30 JR = 1, N
                  TEMP = MAX( TEMP, ABS( VL( JR, JC ) ) )
   30          CONTINUE
            } else {
               DO 40 JR = 1, N
                  TEMP = MAX( TEMP, ABS( VL( JR, JC ) )+ ABS( VL( JR, JC+1 ) ) )
   40          CONTINUE
            }
            IF( TEMP.LT.SMLNUM ) GO TO 70
            TEMP = ONE / TEMP
            if ( ALPHAI( JC ).EQ.ZERO ) {
               DO 50 JR = 1, N
                  VL( JR, JC ) = VL( JR, JC )*TEMP
   50          CONTINUE
            } else {
               DO 60 JR = 1, N
                  VL( JR, JC ) = VL( JR, JC )*TEMP
                  VL( JR, JC+1 ) = VL( JR, JC+1 )*TEMP
   60          CONTINUE
            }
   70    CONTINUE
      }
      if ( ILVR ) {
         sggbak(BALANC, 'R', N, ILO, IHI, LSCALE, RSCALE, N, VR, LDVR, IERR );
         DO 120 JC = 1, N
            IF( ALPHAI( JC ).LT.ZERO ) GO TO 120
            TEMP = ZERO
            if ( ALPHAI( JC ).EQ.ZERO ) {
               DO 80 JR = 1, N
                  TEMP = MAX( TEMP, ABS( VR( JR, JC ) ) )
   80          CONTINUE
            } else {
               DO 90 JR = 1, N
                  TEMP = MAX( TEMP, ABS( VR( JR, JC ) )+ ABS( VR( JR, JC+1 ) ) )
   90          CONTINUE
            }
            IF( TEMP.LT.SMLNUM ) GO TO 120
            TEMP = ONE / TEMP
            if ( ALPHAI( JC ).EQ.ZERO ) {
               DO 100 JR = 1, N
                  VR( JR, JC ) = VR( JR, JC )*TEMP
  100          CONTINUE
            } else {
               DO 110 JR = 1, N
                  VR( JR, JC ) = VR( JR, JC )*TEMP
                  VR( JR, JC+1 ) = VR( JR, JC+1 )*TEMP
  110          CONTINUE
            }
  120    CONTINUE
      }

      // Undo scaling if necessary

  130 CONTINUE

      if ( ILASCL ) {
         slascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR );
         slascl('G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR );
      }

      if ( ILBSCL ) {
         slascl('G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR );
      }

      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
      RETURN

      // End of SGGEVX

      }
