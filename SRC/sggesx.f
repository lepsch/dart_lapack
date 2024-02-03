      SUBROUTINE SGGESX( JOBVSL, JOBVSR, SORT, SELCTG, SENSE, N, A, LDA, B, LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, RCONDE, RCONDV, WORK, LWORK, IWORK, LIWORK, BWORK, INFO )

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
      REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), RCONDE( 2 ), RCONDV( 2 ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * )
      // ..
      // .. Function Arguments ..
      bool               SELCTG;
      // EXTERNAL SELCTG
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL, LQUERY, LST2SL, WANTSB, WANTSE, WANTSN, WANTST, WANTSV;
      int                I, ICOLS, IERR, IHI, IJOB, IJOBVL, IJOBVR, ILEFT, ILO, IP, IRIGHT, IROWS, ITAU, IWRK, LIWMIN, LWRK, MAXWRK, MINWRK;
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PL, PR, SAFMAX, SAFMIN, SMLNUM
      // ..
      // .. Local Arrays ..
      REAL               DIF( 2 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQRF, SGGBAK, SGGBAL, SGGHRD, SHGEQZ, SLACPY, SLASCL, SLASET, SORGQR, SORMQR, STGSEN, XERBLA
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
         INFO = -16
      } else if ( LDVSR.LT.1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) {
         INFO = -18
      }

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.)

      if ( INFO.EQ.0 ) {
         if ( N.GT.0) {
            MINWRK = MAX( 8*N, 6*N + 16 )
            MAXWRK = MINWRK - N + N*ILAENV( 1, 'SGEQRF', ' ', N, 1, N, 0 )             MAXWRK = MAX( MAXWRK, MINWRK - N + N*ILAENV( 1, 'SORMQR', ' ', N, 1, N, -1 ) )
            if ( ILVSL ) {
               MAXWRK = MAX( MAXWRK, MINWRK - N + N*ILAENV( 1, 'SORGQR', ' ', N, 1, N, -1 ) )
            }
            LWRK = MAXWRK
            IF( IJOB.GE.1 ) LWRK = MAX( LWRK, N*N/2 )
         } else {
            MINWRK = 1
            MAXWRK = 1
            LWRK   = 1
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWRK)
         if ( WANTSN .OR. N.EQ.0 ) {
            LIWMIN = 1
         } else {
            LIWMIN = N + 6
         }
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) {
            INFO = -22
         } else if ( LIWORK.LT.LIWMIN  .AND. .NOT.LQUERY ) {
            INFO = -24
         }
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SGGESX', -INFO )
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

      EPS = SLAMCH( 'P' )
      SAFMIN = SLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SMLNUM = SQRT( SAFMIN ) / EPS
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

      // Permute the matrix to make it more nearly triangular
      // (Workspace: need 6*N + 2*N for permutation parameters)

      ILEFT = 1
      IRIGHT = N + 1
      IWRK = IRIGHT + N
      CALL SGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWRK ), IERR )

      // Reduce B to triangular form (QR decomposition of B)
      // (Workspace: need N, prefer N*NB)

      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = IWRK
      IWRK = ITAU + IROWS
      CALL SGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )

      // Apply the orthogonal transformation to matrix A
      // (Workspace: need N, prefer N*NB)

      CALL SORMQR( 'L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR )

      // Initialize VSL
      // (Workspace: need N, prefer N*NB)

      if ( ILVSL ) {
         CALL SLASET( 'Full', N, N, ZERO, ONE, VSL, LDVSL )
         if ( IROWS.GT.1 ) {
            CALL SLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL )
         }
         CALL SORGQR( IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      }

      // Initialize VSR

      IF( ILVSR ) CALL SLASET( 'Full', N, N, ZERO, ONE, VSR, LDVSR )

      // Reduce to generalized Hessenberg form
      // (Workspace: none needed)

      CALL SGGHRD( JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, IERR )

      SDIM = 0

      // Perform QZ algorithm, computing Schur vectors if desired
      // (Workspace: need N)

      IWRK = ITAU
      CALL SHGEQZ( 'S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, IERR )
      if ( IERR.NE.0 ) {
         if ( IERR.GT.0 .AND. IERR.LE.N ) {
            INFO = IERR
         } else if ( IERR.GT.N .AND. IERR.LE.2*N ) {
            INFO = IERR - N
         } else {
            INFO = N + 1
         }
         GO TO 50
      }

      // Sort eigenvalues ALPHA/BETA and compute the reciprocal of
      // condition number(s)
      // (Workspace: If IJOB >= 1, need MAX( 8*(N+1), 2*SDIM*(N-SDIM) )
                  // otherwise, need 8*(N+1) )

      if ( WANTST ) {

         // Undo scaling on eigenvalues before SELCTGing

         if ( ILASCL ) {
            CALL SLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR )             CALL SLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR )
         }
         IF( ILBSCL ) CALL SLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )

         // Select eigenvalues

         DO 10 I = 1, N
            BWORK( I ) = SELCTG( ALPHAR( I ), ALPHAI( I ), BETA( I ) )
   10    CONTINUE

         // Reorder eigenvalues, transform Generalized Schur vectors, and
         // compute reciprocal condition numbers

         CALL STGSEN( IJOB, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PL, PR, DIF, WORK( IWRK ), LWORK-IWRK+1, IWORK, LIWORK, IERR )

         IF( IJOB.GE.1 ) MAXWRK = MAX( MAXWRK, 2*SDIM*( N-SDIM ) )
         if ( IERR.EQ.-22 ) {

             // not enough real workspace

            INFO = -22
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

      IF( ILVSL ) CALL SGGBAK( 'P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSL, LDVSL, IERR )

      IF( ILVSR ) CALL SGGBAK( 'P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSR, LDVSR, IERR )

      // Check if unscaling would cause over/underflow, if so, rescale
      // (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of
      // B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I)

      if ( ILASCL ) {
         DO 20 I = 1, N
            if ( ALPHAI( I ).NE.ZERO ) {
               if ( ( ALPHAR( I ) / SAFMAX ).GT.( ANRMTO / ANRM ) .OR. ( SAFMIN / ALPHAR( I ) ).GT.( ANRM / ANRMTO ) ) {
                  WORK( 1 ) = ABS( A( I, I ) / ALPHAR( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               } else if ( ( ALPHAI( I ) / SAFMAX ).GT.( ANRMTO / ANRM ) .OR. ( SAFMIN / ALPHAI( I ) ).GT.( ANRM / ANRMTO ) ) {
                  WORK( 1 ) = ABS( A( I, I+1 ) / ALPHAI( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               }
            }
   20    CONTINUE
      }

      if ( ILBSCL ) {
         DO 25 I = 1, N
            if ( ALPHAI( I ).NE.ZERO ) {
               if ( ( BETA( I ) / SAFMAX ).GT.( BNRMTO / BNRM ) .OR. ( SAFMIN / BETA( I ) ).GT.( BNRM / BNRMTO ) ) {
                  WORK( 1 ) = ABS( B( I, I ) / BETA( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               }
            }
   25    CONTINUE
      }

      // Undo scaling

      if ( ILASCL ) {
         CALL SLASCL( 'H', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR )
         CALL SLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR )
         CALL SLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR )
      }

      if ( ILBSCL ) {
         CALL SLASCL( 'U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR )
         CALL SLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
      }

      if ( WANTST ) {

         // Check if reordering is correct

         LASTSL = .TRUE.
         LST2SL = .TRUE.
         SDIM = 0
         IP = 0
         DO 40 I = 1, N
            CURSL = SELCTG( ALPHAR( I ), ALPHAI( I ), BETA( I ) )
            if ( ALPHAI( I ).EQ.ZERO ) {
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
   40    CONTINUE

      }

   50 CONTINUE

      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of SGGESX

      }
