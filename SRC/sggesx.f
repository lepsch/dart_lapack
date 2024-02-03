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
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL, LQUERY, LST2SL, WANTSB, WANTSE, WANTSN, WANTST, WANTSV       int                I, ICOLS, IERR, IHI, IJOB, IJOBVL, IJOBVR, ILEFT, ILO, IP, IRIGHT, IROWS, ITAU, IWRK, LIWMIN, LWRK, MAXWRK, MINWRK;;
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

      IF( LSAME( JOBVSL, 'N' ) ) THEN
         IJOBVL = 1
         ILVSL = .FALSE.
      ELSE IF( LSAME( JOBVSL, 'V' ) ) THEN
         IJOBVL = 2
         ILVSL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVSL = .FALSE.
      END IF

      IF( LSAME( JOBVSR, 'N' ) ) THEN
         IJOBVR = 1
         ILVSR = .FALSE.
      ELSE IF( LSAME( JOBVSR, 'V' ) ) THEN
         IJOBVR = 2
         ILVSR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVSR = .FALSE.
      END IF

      WANTST = LSAME( SORT, 'S' )
      WANTSN = LSAME( SENSE, 'N' )
      WANTSE = LSAME( SENSE, 'E' )
      WANTSV = LSAME( SENSE, 'V' )
      WANTSB = LSAME( SENSE, 'B' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
      IF( WANTSN ) THEN
         IJOB = 0
      ELSE IF( WANTSE ) THEN
         IJOB = 1
      ELSE IF( WANTSV ) THEN
         IJOB = 2
      ELSE IF( WANTSB ) THEN
         IJOB = 4
      END IF

      // Test the input arguments

      INFO = 0
      IF( IJOBVL.LE.0 ) THEN
         INFO = -1
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -2
      ELSE IF( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( .NOT.( WANTSN .OR. WANTSE .OR. WANTSV .OR. WANTSB ) .OR. ( .NOT.WANTST .AND. .NOT.WANTSN ) ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDVSL.LT.1 .OR. ( ILVSL .AND. LDVSL.LT.N ) ) THEN
         INFO = -16
      ELSE IF( LDVSR.LT.1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) THEN
         INFO = -18
      END IF

      // Compute workspace
       // (Note: Comments in the code beginning "Workspace:" describe the
        // minimal amount of workspace needed at that point in the code,
        // as well as the preferred amount for good performance.
        // NB refers to the optimal block size for the immediately
        // following subroutine, as returned by ILAENV.)

      IF( INFO.EQ.0 ) THEN
         IF( N.GT.0) THEN
            MINWRK = MAX( 8*N, 6*N + 16 )
            MAXWRK = MINWRK - N + N*ILAENV( 1, 'SGEQRF', ' ', N, 1, N, 0 )             MAXWRK = MAX( MAXWRK, MINWRK - N + N*ILAENV( 1, 'SORMQR', ' ', N, 1, N, -1 ) )
            IF( ILVSL ) THEN
               MAXWRK = MAX( MAXWRK, MINWRK - N + N*ILAENV( 1, 'SORGQR', ' ', N, 1, N, -1 ) )
            END IF
            LWRK = MAXWRK
            IF( IJOB.GE.1 ) LWRK = MAX( LWRK, N*N/2 )
         ELSE
            MINWRK = 1
            MAXWRK = 1
            LWRK   = 1
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(LWRK)
         IF( WANTSN .OR. N.EQ.0 ) THEN
            LIWMIN = 1
         ELSE
            LIWMIN = N + 6
         END IF
         IWORK( 1 ) = LIWMIN

         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -22
         ELSE IF( LIWORK.LT.LIWMIN  .AND. .NOT.LQUERY ) THEN
            INFO = -24
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGGESX', -INFO )
         RETURN
      ELSE IF (LQUERY) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) THEN
         SDIM = 0
         RETURN
      END IF

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SAFMIN = SLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SMLNUM = SQRT( SAFMIN ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = SLANGE( 'M', N, N, A, LDA, WORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
      IF( ILASCL ) CALL SLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = SLANGE( 'M', N, N, B, LDB, WORK )
      ILBSCL = .FALSE.
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      ELSE IF( BNRM.GT.BIGNUM ) THEN
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      END IF
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

      IF( ILVSL ) THEN
         CALL SLASET( 'Full', N, N, ZERO, ONE, VSL, LDVSL )
         IF( IROWS.GT.1 ) THEN
            CALL SLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL )
         END IF
         CALL SORGQR( IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF

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
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 50
      END IF

      // Sort eigenvalues ALPHA/BETA and compute the reciprocal of
      // condition number(s)
      // (Workspace: If IJOB >= 1, need MAX( 8*(N+1), 2*SDIM*(N-SDIM) )
                  // otherwise, need 8*(N+1) )

      IF( WANTST ) THEN

         // Undo scaling on eigenvalues before SELCTGing

         IF( ILASCL ) THEN
            CALL SLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR )             CALL SLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR )
         END IF
         IF( ILBSCL ) CALL SLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )

         // Select eigenvalues

         DO 10 I = 1, N
            BWORK( I ) = SELCTG( ALPHAR( I ), ALPHAI( I ), BETA( I ) )
   10    CONTINUE

         // Reorder eigenvalues, transform Generalized Schur vectors, and
         // compute reciprocal condition numbers

         CALL STGSEN( IJOB, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PL, PR, DIF, WORK( IWRK ), LWORK-IWRK+1, IWORK, LIWORK, IERR )

         IF( IJOB.GE.1 ) MAXWRK = MAX( MAXWRK, 2*SDIM*( N-SDIM ) )
         IF( IERR.EQ.-22 ) THEN

             // not enough real workspace

            INFO = -22
         ELSE
            IF( IJOB.EQ.1 .OR. IJOB.EQ.4 ) THEN
               RCONDE( 1 ) = PL
               RCONDE( 2 ) = PR
            END IF
            IF( IJOB.EQ.2 .OR. IJOB.EQ.4 ) THEN
               RCONDV( 1 ) = DIF( 1 )
               RCONDV( 2 ) = DIF( 2 )
            END IF
            IF( IERR.EQ.1 ) INFO = N + 3
         END IF

      END IF

      // Apply permutation to VSL and VSR
      // (Workspace: none needed)

      IF( ILVSL ) CALL SGGBAK( 'P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSL, LDVSL, IERR )

      IF( ILVSR ) CALL SGGBAK( 'P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSR, LDVSR, IERR )

      // Check if unscaling would cause over/underflow, if so, rescale
      // (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of
      // B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I)

      IF( ILASCL ) THEN
         DO 20 I = 1, N
            IF( ALPHAI( I ).NE.ZERO ) THEN
               IF( ( ALPHAR( I ) / SAFMAX ).GT.( ANRMTO / ANRM ) .OR. ( SAFMIN / ALPHAR( I ) ).GT.( ANRM / ANRMTO ) ) THEN
                  WORK( 1 ) = ABS( A( I, I ) / ALPHAR( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               ELSE IF( ( ALPHAI( I ) / SAFMAX ).GT.( ANRMTO / ANRM ) .OR. ( SAFMIN / ALPHAI( I ) ).GT.( ANRM / ANRMTO ) ) THEN
                  WORK( 1 ) = ABS( A( I, I+1 ) / ALPHAI( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               END IF
            END IF
   20    CONTINUE
      END IF

      IF( ILBSCL ) THEN
         DO 25 I = 1, N
            IF( ALPHAI( I ).NE.ZERO ) THEN
               IF( ( BETA( I ) / SAFMAX ).GT.( BNRMTO / BNRM ) .OR. ( SAFMIN / BETA( I ) ).GT.( BNRM / BNRMTO ) ) THEN
                  WORK( 1 ) = ABS( B( I, I ) / BETA( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               END IF
            END IF
   25    CONTINUE
      END IF

      // Undo scaling

      IF( ILASCL ) THEN
         CALL SLASCL( 'H', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR )
         CALL SLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR )
         CALL SLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR )
      END IF

      IF( ILBSCL ) THEN
         CALL SLASCL( 'U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR )
         CALL SLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
      END IF

      IF( WANTST ) THEN

         // Check if reordering is correct

         LASTSL = .TRUE.
         LST2SL = .TRUE.
         SDIM = 0
         IP = 0
         DO 40 I = 1, N
            CURSL = SELCTG( ALPHAR( I ), ALPHAI( I ), BETA( I ) )
            IF( ALPHAI( I ).EQ.ZERO ) THEN
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
   40    CONTINUE

      END IF

   50 CONTINUE

      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of SGGESX

      }
