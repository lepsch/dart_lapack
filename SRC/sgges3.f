      SUBROUTINE SGGES3( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, BWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVSL, JOBVSR, SORT;
      int                INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * )
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
      bool               CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL, LQUERY, LST2SL, WANTST       int                I, ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO, IP, IRIGHT, IROWS, ITAU, IWRK, LWKOPT, LWKMIN;;
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PVSL, PVSR, SAFMAX, SAFMIN, SMLNUM
      // ..
      // .. Local Arrays ..
      int                IDUM( 1 );
      REAL               DIF( 2 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQRF, SGGBAK, SGGBAL, SGGHD3, SLAQZ0, SLACPY, SLASCL, SLASET, SORGQR, SORMQR, STGSEN, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANGE, SROUNDUP_LWORK
      // EXTERNAL LSAME, SLAMCH, SLANGE, SROUNDUP_LWORK
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

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.EQ.0 ) THEN
         LWKMIN = 1
      ELSE
         LWKMIN = 6*N+16
      END IF

      IF( IJOBVL.LE.0 ) THEN
         INFO = -1
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -2
      ELSE IF( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDVSL.LT.1 .OR. ( ILVSL .AND. LDVSL.LT.N ) ) THEN
         INFO = -15
      ELSE IF( LDVSR.LT.1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) THEN
         INFO = -17
      ELSE IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -19
      END IF

      // Compute workspace

      IF( INFO.EQ.0 ) THEN
         CALL SGEQRF( N, N, B, LDB, WORK, WORK, -1, IERR )
         LWKOPT = MAX( LWKMIN, 3*N+INT( WORK( 1 ) ) )
         CALL SORMQR( 'L', 'T', N, N, N, B, LDB, WORK, A, LDA, WORK, -1, IERR )
         LWKOPT = MAX( LWKOPT, 3*N+INT( WORK( 1 ) ) )
         IF( ILVSL ) THEN
            CALL SORGQR( N, N, N, VSL, LDVSL, WORK, WORK, -1, IERR )
            LWKOPT = MAX( LWKOPT, 3*N+INT( WORK( 1 ) ) )
         END IF
         CALL SGGHD3( JOBVSL, JOBVSR, N, 1, N, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, WORK, -1, IERR )
         LWKOPT = MAX( LWKOPT, 3*N+INT( WORK( 1 ) ) )
         CALL SLAQZ0( 'S', JOBVSL, JOBVSR, N, 1, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK, -1, 0, IERR )
         LWKOPT = MAX( LWKOPT, 2*N+INT( WORK( 1 ) ) )
         IF( WANTST ) THEN
            CALL STGSEN( 0, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PVSL, PVSR, DIF, WORK, -1, IDUM, 1, IERR )
            LWKOPT = MAX( LWKOPT, 2*N+INT( WORK( 1 ) ) )
         END IF
         IF( N.EQ.0 ) THEN
            WORK( 1 ) = 1
         ELSE
            WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
         END IF
      END IF

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGGES3 ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) THEN
         SDIM = 0
         WORK( 1 ) = 1
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

      ILEFT = 1
      IRIGHT = N + 1
      IWRK = IRIGHT + N
      CALL SGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), WORK( IWRK ), IERR )

      // Reduce B to triangular form (QR decomposition of B)

      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = IWRK
      IWRK = ITAU + IROWS
      CALL SGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )

      // Apply the orthogonal transformation to matrix A

      CALL SORMQR( 'L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR )

      // Initialize VSL

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

      CALL SGGHD3( JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, IERR )

      // Perform QZ algorithm, computing Schur vectors if desired

      IWRK = ITAU
      CALL SLAQZ0( 'S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, 0, IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 40
      END IF

      // Sort eigenvalues ALPHA/BETA if desired

      SDIM = 0
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

         CALL STGSEN( 0, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PVSL, PVSR, DIF, WORK( IWRK ), LWORK-IWRK+1, IDUM, 1, IERR )
         IF( IERR.EQ.1 ) INFO = N + 3

      END IF

      // Apply back-permutation to VSL and VSR

      IF( ILVSL ) CALL SGGBAK( 'P', 'L', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSL, LDVSL, IERR )

      IF( ILVSR ) CALL SGGBAK( 'P', 'R', N, ILO, IHI, WORK( ILEFT ), WORK( IRIGHT ), N, VSR, LDVSR, IERR )

      // Check if unscaling would cause over/underflow, if so, rescale
      // (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of
      // B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I)

      IF( ILASCL )THEN
         DO 50 I = 1, N
            IF( ALPHAI( I ).NE.ZERO ) THEN
               IF( ( ALPHAR( I )/SAFMAX ).GT.( ANRMTO/ANRM ) .OR. ( SAFMIN/ALPHAR( I ) ).GT.( ANRM/ANRMTO ) ) THEN
                  WORK( 1 ) = ABS( A( I, I )/ALPHAR( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               ELSE IF( ( ALPHAI( I )/SAFMAX ).GT.( ANRMTO/ANRM ) .OR. ( SAFMIN/ALPHAI( I ) ).GT.( ANRM/ANRMTO ) ) THEN
                  WORK( 1 ) = ABS( A( I, I+1 )/ALPHAI( I ) )
                  BETA( I ) = BETA( I )*WORK( 1 )
                  ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                  ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
               END IF
            END IF
   50    CONTINUE
      END IF

      IF( ILBSCL )THEN
         DO 60 I = 1, N
            IF( ALPHAI( I ).NE.ZERO ) THEN
                IF( ( BETA( I )/SAFMAX ).GT.( BNRMTO/BNRM ) .OR. ( SAFMIN/BETA( I ) ).GT.( BNRM/BNRMTO ) ) THEN
                   WORK( 1 ) = ABS(B( I, I )/BETA( I ))
                   BETA( I ) = BETA( I )*WORK( 1 )
                   ALPHAR( I ) = ALPHAR( I )*WORK( 1 )
                   ALPHAI( I ) = ALPHAI( I )*WORK( 1 )
                END IF
             END IF
   60    CONTINUE
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
         DO 30 I = 1, N
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
   30    CONTINUE

      END IF

   40 CONTINUE

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

      RETURN

      // End of SGGES3

      }
