      SUBROUTINE ZGGES3( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK, BWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             JOBVSL, JOBVSR, SORT;
      int                INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * )
      // ..
      // .. Function Arguments ..
      bool               SELCTG;
      // EXTERNAL SELCTG
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ), CONE = ( 1.0D0, 0.0D0 ) )
      // ..
      // .. Local Scalars ..
      bool               CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL, LQUERY, WANTST       int                I, ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO, IRIGHT, IROWS, IRWRK, ITAU, IWRK, LWKOPT, LWKMIN;;
      double             ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PVSL, PVSR, SMLNUM;
      // ..
      // .. Local Arrays ..
      int                IDUM( 1 );
      double             DIF( 2 );
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQRF, ZGGBAK, ZGGBAL, ZGGHD3, ZLAQZ0, ZLACPY, ZLASCL, ZLASET, ZTGSEN, ZUNGQR, ZUNMQR
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, DLAMCH, ZLANGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..
*
      // Decode the input arguments
*
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
*
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
*
      WANTST = LSAME( SORT, 'S' )
*
      // Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      LWKMIN = MAX( 1, 2*N )
*
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
         INFO = -14
      ELSE IF( LDVSR.LT.1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) THEN
         INFO = -16
      ELSE IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -18
      END IF
*
      // Compute workspace
*
      IF( INFO.EQ.0 ) THEN
         CALL ZGEQRF( N, N, B, LDB, WORK, WORK, -1, IERR )
         LWKOPT = MAX( LWKMIN,  N + INT( WORK( 1 ) ) )
         CALL ZUNMQR( 'L', 'C', N, N, N, B, LDB, WORK, A, LDA, WORK, -1, IERR )
         LWKOPT = MAX( LWKOPT, N + INT( WORK( 1 ) ) )
         IF( ILVSL ) THEN
            CALL ZUNGQR( N, N, N, VSL, LDVSL, WORK, WORK, -1, IERR )
            LWKOPT = MAX( LWKOPT, N + INT ( WORK( 1 ) ) )
         END IF
         CALL ZGGHD3( JOBVSL, JOBVSR, N, 1, N, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, WORK, -1, IERR )
         LWKOPT = MAX( LWKOPT, N + INT( WORK( 1 ) ) )
         CALL ZLAQZ0( 'S', JOBVSL, JOBVSR, N, 1, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, -1, RWORK, 0, IERR )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
         IF( WANTST ) THEN
            CALL ZTGSEN( 0, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PVSL, PVSR, DIF, WORK, -1, IDUM, 1, IERR )
            LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
         END IF
         IF( N.EQ.0 ) THEN
            WORK( 1 ) = 1
         ELSE
            WORK( 1 ) = DCMPLX( LWKOPT )
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGES3 ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 ) THEN
         SDIM = 0
         RETURN
      END IF
*
      // Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
      // Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = ZLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
*
      IF( ILASCL ) CALL ZLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )
*
      // Scale B if max element outside range [SMLNUM,BIGNUM]
*
      BNRM = ZLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      ELSE IF( BNRM.GT.BIGNUM ) THEN
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      END IF
*
      IF( ILBSCL ) CALL ZLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )
*
      // Permute the matrix to make it more nearly triangular
*
      ILEFT = 1
      IRIGHT = N + 1
      IRWRK = IRIGHT + N
      CALL ZGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), RWORK( IRWRK ), IERR )
*
      // Reduce B to triangular form (QR decomposition of B)
*
      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = 1
      IWRK = ITAU + IROWS
      CALL ZGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
*
      // Apply the orthogonal transformation to matrix A
*
      CALL ZUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR )
*
      // Initialize VSL
*
      IF( ILVSL ) THEN
         CALL ZLASET( 'Full', N, N, CZERO, CONE, VSL, LDVSL )
         IF( IROWS.GT.1 ) THEN
            CALL ZLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL )
         END IF
         CALL ZUNGQR( IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF
*
      // Initialize VSR
*
      IF( ILVSR ) CALL ZLASET( 'Full', N, N, CZERO, CONE, VSR, LDVSR )
*
      // Reduce to generalized Hessenberg form
*
      CALL ZGGHD3( JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, IERR )
*
      SDIM = 0
*
      // Perform QZ algorithm, computing Schur vectors if desired
*
      IWRK = ITAU
      CALL ZLAQZ0( 'S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, RWORK( IRWRK ), 0, IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 30
      END IF
*
      // Sort eigenvalues ALPHA/BETA if desired
*
      IF( WANTST ) THEN
*
         // Undo scaling on eigenvalues before selecting
*
         IF( ILASCL ) CALL ZLASCL( 'G', 0, 0, ANRM, ANRMTO, N, 1, ALPHA, N, IERR )          IF( ILBSCL ) CALL ZLASCL( 'G', 0, 0, BNRM, BNRMTO, N, 1, BETA, N, IERR )
*
         // Select eigenvalues
*
         DO 10 I = 1, N
            BWORK( I ) = SELCTG( ALPHA( I ), BETA( I ) )
   10    CONTINUE
*
         CALL ZTGSEN( 0, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PVSL, PVSR, DIF, WORK( IWRK ), LWORK-IWRK+1, IDUM, 1, IERR )
         IF( IERR.EQ.1 ) INFO = N + 3
*
      END IF
*
      // Apply back-permutation to VSL and VSR
*
      IF( ILVSL ) CALL ZGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSL, LDVSL, IERR )       IF( ILVSR ) CALL ZGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSR, LDVSR, IERR )
*
      // Undo scaling
*
      IF( ILASCL ) THEN
         CALL ZLASCL( 'U', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR )
         CALL ZLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )
      END IF
*
      IF( ILBSCL ) THEN
         CALL ZLASCL( 'U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR )
         CALL ZLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
      END IF
*
      IF( WANTST ) THEN
*
         // Check if reordering is correct
*
         LASTSL = .TRUE.
         SDIM = 0
         DO 20 I = 1, N
            CURSL = SELCTG( ALPHA( I ), BETA( I ) )
            IF( CURSL ) SDIM = SDIM + 1             IF( CURSL .AND. .NOT.LASTSL ) INFO = N + 2
            LASTSL = CURSL
   20    CONTINUE
*
      END IF
*
   30 CONTINUE
*
      WORK( 1 ) = DCMPLX( LWKOPT )
*
      RETURN
*
      // End of ZGGES3
*
      END
