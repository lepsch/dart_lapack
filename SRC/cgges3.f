      SUBROUTINE CGGES3( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK, BWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBVSL, JOBVSR, SORT;
      int                INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM;
      // ..
      // .. Array Arguments ..
      bool               BWORK( * );
      REAL               RWORK( * )
      COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ), WORK( * )
      // ..
      // .. Function Arguments ..
      bool               SELCTG;
      // EXTERNAL SELCTG
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E0, 0.0E0 ), CONE = ( 1.0E0, 0.0E0 ) ;
      // ..
      // .. Local Scalars ..
      bool               CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL, LQUERY, WANTST;
      int                I, ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO, IRIGHT, IROWS, IRWRK, ITAU, IWRK, LWKOPT, LWKMIN;
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PVSL, PVSR, SMLNUM
      // ..
      // .. Local Arrays ..
      int                IDUM( 1 );
      REAL               DIF( 2 )
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQRF, CGGBAK, CGGBAL, CGGHD3, CLAQZ0, CLACPY, CLASCL, CLASET, CTGSEN, CUNGQR, CUNMQR, XERBLA
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGE, SLAMCH, SROUNDUP_LWORK
      // EXTERNAL LSAME, CLANGE, SLAMCH, SROUNDUP_LWORK
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

      // Test the input arguments

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      LWKMIN = MAX( 1, 2*N )

      if ( IJOBVL.LE.0 ) {
         INFO = -1
      } else if ( IJOBVR.LE.0 ) {
         INFO = -2
      } else if ( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -9
      } else if ( LDVSL.LT.1 .OR. ( ILVSL .AND. LDVSL.LT.N ) ) {
         INFO = -14
      } else if ( LDVSR.LT.1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) {
         INFO = -16
      } else if ( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) {
         INFO = -18
      }

      // Compute workspace

      if ( INFO.EQ.0 ) {
         CALL CGEQRF( N, N, B, LDB, WORK, WORK, -1, IERR )
         LWKOPT = MAX( LWKMIN, N + INT( WORK( 1 ) ) )
         CALL CUNMQR( 'L', 'C', N, N, N, B, LDB, WORK, A, LDA, WORK, -1, IERR )
         LWKOPT = MAX( LWKOPT, N + INT( WORK( 1 ) ) )
         if ( ILVSL ) {
            CALL CUNGQR( N, N, N, VSL, LDVSL, WORK, WORK, -1, IERR )
            LWKOPT = MAX( LWKOPT, N + INT( WORK( 1 ) ) )
         }
         CALL CGGHD3( JOBVSL, JOBVSR, N, 1, N, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, WORK, -1, IERR )
         LWKOPT = MAX( LWKOPT, N + INT( WORK( 1 ) ) )
         CALL CLAQZ0( 'S', JOBVSL, JOBVSR, N, 1, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK, -1, RWORK, 0, IERR )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
         if ( WANTST ) {
            CALL CTGSEN( 0, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PVSL, PVSR, DIF, WORK, -1, IDUM, 1, IERR )
            LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
         }
         if ( N.EQ.0 ) {
            WORK( 1 ) = 1
         } else {
            WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
         }
      }


      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CGGES3 ', -INFO )
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         SDIM = 0
         RETURN
      }

      // Get machine constants

      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM

      // Scale A if max element outside range [SMLNUM,BIGNUM]

      ANRM = CLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      if ( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) {
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      } else if ( ANRM.GT.BIGNUM ) {
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      }

      IF( ILASCL ) CALL CLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )

      // Scale B if max element outside range [SMLNUM,BIGNUM]

      BNRM = CLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      if ( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) {
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      } else if ( BNRM.GT.BIGNUM ) {
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      }

      IF( ILBSCL ) CALL CLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )

      // Permute the matrix to make it more nearly triangular

      ILEFT = 1
      IRIGHT = N + 1
      IRWRK = IRIGHT + N
      CALL CGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), RWORK( IRWRK ), IERR )

      // Reduce B to triangular form (QR decomposition of B)

      IROWS = IHI + 1 - ILO
      ICOLS = N + 1 - ILO
      ITAU = 1
      IWRK = ITAU + IROWS
      CALL CGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )

      // Apply the orthogonal transformation to matrix A

      CALL CUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR )

      // Initialize VSL

      if ( ILVSL ) {
         CALL CLASET( 'Full', N, N, CZERO, CONE, VSL, LDVSL )
         if ( IROWS.GT.1 ) {
            CALL CLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VSL( ILO+1, ILO ), LDVSL )
         }
         CALL CUNGQR( IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      }

      // Initialize VSR

      IF( ILVSR ) CALL CLASET( 'Full', N, N, CZERO, CONE, VSR, LDVSR )

      // Reduce to generalized Hessenberg form

      CALL CGGHD3( JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, IERR )

      SDIM = 0

      // Perform QZ algorithm, computing Schur vectors if desired

      IWRK = ITAU
      CALL CLAQZ0( 'S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ), LWORK+1-IWRK, RWORK( IRWRK ), 0, IERR )
      if ( IERR.NE.0 ) {
         if ( IERR.GT.0 .AND. IERR.LE.N ) {
            INFO = IERR
         } else if ( IERR.GT.N .AND. IERR.LE.2*N ) {
            INFO = IERR - N
         } else {
            INFO = N + 1
         }
         GO TO 30
      }

      // Sort eigenvalues ALPHA/BETA if desired

      if ( WANTST ) {

         // Undo scaling on eigenvalues before selecting

         IF( ILASCL ) CALL CLASCL( 'G', 0, 0, ANRM, ANRMTO, N, 1, ALPHA, N, IERR )          IF( ILBSCL ) CALL CLASCL( 'G', 0, 0, BNRM, BNRMTO, N, 1, BETA, N, IERR )

         // Select eigenvalues

         DO 10 I = 1, N
            BWORK( I ) = SELCTG( ALPHA( I ), BETA( I ) )
   10    CONTINUE

         CALL CTGSEN( 0, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PVSL, PVSR, DIF, WORK( IWRK ), LWORK-IWRK+1, IDUM, 1, IERR )
         IF( IERR.EQ.1 ) INFO = N + 3

      }

      // Apply back-permutation to VSL and VSR

      IF( ILVSL ) CALL CGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSL, LDVSL, IERR )       IF( ILVSR ) CALL CGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ), RWORK( IRIGHT ), N, VSR, LDVSR, IERR )

      // Undo scaling

      if ( ILASCL ) {
         CALL CLASCL( 'U', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR )
         CALL CLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )
      }

      if ( ILBSCL ) {
         CALL CLASCL( 'U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR )
         CALL CLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
      }

      if ( WANTST ) {

         // Check if reordering is correct

         LASTSL = .TRUE.
         SDIM = 0
         DO 20 I = 1, N
            CURSL = SELCTG( ALPHA( I ), BETA( I ) )
            IF( CURSL ) SDIM = SDIM + 1             IF( CURSL .AND. .NOT.LASTSL ) INFO = N + 2
            LASTSL = CURSL
   20    CONTINUE

      }

   30 CONTINUE

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

      RETURN

      // End of CGGES3

      }
