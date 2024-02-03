      SUBROUTINE SGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO, IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, WORK, LWORK, IWORK, BWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             BALANC, JOBVL, JOBVR, SENSE;
      int                IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N
      REAL               ABNRM, BBNRM
*     ..
*     .. Array Arguments ..
      bool               BWORK( * );
      int                IWORK( * )
      REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), BETA( * ), LSCALE( * ), RCONDE( * ), RCONDV( * ), RSCALE( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      bool               ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY, NOSCL, PAIR, WANTSB, WANTSE, WANTSN, WANTSV;
      String             CHTEMP;
      int                I, ICOLS, IERR, IJOBVL, IJOBVR, IN, IROWS, ITAU, IWRK, IWRK1, J, JC, JR, M, MAXWRK, MINWRK, MM
      REAL               ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, SMLNUM, TEMP
*     ..
*     .. Local Arrays ..
      bool               LDUMMA( 1 );
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEQRF, SGGBAK, SGGBAL, SGGHRD, SHGEQZ, SLACPY, SLASCL, SLASET, SORGQR, SORMQR, STGEVC, STGSNA, XERBLA
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV
      REAL               SLAMCH, SLANGE, SROUNDUP_LWORK
      EXTERNAL           LSAME, ILAENV, SLAMCH, SLANGE, SROUNDUP_LWORK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( LSAME( JOBVL, 'N' ) ) THEN
         IJOBVL = 1
         ILVL = .FALSE.
      ELSE IF( LSAME( JOBVL, 'V' ) ) THEN
         IJOBVL = 2
         ILVL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVL = .FALSE.
      END IF
*
      IF( LSAME( JOBVR, 'N' ) ) THEN
         IJOBVR = 1
         ILVR = .FALSE.
      ELSE IF( LSAME( JOBVR, 'V' ) ) THEN
         IJOBVR = 2
         ILVR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVR = .FALSE.
      END IF
      ILV = ILVL .OR. ILVR
*
      NOSCL  = LSAME( BALANC, 'N' ) .OR. LSAME( BALANC, 'P' )
      WANTSN = LSAME( SENSE, 'N' )
      WANTSE = LSAME( SENSE, 'E' )
      WANTSV = LSAME( SENSE, 'V' )
      WANTSB = LSAME( SENSE, 'B' )
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.( NOSCL .OR. LSAME( BALANC, 'S' ) .OR. LSAME( BALANC, 'B' ) ) ) THEN
         INFO = -1
      ELSE IF( IJOBVL.LE.0 ) THEN
         INFO = -2
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -3
      ELSE IF( .NOT.( WANTSN .OR. WANTSE .OR. WANTSB .OR. WANTSV ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDVL.LT.1 .OR. ( ILVL .AND. LDVL.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) THEN
         INFO = -16
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV. The workspace is
*       computed assuming ILO = 1 and IHI = N, the worst case.)
*
      IF( INFO.EQ.0 ) THEN
         IF( N.EQ.0 ) THEN
            MINWRK = 1
            MAXWRK = 1
         ELSE
            IF( NOSCL .AND. .NOT.ILV ) THEN
               MINWRK = 2*N
            ELSE
               MINWRK = 6*N
            END IF
            IF( WANTSE ) THEN
               MINWRK = 10*N
            ELSE IF( WANTSV .OR. WANTSB ) THEN
               MINWRK = 2*N*( N + 4 ) + 16
            END IF
            MAXWRK = MINWRK
            MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'SGEQRF', ' ', N, 1, N, 0 ) )             MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'SORMQR', ' ', N, 1, N, 0 ) )
            IF( ILVL ) THEN
               MAXWRK = MAX( MAXWRK, N + N*ILAENV( 1, 'SORGQR', ' ', N, 1, N, 0 ) )
            END IF
         END IF
         WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
*
         IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
            INFO = -26
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGGEVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*
*     Get machine constants
*
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
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
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
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
*
*     Permute and/or balance the matrix pair (A,B)
*     (Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise)
*
      CALL SGGBAL( BALANC, N, A, LDA, B, LDB, ILO, IHI, LSCALE, RSCALE, WORK, IERR )
*
*     Compute ABNRM and BBNRM
*
      ABNRM = SLANGE( '1', N, N, A, LDA, WORK( 1 ) )
      IF( ILASCL ) THEN
         WORK( 1 ) = ABNRM
         CALL SLASCL( 'G', 0, 0, ANRMTO, ANRM, 1, 1, WORK( 1 ), 1, IERR )
         ABNRM = WORK( 1 )
      END IF
*
      BBNRM = SLANGE( '1', N, N, B, LDB, WORK( 1 ) )
      IF( ILBSCL ) THEN
         WORK( 1 ) = BBNRM
         CALL SLASCL( 'G', 0, 0, BNRMTO, BNRM, 1, 1, WORK( 1 ), 1, IERR )
         BBNRM = WORK( 1 )
      END IF
*
*     Reduce B to triangular form (QR decomposition of B)
*     (Workspace: need N, prefer N*NB )
*
      IROWS = IHI + 1 - ILO
      IF( ILV .OR. .NOT.WANTSN ) THEN
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      END IF
      ITAU = 1
      IWRK = ITAU + IROWS
      CALL SGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
*
*     Apply the orthogonal transformation to A
*     (Workspace: need N, prefer N*NB)
*
      CALL SORMQR( 'L', 'T', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB, WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ), LWORK+1-IWRK, IERR )
*
*     Initialize VL and/or VR
*     (Workspace: need N, prefer N*NB)
*
      IF( ILVL ) THEN
         CALL SLASET( 'Full', N, N, ZERO, ONE, VL, LDVL )
         IF( IROWS.GT.1 ) THEN
            CALL SLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB, VL( ILO+1, ILO ), LDVL )
         END IF
         CALL SORGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL, WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF
*
      IF( ILVR ) CALL SLASET( 'Full', N, N, ZERO, ONE, VR, LDVR )
*
*     Reduce to generalized Hessenberg form
*     (Workspace: none needed)
*
      IF( ILV .OR. .NOT.WANTSN ) THEN
*
*        Eigenvectors requested -- work on whole matrix.
*
         CALL SGGHRD( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL, LDVL, VR, LDVR, IERR )
      ELSE
         CALL SGGHRD( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA, B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR )
      END IF
*
*     Perform QZ algorithm (Compute eigenvalues, and optionally, the
*     Schur forms and Schur vectors)
*     (Workspace: need N)
*
      IF( ILV .OR. .NOT.WANTSN ) THEN
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      END IF
*
      CALL SHGEQZ( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 130
      END IF
*
*     Compute Eigenvectors and estimate condition numbers if desired
*     (Workspace: STGEVC: need 6*N
*                 STGSNA: need 2*N*(N+2)+16 if SENSE = 'V' or 'B',
*                         need N otherwise )
*
      IF( ILV .OR. .NOT.WANTSN ) THEN
         IF( ILV ) THEN
            IF( ILVL ) THEN
               IF( ILVR ) THEN
                  CHTEMP = 'B'
               ELSE
                  CHTEMP = 'L'
               END IF
            ELSE
               CHTEMP = 'R'
            END IF
*
            CALL STGEVC( CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL, VR, LDVR, N, IN, WORK, IERR )
            IF( IERR.NE.0 ) THEN
               INFO = N + 2
               GO TO 130
            END IF
         END IF
*
         IF( .NOT.WANTSN ) THEN
*
*           compute eigenvectors (STGEVC) and estimate condition
*           numbers (STGSNA). Note that the definition of the condition
*           number is not invariant under transformation (u,v) to
*           (Q*u, Z*v), where (u,v) are eigenvectors of the generalized
*           Schur form (S,T), Q and Z are orthogonal matrices. In order
*           to avoid using extra 2*N*N workspace, we have to recalculate
*           eigenvectors and estimate one condition numbers at a time.
*
            PAIR = .FALSE.
            DO 20 I = 1, N
*
               IF( PAIR ) THEN
                  PAIR = .FALSE.
                  GO TO 20
               END IF
               MM = 1
               IF( I.LT.N ) THEN
                  IF( A( I+1, I ).NE.ZERO ) THEN
                     PAIR = .TRUE.
                     MM = 2
                  END IF
               END IF
*
               DO 10 J = 1, N
                  BWORK( J ) = .FALSE.
   10          CONTINUE
               IF( MM.EQ.1 ) THEN
                  BWORK( I ) = .TRUE.
               ELSE IF( MM.EQ.2 ) THEN
                  BWORK( I ) = .TRUE.
                  BWORK( I+1 ) = .TRUE.
               END IF
*
               IWRK = MM*N + 1
               IWRK1 = IWRK + MM*N
*
*              Compute a pair of left and right eigenvectors.
*              (compute workspace: need up to 4*N + 6*N)
*
               IF( WANTSE .OR. WANTSB ) THEN
                  CALL STGEVC( 'B', 'S', BWORK, N, A, LDA, B, LDB, WORK( 1 ), N, WORK( IWRK ), N, MM, M, WORK( IWRK1 ), IERR )
                  IF( IERR.NE.0 ) THEN
                     INFO = N + 2
                     GO TO 130
                  END IF
               END IF
*
               CALL STGSNA( SENSE, 'S', BWORK, N, A, LDA, B, LDB, WORK( 1 ), N, WORK( IWRK ), N, RCONDE( I ), RCONDV( I ), MM, M, WORK( IWRK1 ), LWORK-IWRK1+1, IWORK, IERR )
*
   20       CONTINUE
         END IF
      END IF
*
*     Undo balancing on VL and VR and normalization
*     (Workspace: none needed)
*
      IF( ILVL ) THEN
         CALL SGGBAK( BALANC, 'L', N, ILO, IHI, LSCALE, RSCALE, N, VL, LDVL, IERR )
*
         DO 70 JC = 1, N
            IF( ALPHAI( JC ).LT.ZERO ) GO TO 70
            TEMP = ZERO
            IF( ALPHAI( JC ).EQ.ZERO ) THEN
               DO 30 JR = 1, N
                  TEMP = MAX( TEMP, ABS( VL( JR, JC ) ) )
   30          CONTINUE
            ELSE
               DO 40 JR = 1, N
                  TEMP = MAX( TEMP, ABS( VL( JR, JC ) )+ ABS( VL( JR, JC+1 ) ) )
   40          CONTINUE
            END IF
            IF( TEMP.LT.SMLNUM ) GO TO 70
            TEMP = ONE / TEMP
            IF( ALPHAI( JC ).EQ.ZERO ) THEN
               DO 50 JR = 1, N
                  VL( JR, JC ) = VL( JR, JC )*TEMP
   50          CONTINUE
            ELSE
               DO 60 JR = 1, N
                  VL( JR, JC ) = VL( JR, JC )*TEMP
                  VL( JR, JC+1 ) = VL( JR, JC+1 )*TEMP
   60          CONTINUE
            END IF
   70    CONTINUE
      END IF
      IF( ILVR ) THEN
         CALL SGGBAK( BALANC, 'R', N, ILO, IHI, LSCALE, RSCALE, N, VR, LDVR, IERR )
         DO 120 JC = 1, N
            IF( ALPHAI( JC ).LT.ZERO ) GO TO 120
            TEMP = ZERO
            IF( ALPHAI( JC ).EQ.ZERO ) THEN
               DO 80 JR = 1, N
                  TEMP = MAX( TEMP, ABS( VR( JR, JC ) ) )
   80          CONTINUE
            ELSE
               DO 90 JR = 1, N
                  TEMP = MAX( TEMP, ABS( VR( JR, JC ) )+ ABS( VR( JR, JC+1 ) ) )
   90          CONTINUE
            END IF
            IF( TEMP.LT.SMLNUM ) GO TO 120
            TEMP = ONE / TEMP
            IF( ALPHAI( JC ).EQ.ZERO ) THEN
               DO 100 JR = 1, N
                  VR( JR, JC ) = VR( JR, JC )*TEMP
  100          CONTINUE
            ELSE
               DO 110 JR = 1, N
                  VR( JR, JC ) = VR( JR, JC )*TEMP
                  VR( JR, JC+1 ) = VR( JR, JC+1 )*TEMP
  110          CONTINUE
            END IF
  120    CONTINUE
      END IF
*
*     Undo scaling if necessary
*
  130 CONTINUE
*
      IF( ILASCL ) THEN
         CALL SLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAR, N, IERR )
         CALL SLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHAI, N, IERR )
      END IF
*
      IF( ILBSCL ) THEN
         CALL SLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
      END IF
*
      WORK( 1 ) = SROUNDUP_LWORK(MAXWRK)
      RETURN
*
*     End of SGGEVX
*
      END
