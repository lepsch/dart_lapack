      SUBROUTINE CLATMT( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RANK, KL, KU, PACK, A, LDA, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL               COND, DMAX
      INTEGER            INFO, KL, KU, LDA, M, MODE, N, RANK
      CHARACTER          DIST, PACK, SYM
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), WORK( * )
      REAL               D( * )
      INTEGER            ISEED( 4 )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      COMPLEX            CZERO
      PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
      REAL               TWOPI
      PARAMETER  ( TWOPI = 6.28318530717958647692528676655900576839E+0 )
*     ..
*     .. Local Scalars ..
      COMPLEX            C, CT, CTEMP, DUMMY, EXTRA, S, ST
      REAL               ALPHA, ANGLE, REALC, TEMP
      INTEGER            I, IC, ICOL, IDIST, IENDCH, IINFO, IL, ILDA, IOFFG, IOFFST, IPACK, IPACKG, IR, IR1, IR2, IROW, IRSIGN, ISKEW, ISYM, ISYMPK, J, JC, JCH, JKL, JKU, JR, K, LLB, MINLDA, MNMIN, MR, NC, UUB
      LOGICAL            CSYM, GIVENS, ILEXTR, ILTEMP, TOPDWN
*     ..
*     .. External Functions ..
      COMPLEX            CLARND
      REAL               SLARND
      LOGICAL            LSAME
      EXTERNAL           CLARND, SLARND, LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLAGGE, CLAGHE, CLAGSY, CLAROT, CLARTG, CLASET, SLATM7, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, CONJG, COS, MAX, MIN, MOD, REAL, SIN
*     ..
*     .. Executable Statements ..
*
*     1)      Decode and Test the input parameters.
*             Initialize flags & seed.
*
      INFO = 0
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) RETURN
*
*     Decode DIST
*
      IF( LSAME( DIST, 'U' ) ) THEN
         IDIST = 1
      ELSE IF( LSAME( DIST, 'S' ) ) THEN
         IDIST = 2
      ELSE IF( LSAME( DIST, 'N' ) ) THEN
         IDIST = 3
      ELSE
         IDIST = -1
      END IF
*
*     Decode SYM
*
      IF( LSAME( SYM, 'N' ) ) THEN
         ISYM = 1
         IRSIGN = 0
         CSYM = .FALSE.
      ELSE IF( LSAME( SYM, 'P' ) ) THEN
         ISYM = 2
         IRSIGN = 0
         CSYM = .FALSE.
      ELSE IF( LSAME( SYM, 'S' ) ) THEN
         ISYM = 2
         IRSIGN = 0
         CSYM = .TRUE.
      ELSE IF( LSAME( SYM, 'H' ) ) THEN
         ISYM = 2
         IRSIGN = 1
         CSYM = .FALSE.
      ELSE
         ISYM = -1
      END IF
*
*     Decode PACK
*
      ISYMPK = 0
      IF( LSAME( PACK, 'N' ) ) THEN
         IPACK = 0
      ELSE IF( LSAME( PACK, 'U' ) ) THEN
         IPACK = 1
         ISYMPK = 1
      ELSE IF( LSAME( PACK, 'L' ) ) THEN
         IPACK = 2
         ISYMPK = 1
      ELSE IF( LSAME( PACK, 'C' ) ) THEN
         IPACK = 3
         ISYMPK = 2
      ELSE IF( LSAME( PACK, 'R' ) ) THEN
         IPACK = 4
         ISYMPK = 3
      ELSE IF( LSAME( PACK, 'B' ) ) THEN
         IPACK = 5
         ISYMPK = 3
      ELSE IF( LSAME( PACK, 'Q' ) ) THEN
         IPACK = 6
         ISYMPK = 2
      ELSE IF( LSAME( PACK, 'Z' ) ) THEN
         IPACK = 7
      ELSE
         IPACK = -1
      END IF
*
*     Set certain internal parameters
*
      MNMIN = MIN( M, N )
      LLB = MIN( KL, M-1 )
      UUB = MIN( KU, N-1 )
      MR = MIN( M, N+LLB )
      NC = MIN( N, M+UUB )
*
      IF( IPACK.EQ.5 .OR. IPACK.EQ.6 ) THEN
         MINLDA = UUB + 1
      ELSE IF( IPACK.EQ.7 ) THEN
         MINLDA = LLB + UUB + 1
      ELSE
         MINLDA = M
      END IF
*
*     Use Givens rotation method if bandwidth small enough,
*     or if LDA is too small to store the matrix unpacked.
*
      GIVENS = .FALSE.
      IF( ISYM.EQ.1 ) THEN
         IF( REAL( LLB+UUB ).LT.0.3*REAL( MAX( 1, MR+NC ) ) ) GIVENS = .TRUE.
      ELSE
         IF( 2*LLB.LT.M ) GIVENS = .TRUE.
      END IF
      IF( LDA.LT.M .AND. LDA.GE.MINLDA ) GIVENS = .TRUE.
*
*     Set INFO if an error
*
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.NE.N .AND. ISYM.NE.1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( IDIST.EQ.-1 ) THEN
         INFO = -3
      ELSE IF( ISYM.EQ.-1 ) THEN
         INFO = -5
      ELSE IF( ABS( MODE ).GT.6 ) THEN
         INFO = -7
      ELSE IF( ( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) .AND. COND.LT.ONE ) THEN
         INFO = -8
      ELSE IF( KL.LT.0 ) THEN
         INFO = -10
      ELSE IF( KU.LT.0 .OR. ( ISYM.NE.1 .AND. KL.NE.KU ) ) THEN
         INFO = -11
      ELSE IF( IPACK.EQ.-1 .OR. ( ISYMPK.EQ.1 .AND. ISYM.EQ.1 ) .OR. ( ISYMPK.EQ.2 .AND. ISYM.EQ.1 .AND. KL.GT.0 ) .OR. ( ISYMPK.EQ.3 .AND. ISYM.EQ.1 .AND. KU.GT.0 ) .OR. ( ISYMPK.NE.0 .AND. M.NE.N ) ) THEN
         INFO = -12
      ELSE IF( LDA.LT.MAX( 1, MINLDA ) ) THEN
         INFO = -14
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLATMT', -INFO )
         RETURN
      END IF
*
*     Initialize random number generator
*
      DO 100 I = 1, 4
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
  100 CONTINUE
*
      IF( MOD( ISEED( 4 ), 2 ).NE.1 ) ISEED( 4 ) = ISEED( 4 ) + 1
*
*     2)      Set up D  if indicated.
*
*             Compute D according to COND and MODE
*
      CALL SLATM7( MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, RANK, IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
*
*     Choose Top-Down if D is (apparently) increasing,
*     Bottom-Up if D is (apparently) decreasing.
*
      IF( ABS( D( 1 ) ).LE.ABS( D( RANK ) ) ) THEN
         TOPDWN = .TRUE.
      ELSE
         TOPDWN = .FALSE.
      END IF
*
      IF( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) THEN
*
*        Scale by DMAX
*
         TEMP = ABS( D( 1 ) )
         DO 110 I = 2, RANK
            TEMP = MAX( TEMP, ABS( D( I ) ) )
  110    CONTINUE
*
         IF( TEMP.GT.ZERO ) THEN
            ALPHA = DMAX / TEMP
         ELSE
            INFO = 2
            RETURN
         END IF
*
         CALL SSCAL( RANK, ALPHA, D, 1 )
*
      END IF
*
      CALL CLASET( 'Full', LDA, N, CZERO, CZERO, A, LDA )
*
*     3)      Generate Banded Matrix using Givens rotations.
*             Also the special case of UUB=LLB=0
*
*               Compute Addressing constants to cover all
*               storage formats.  Whether GE, HE, SY, GB, HB, or SB,
*               upper or lower triangle or both,
*               the (i,j)-th element is in
*               A( i - ISKEW*j + IOFFST, j )
*
      IF( IPACK.GT.4 ) THEN
         ILDA = LDA - 1
         ISKEW = 1
         IF( IPACK.GT.5 ) THEN
            IOFFST = UUB + 1
         ELSE
            IOFFST = 1
         END IF
      ELSE
         ILDA = LDA
         ISKEW = 0
         IOFFST = 0
      END IF
*
*     IPACKG is the format that the matrix is generated in. If this is
*     different from IPACK, then the matrix must be repacked at the
*     end.  It also signals how to compute the norm, for scaling.
*
      IPACKG = 0
*
*     Diagonal Matrix -- We are done, unless it
*     is to be stored HP/SP/PP/TP (PACK='R' or 'C')
*
      IF( LLB.EQ.0 .AND. UUB.EQ.0 ) THEN
         DO 120 J = 1, MNMIN
            A( ( 1-ISKEW )*J+IOFFST, J ) = CMPLX( D( J ) )
  120    CONTINUE
*
         IF( IPACK.LE.2 .OR. IPACK.GE.5 ) IPACKG = IPACK
*
      ELSE IF( GIVENS ) THEN
*
*        Check whether to use Givens rotations,
*        Householder transformations, or nothing.
*
         IF( ISYM.EQ.1 ) THEN
*
*           Non-symmetric -- A = U D V
*
            IF( IPACK.GT.4 ) THEN
               IPACKG = IPACK
            ELSE
               IPACKG = 0
            END IF
*
            DO 130 J = 1, MNMIN
               A( ( 1-ISKEW )*J+IOFFST, J ) = CMPLX( D( J ) )
  130       CONTINUE
*
            IF( TOPDWN ) THEN
               JKL = 0
               DO 160 JKU = 1, UUB
*
*                 Transform from bandwidth JKL, JKU-1 to JKL, JKU
*
*                 Last row actually rotated is M
*                 Last column actually rotated is MIN( M+JKU, N )
*
                  DO 150 JR = 1, MIN( M+JKU, N ) + JKL - 1
                     EXTRA = CZERO
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )*CLARND( 5, ISEED )
                     S = SIN( ANGLE )*CLARND( 5, ISEED )
                     ICOL = MAX( 1, JR-JKL )
                     IF( JR.LT.M ) THEN
                        IL = MIN( N, JR+JKU ) + 1 - ICOL
                        CALL CLAROT( .TRUE., JR.GT.JKL, .FALSE., IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, EXTRA, DUMMY )
                     END IF
*
*                    Chase "EXTRA" back up
*
                     IR = JR
                     IC = ICOL
                     DO 140 JCH = JR - JKL, 1, -JKL - JKU
                        IF( IR.LT.M ) THEN
                           CALL CLARTG( A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, REALC, S, DUMMY )
                           DUMMY = CLARND( 5, ISEED )
                           C = CONJG( REALC*DUMMY )
                           S = CONJG( -S*DUMMY )
                        END IF
                        IROW = MAX( 1, JCH-JKU )
                        IL = IR + 2 - IROW
                        CTEMP = CZERO
                        ILTEMP = JCH.GT.JKU
                        CALL CLAROT( .FALSE., ILTEMP, .TRUE., IL, C, S, A( IROW-ISKEW*IC+IOFFST, IC ), ILDA, CTEMP, EXTRA )
                        IF( ILTEMP ) THEN
                           CALL CLARTG( A( IROW+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), CTEMP, REALC, S, DUMMY )
                           DUMMY = CLARND( 5, ISEED )
                           C = CONJG( REALC*DUMMY )
                           S = CONJG( -S*DUMMY )
*
                           ICOL = MAX( 1, JCH-JKU-JKL )
                           IL = IC + 2 - ICOL
                           EXTRA = CZERO
                           CALL CLAROT( .TRUE., JCH.GT.JKU+JKL, .TRUE., IL, C, S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, CTEMP )
                           IC = ICOL
                           IR = IROW
                        END IF
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
*
               JKU = UUB
               DO 190 JKL = 1, LLB
*
*                 Transform from bandwidth JKL-1, JKU to JKL, JKU
*
                  DO 180 JC = 1, MIN( N+JKL, M ) + JKU - 1
                     EXTRA = CZERO
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )*CLARND( 5, ISEED )
                     S = SIN( ANGLE )*CLARND( 5, ISEED )
                     IROW = MAX( 1, JC-JKU )
                     IF( JC.LT.N ) THEN
                        IL = MIN( M, JC+JKL ) + 1 - IROW
                        CALL CLAROT( .FALSE., JC.GT.JKU, .FALSE., IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, EXTRA, DUMMY )
                     END IF
*
*                    Chase "EXTRA" back up
*
                     IC = JC
                     IR = IROW
                     DO 170 JCH = JC - JKU, 1, -JKL - JKU
                        IF( IC.LT.N ) THEN
                           CALL CLARTG( A( IR+1-ISKEW*( IC+1 )+IOFFST, IC+1 ), EXTRA, REALC, S, DUMMY )
                           DUMMY = CLARND( 5, ISEED )
                           C = CONJG( REALC*DUMMY )
                           S = CONJG( -S*DUMMY )
                        END IF
                        ICOL = MAX( 1, JCH-JKL )
                        IL = IC + 2 - ICOL
                        CTEMP = CZERO
                        ILTEMP = JCH.GT.JKL
                        CALL CLAROT( .TRUE., ILTEMP, .TRUE., IL, C, S, A( IR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, CTEMP, EXTRA )
                        IF( ILTEMP ) THEN
                           CALL CLARTG( A( IR+1-ISKEW*( ICOL+1 )+IOFFST, ICOL+1 ), CTEMP, REALC, S, DUMMY )
                           DUMMY = CLARND( 5, ISEED )
                           C = CONJG( REALC*DUMMY )
                           S = CONJG( -S*DUMMY )
                           IROW = MAX( 1, JCH-JKL-JKU )
                           IL = IR + 2 - IROW
                           EXTRA = CZERO
                           CALL CLAROT( .FALSE., JCH.GT.JKL+JKU, .TRUE., IL, C, S, A( IROW-ISKEW*ICOL+ IOFFST, ICOL ), ILDA, EXTRA, CTEMP )
                           IC = ICOL
                           IR = IROW
                        END IF
  170                CONTINUE
  180             CONTINUE
  190          CONTINUE
*
            ELSE
*
*              Bottom-Up -- Start at the bottom right.
*
               JKL = 0
               DO 220 JKU = 1, UUB
*
*                 Transform from bandwidth JKL, JKU-1 to JKL, JKU
*
*                 First row actually rotated is M
*                 First column actually rotated is MIN( M+JKU, N )
*
                  IENDCH = MIN( M, N+JKL ) - 1
                  DO 210 JC = MIN( M+JKU, N ) - 1, 1 - JKL, -1
                     EXTRA = CZERO
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )*CLARND( 5, ISEED )
                     S = SIN( ANGLE )*CLARND( 5, ISEED )
                     IROW = MAX( 1, JC-JKU+1 )
                     IF( JC.GT.0 ) THEN
                        IL = MIN( M, JC+JKL+1 ) + 1 - IROW
                        CALL CLAROT( .FALSE., .FALSE., JC+JKL.LT.M, IL, C, S, A( IROW-ISKEW*JC+IOFFST, JC ), ILDA, DUMMY, EXTRA )
                     END IF
*
*                    Chase "EXTRA" back down
*
                     IC = JC
                     DO 200 JCH = JC + JKL, IENDCH, JKL + JKU
                        ILEXTR = IC.GT.0
                        IF( ILEXTR ) THEN
                           CALL CLARTG( A( JCH-ISKEW*IC+IOFFST, IC ), EXTRA, REALC, S, DUMMY )
                           DUMMY = CLARND( 5, ISEED )
                           C = REALC*DUMMY
                           S = S*DUMMY
                        END IF
                        IC = MAX( 1, IC )
                        ICOL = MIN( N-1, JCH+JKU )
                        ILTEMP = JCH + JKU.LT.N
                        CTEMP = CZERO
                        CALL CLAROT( .TRUE., ILEXTR, ILTEMP, ICOL+2-IC, C, S, A( JCH-ISKEW*IC+IOFFST, IC ), ILDA, EXTRA, CTEMP )
                        IF( ILTEMP ) THEN
                           CALL CLARTG( A( JCH-ISKEW*ICOL+IOFFST, ICOL ), CTEMP, REALC, S, DUMMY )
                           DUMMY = CLARND( 5, ISEED )
                           C = REALC*DUMMY
                           S = S*DUMMY
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = CZERO
                           CALL CLAROT( .FALSE., .TRUE., JCH+JKL+JKU.LE.IENDCH, IL, C, S, A( JCH-ISKEW*ICOL+IOFFST, ICOL ), ILDA, CTEMP, EXTRA )
                           IC = ICOL
                        END IF
  200                CONTINUE
  210             CONTINUE
  220          CONTINUE
*
               JKU = UUB
               DO 250 JKL = 1, LLB
*
*                 Transform from bandwidth JKL-1, JKU to JKL, JKU
*
*                 First row actually rotated is MIN( N+JKL, M )
*                 First column actually rotated is N
*
                  IENDCH = MIN( N, M+JKU ) - 1
                  DO 240 JR = MIN( N+JKL, M ) - 1, 1 - JKU, -1
                     EXTRA = CZERO
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )*CLARND( 5, ISEED )
                     S = SIN( ANGLE )*CLARND( 5, ISEED )
                     ICOL = MAX( 1, JR-JKL+1 )
                     IF( JR.GT.0 ) THEN
                        IL = MIN( N, JR+JKU+1 ) + 1 - ICOL
                        CALL CLAROT( .TRUE., .FALSE., JR+JKU.LT.N, IL, C, S, A( JR-ISKEW*ICOL+IOFFST, ICOL ), ILDA, DUMMY, EXTRA )
                     END IF
*
*                    Chase "EXTRA" back down
*
                     IR = JR
                     DO 230 JCH = JR + JKU, IENDCH, JKL + JKU
                        ILEXTR = IR.GT.0
                        IF( ILEXTR ) THEN
                           CALL CLARTG( A( IR-ISKEW*JCH+IOFFST, JCH ), EXTRA, REALC, S, DUMMY )
                           DUMMY = CLARND( 5, ISEED )
                           C = REALC*DUMMY
                           S = S*DUMMY
                        END IF
                        IR = MAX( 1, IR )
                        IROW = MIN( M-1, JCH+JKL )
                        ILTEMP = JCH + JKL.LT.M
                        CTEMP = CZERO
                        CALL CLAROT( .FALSE., ILEXTR, ILTEMP, IROW+2-IR, C, S, A( IR-ISKEW*JCH+IOFFST, JCH ), ILDA, EXTRA, CTEMP )
                        IF( ILTEMP ) THEN
                           CALL CLARTG( A( IROW-ISKEW*JCH+IOFFST, JCH ), CTEMP, REALC, S, DUMMY )
                           DUMMY = CLARND( 5, ISEED )
                           C = REALC*DUMMY
                           S = S*DUMMY
                           IL = MIN( IENDCH, JCH+JKL+JKU ) + 2 - JCH
                           EXTRA = CZERO
                           CALL CLAROT( .TRUE., .TRUE., JCH+JKL+JKU.LE.IENDCH, IL, C, S, A( IROW-ISKEW*JCH+IOFFST, JCH ), ILDA, CTEMP, EXTRA )
                           IR = IROW
                        END IF
  230                CONTINUE
  240             CONTINUE
  250          CONTINUE
*
            END IF
*
         ELSE
*
*           Symmetric -- A = U D U'
*           Hermitian -- A = U D U*
*
            IPACKG = IPACK
            IOFFG = IOFFST
*
            IF( TOPDWN ) THEN
*
*              Top-Down -- Generate Upper triangle only
*
               IF( IPACK.GE.5 ) THEN
                  IPACKG = 6
                  IOFFG = UUB + 1
               ELSE
                  IPACKG = 1
               END IF
*
               DO 260 J = 1, MNMIN
                  A( ( 1-ISKEW )*J+IOFFG, J ) = CMPLX( D( J ) )
  260          CONTINUE
*
               DO 290 K = 1, UUB
                  DO 280 JC = 1, N - 1
                     IROW = MAX( 1, JC-K )
                     IL = MIN( JC+1, K+2 )
                     EXTRA = CZERO
                     CTEMP = A( JC-ISKEW*( JC+1 )+IOFFG, JC+1 )
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )*CLARND( 5, ISEED )
                     S = SIN( ANGLE )*CLARND( 5, ISEED )
                     IF( CSYM ) THEN
                        CT = C
                        ST = S
                     ELSE
                        CTEMP = CONJG( CTEMP )
                        CT = CONJG( C )
                        ST = CONJG( S )
                     END IF
                     CALL CLAROT( .FALSE., JC.GT.K, .TRUE., IL, C, S, A( IROW-ISKEW*JC+IOFFG, JC ), ILDA, EXTRA, CTEMP )                      CALL CLAROT( .TRUE., .TRUE., .FALSE., MIN( K, N-JC )+1, CT, ST, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, CTEMP, DUMMY )
*
*                    Chase EXTRA back up the matrix
*
                     ICOL = JC
                     DO 270 JCH = JC - K, 1, -K
                        CALL CLARTG( A( JCH+1-ISKEW*( ICOL+1 )+IOFFG, ICOL+1 ), EXTRA, REALC, S, DUMMY )
                        DUMMY = CLARND( 5, ISEED )
                        C = CONJG( REALC*DUMMY )
                        S = CONJG( -S*DUMMY )
                        CTEMP = A( JCH-ISKEW*( JCH+1 )+IOFFG, JCH+1 )
                        IF( CSYM ) THEN
                           CT = C
                           ST = S
                        ELSE
                           CTEMP = CONJG( CTEMP )
                           CT = CONJG( C )
                           ST = CONJG( S )
                        END IF
                        CALL CLAROT( .TRUE., .TRUE., .TRUE., K+2, C, S, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, CTEMP, EXTRA )
                        IROW = MAX( 1, JCH-K )
                        IL = MIN( JCH+1, K+2 )
                        EXTRA = CZERO
                        CALL CLAROT( .FALSE., JCH.GT.K, .TRUE., IL, CT, ST, A( IROW-ISKEW*JCH+IOFFG, JCH ), ILDA, EXTRA, CTEMP )
                        ICOL = JCH
  270                CONTINUE
  280             CONTINUE
  290          CONTINUE
*
*              If we need lower triangle, copy from upper. Note that
*              the order of copying is chosen to work for 'q' -> 'b'
*
               IF( IPACK.NE.IPACKG .AND. IPACK.NE.3 ) THEN
                  DO 320 JC = 1, N
                     IROW = IOFFST - ISKEW*JC
                     IF( CSYM ) THEN
                        DO 300 JR = JC, MIN( N, JC+UUB )
                           A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
  300                   CONTINUE
                     ELSE
                        DO 310 JR = JC, MIN( N, JC+UUB )
                           A( JR+IROW, JC ) = CONJG( A( JC-ISKEW*JR+ IOFFG, JR ) )
  310                   CONTINUE
                     END IF
  320             CONTINUE
                  IF( IPACK.EQ.5 ) THEN
                     DO 340 JC = N - UUB + 1, N
                        DO 330 JR = N + 2 - JC, UUB + 1
                           A( JR, JC ) = CZERO
  330                   CONTINUE
  340                CONTINUE
                  END IF
                  IF( IPACKG.EQ.6 ) THEN
                     IPACKG = IPACK
                  ELSE
                     IPACKG = 0
                  END IF
               END IF
            ELSE
*
*              Bottom-Up -- Generate Lower triangle only
*
               IF( IPACK.GE.5 ) THEN
                  IPACKG = 5
                  IF( IPACK.EQ.6 ) IOFFG = 1
               ELSE
                  IPACKG = 2
               END IF
*
               DO 350 J = 1, MNMIN
                  A( ( 1-ISKEW )*J+IOFFG, J ) = CMPLX( D( J ) )
  350          CONTINUE
*
               DO 380 K = 1, UUB
                  DO 370 JC = N - 1, 1, -1
                     IL = MIN( N+1-JC, K+2 )
                     EXTRA = CZERO
                     CTEMP = A( 1+( 1-ISKEW )*JC+IOFFG, JC )
                     ANGLE = TWOPI*SLARND( 1, ISEED )
                     C = COS( ANGLE )*CLARND( 5, ISEED )
                     S = SIN( ANGLE )*CLARND( 5, ISEED )
                     IF( CSYM ) THEN
                        CT = C
                        ST = S
                     ELSE
                        CTEMP = CONJG( CTEMP )
                        CT = CONJG( C )
                        ST = CONJG( S )
                     END IF
                     CALL CLAROT( .FALSE., .TRUE., N-JC.GT.K, IL, C, S, A( ( 1-ISKEW )*JC+IOFFG, JC ), ILDA, CTEMP, EXTRA )
                     ICOL = MAX( 1, JC-K+1 )
                     CALL CLAROT( .TRUE., .FALSE., .TRUE., JC+2-ICOL, CT, ST, A( JC-ISKEW*ICOL+IOFFG, ICOL ), ILDA, DUMMY, CTEMP )
*
*                    Chase EXTRA back down the matrix
*
                     ICOL = JC
                     DO 360 JCH = JC + K, N - 1, K
                        CALL CLARTG( A( JCH-ISKEW*ICOL+IOFFG, ICOL ), EXTRA, REALC, S, DUMMY )
                        DUMMY = CLARND( 5, ISEED )
                        C = REALC*DUMMY
                        S = S*DUMMY
                        CTEMP = A( 1+( 1-ISKEW )*JCH+IOFFG, JCH )
                        IF( CSYM ) THEN
                           CT = C
                           ST = S
                        ELSE
                           CTEMP = CONJG( CTEMP )
                           CT = CONJG( C )
                           ST = CONJG( S )
                        END IF
                        CALL CLAROT( .TRUE., .TRUE., .TRUE., K+2, C, S, A( JCH-ISKEW*ICOL+IOFFG, ICOL ), ILDA, EXTRA, CTEMP )
                        IL = MIN( N+1-JCH, K+2 )
                        EXTRA = CZERO
                        CALL CLAROT( .FALSE., .TRUE., N-JCH.GT.K, IL, CT, ST, A( ( 1-ISKEW )*JCH+IOFFG, JCH ), ILDA, CTEMP, EXTRA )
                        ICOL = JCH
  360                CONTINUE
  370             CONTINUE
  380          CONTINUE
*
*              If we need upper triangle, copy from lower. Note that
*              the order of copying is chosen to work for 'b' -> 'q'
*
               IF( IPACK.NE.IPACKG .AND. IPACK.NE.4 ) THEN
                  DO 410 JC = N, 1, -1
                     IROW = IOFFST - ISKEW*JC
                     IF( CSYM ) THEN
                        DO 390 JR = JC, MAX( 1, JC-UUB ), -1
                           A( JR+IROW, JC ) = A( JC-ISKEW*JR+IOFFG, JR )
  390                   CONTINUE
                     ELSE
                        DO 400 JR = JC, MAX( 1, JC-UUB ), -1
                           A( JR+IROW, JC ) = CONJG( A( JC-ISKEW*JR+ IOFFG, JR ) )
  400                   CONTINUE
                     END IF
  410             CONTINUE
                  IF( IPACK.EQ.6 ) THEN
                     DO 430 JC = 1, UUB
                        DO 420 JR = 1, UUB + 1 - JC
                           A( JR, JC ) = CZERO
  420                   CONTINUE
  430                CONTINUE
                  END IF
                  IF( IPACKG.EQ.5 ) THEN
                     IPACKG = IPACK
                  ELSE
                     IPACKG = 0
                  END IF
               END IF
            END IF
*
*           Ensure that the diagonal is real if Hermitian
*
            IF( .NOT.CSYM ) THEN
               DO 440 JC = 1, N
                  IROW = IOFFST + ( 1-ISKEW )*JC
                  A( IROW, JC ) = CMPLX( REAL( A( IROW, JC ) ) )
  440          CONTINUE
            END IF
*
         END IF
*
      ELSE
*
*        4)      Generate Banded Matrix by first
*                Rotating by random Unitary matrices,
*                then reducing the bandwidth using Householder
*                transformations.
*
*                Note: we should get here only if LDA .ge. N
*
         IF( ISYM.EQ.1 ) THEN
*
*           Non-symmetric -- A = U D V
*
            CALL CLAGGE( MR, NC, LLB, UUB, D, A, LDA, ISEED, WORK, IINFO )
         ELSE
*
*           Symmetric -- A = U D U' or
*           Hermitian -- A = U D U*
*
            IF( CSYM ) THEN
               CALL CLAGSY( M, LLB, D, A, LDA, ISEED, WORK, IINFO )
            ELSE
               CALL CLAGHE( M, LLB, D, A, LDA, ISEED, WORK, IINFO )
            END IF
         END IF
*
         IF( IINFO.NE.0 ) THEN
            INFO = 3
            RETURN
         END IF
      END IF
*
*     5)      Pack the matrix
*
      IF( IPACK.NE.IPACKG ) THEN
         IF( IPACK.EQ.1 ) THEN
*
*           'U' -- Upper triangular, not packed
*
            DO 460 J = 1, M
               DO 450 I = J + 1, M
                  A( I, J ) = CZERO
  450          CONTINUE
  460       CONTINUE
*
         ELSE IF( IPACK.EQ.2 ) THEN
*
*           'L' -- Lower triangular, not packed
*
            DO 480 J = 2, M
               DO 470 I = 1, J - 1
                  A( I, J ) = CZERO
  470          CONTINUE
  480       CONTINUE
*
         ELSE IF( IPACK.EQ.3 ) THEN
*
*           'C' -- Upper triangle packed Columnwise.
*
            ICOL = 1
            IROW = 0
            DO 500 J = 1, M
               DO 490 I = 1, J
                  IROW = IROW + 1
                  IF( IROW.GT.LDA ) THEN
                     IROW = 1
                     ICOL = ICOL + 1
                  END IF
                  A( IROW, ICOL ) = A( I, J )
  490          CONTINUE
  500       CONTINUE
*
         ELSE IF( IPACK.EQ.4 ) THEN
*
*           'R' -- Lower triangle packed Columnwise.
*
            ICOL = 1
            IROW = 0
            DO 520 J = 1, M
               DO 510 I = J, M
                  IROW = IROW + 1
                  IF( IROW.GT.LDA ) THEN
                     IROW = 1
                     ICOL = ICOL + 1
                  END IF
                  A( IROW, ICOL ) = A( I, J )
  510          CONTINUE
  520       CONTINUE
*
         ELSE IF( IPACK.GE.5 ) THEN
*
*           'B' -- The lower triangle is packed as a band matrix.
*           'Q' -- The upper triangle is packed as a band matrix.
*           'Z' -- The whole matrix is packed as a band matrix.
*
            IF( IPACK.EQ.5 ) UUB = 0             IF( IPACK.EQ.6 ) LLB = 0
*
            DO 540 J = 1, UUB
               DO 530 I = MIN( J+LLB, M ), 1, -1
                  A( I-J+UUB+1, J ) = A( I, J )
  530          CONTINUE
  540       CONTINUE
*
            DO 560 J = UUB + 2, N
               DO 550 I = J - UUB, MIN( J+LLB, M )
                  A( I-J+UUB+1, J ) = A( I, J )
  550          CONTINUE
  560       CONTINUE
         END IF
*
*        If packed, zero out extraneous elements.
*
*        Symmetric/Triangular Packed --
*        zero out everything after A(IROW,ICOL)
*
         IF( IPACK.EQ.3 .OR. IPACK.EQ.4 ) THEN
            DO 580 JC = ICOL, M
               DO 570 JR = IROW + 1, LDA
                  A( JR, JC ) = CZERO
  570          CONTINUE
               IROW = 0
  580       CONTINUE
*
         ELSE IF( IPACK.GE.5 ) THEN
*
*           Packed Band --
*              1st row is now in A( UUB+2-j, j), zero above it
*              m-th row is now in A( M+UUB-j,j), zero below it
*              last non-zero diagonal is now in A( UUB+LLB+1,j ),
*                 zero below it, too.
*
            IR1 = UUB + LLB + 2
            IR2 = UUB + M + 2
            DO 610 JC = 1, N
               DO 590 JR = 1, UUB + 1 - JC
                  A( JR, JC ) = CZERO
  590          CONTINUE
               DO 600 JR = MAX( 1, MIN( IR1, IR2-JC ) ), LDA
                  A( JR, JC ) = CZERO
  600          CONTINUE
  610       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of CLATMT
*
      END
