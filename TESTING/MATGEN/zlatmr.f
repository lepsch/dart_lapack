      SUBROUTINE ZLATMR( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RSIGN, GRADE, DL, MODEL, CONDL, DR, MODER, CONDR, PIVTNG, IPIVOT, KL, KU, SPARSE, ANORM, PACK, A, LDA, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             DIST, GRADE, PACK, PIVTNG, RSIGN, SYM;
      int                INFO, KL, KU, LDA, M, MODE, MODEL, MODER, N
      DOUBLE PRECISION   ANORM, COND, CONDL, CONDR, SPARSE
      COMPLEX*16         DMAX
*     ..
*     .. Array Arguments ..
      int                IPIVOT( * ), ISEED( 4 ), IWORK( * )
      COMPLEX*16         A( LDA, * ), D( * ), DL( * ), DR( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            BADPVT, DZERO, FULBND
      int                I, IDIST, IGRADE, IISUB, IPACK, IPVTNG, IRSIGN, ISUB, ISYM, J, JJSUB, JSUB, K, KLL, KUU, MNMIN, MNSUB, MXSUB, NPVTS
      DOUBLE PRECISION   ONORM, TEMP
      COMPLEX*16         CALPHA, CTEMP
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   TEMPA( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   ZLANGB, ZLANGE, ZLANSB, ZLANSP, ZLANSY
      COMPLEX*16         ZLATM2, ZLATM3
      EXTERNAL           LSAME, ZLANGB, ZLANGE, ZLANSB, ZLANSP, ZLANSY, ZLATM2, ZLATM3
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDSCAL, ZLATM1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, MAX, MIN, MOD
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
      ELSE IF( LSAME( DIST, 'D' ) ) THEN
         IDIST = 4
      ELSE
         IDIST = -1
      END IF
*
*     Decode SYM
*
      IF( LSAME( SYM, 'H' ) ) THEN
         ISYM = 0
      ELSE IF( LSAME( SYM, 'N' ) ) THEN
         ISYM = 1
      ELSE IF( LSAME( SYM, 'S' ) ) THEN
         ISYM = 2
      ELSE
         ISYM = -1
      END IF
*
*     Decode RSIGN
*
      IF( LSAME( RSIGN, 'F' ) ) THEN
         IRSIGN = 0
      ELSE IF( LSAME( RSIGN, 'T' ) ) THEN
         IRSIGN = 1
      ELSE
         IRSIGN = -1
      END IF
*
*     Decode PIVTNG
*
      IF( LSAME( PIVTNG, 'N' ) ) THEN
         IPVTNG = 0
      ELSE IF( LSAME( PIVTNG, ' ' ) ) THEN
         IPVTNG = 0
      ELSE IF( LSAME( PIVTNG, 'L' ) ) THEN
         IPVTNG = 1
         NPVTS = M
      ELSE IF( LSAME( PIVTNG, 'R' ) ) THEN
         IPVTNG = 2
         NPVTS = N
      ELSE IF( LSAME( PIVTNG, 'B' ) ) THEN
         IPVTNG = 3
         NPVTS = MIN( N, M )
      ELSE IF( LSAME( PIVTNG, 'F' ) ) THEN
         IPVTNG = 3
         NPVTS = MIN( N, M )
      ELSE
         IPVTNG = -1
      END IF
*
*     Decode GRADE
*
      IF( LSAME( GRADE, 'N' ) ) THEN
         IGRADE = 0
      ELSE IF( LSAME( GRADE, 'L' ) ) THEN
         IGRADE = 1
      ELSE IF( LSAME( GRADE, 'R' ) ) THEN
         IGRADE = 2
      ELSE IF( LSAME( GRADE, 'B' ) ) THEN
         IGRADE = 3
      ELSE IF( LSAME( GRADE, 'E' ) ) THEN
         IGRADE = 4
      ELSE IF( LSAME( GRADE, 'H' ) ) THEN
         IGRADE = 5
      ELSE IF( LSAME( GRADE, 'S' ) ) THEN
         IGRADE = 6
      ELSE
         IGRADE = -1
      END IF
*
*     Decode PACK
*
      IF( LSAME( PACK, 'N' ) ) THEN
         IPACK = 0
      ELSE IF( LSAME( PACK, 'U' ) ) THEN
         IPACK = 1
      ELSE IF( LSAME( PACK, 'L' ) ) THEN
         IPACK = 2
      ELSE IF( LSAME( PACK, 'C' ) ) THEN
         IPACK = 3
      ELSE IF( LSAME( PACK, 'R' ) ) THEN
         IPACK = 4
      ELSE IF( LSAME( PACK, 'B' ) ) THEN
         IPACK = 5
      ELSE IF( LSAME( PACK, 'Q' ) ) THEN
         IPACK = 6
      ELSE IF( LSAME( PACK, 'Z' ) ) THEN
         IPACK = 7
      ELSE
         IPACK = -1
      END IF
*
*     Set certain internal parameters
*
      MNMIN = MIN( M, N )
      KLL = MIN( KL, M-1 )
      KUU = MIN( KU, N-1 )
*
*     If inv(DL) is used, check to see if DL has a zero entry.
*
      DZERO = .FALSE.
      IF( IGRADE.EQ.4 .AND. MODEL.EQ.0 ) THEN
         DO 10 I = 1, M
            IF( DL( I ).EQ.CZERO ) DZERO = .TRUE.
   10    CONTINUE
      END IF
*
*     Check values in IPIVOT
*
      BADPVT = .FALSE.
      IF( IPVTNG.GT.0 ) THEN
         DO 20 J = 1, NPVTS
            IF( IPIVOT( J ).LE.0 .OR. IPIVOT( J ).GT.NPVTS ) BADPVT = .TRUE.
   20    CONTINUE
      END IF
*
*     Set INFO if an error
*
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.NE.N .AND. ( ISYM.EQ.0 .OR. ISYM.EQ.2 ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( IDIST.EQ.-1 ) THEN
         INFO = -3
      ELSE IF( ISYM.EQ.-1 ) THEN
         INFO = -5
      ELSE IF( MODE.LT.-6 .OR. MODE.GT.6 ) THEN
         INFO = -7
      ELSE IF( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. COND.LT.ONE ) THEN
         INFO = -8
      ELSE IF( ( MODE.NE.-6 .AND. MODE.NE.0 .AND. MODE.NE.6 ) .AND. IRSIGN.EQ.-1 ) THEN
         INFO = -10
      ELSE IF( IGRADE.EQ.-1 .OR. ( IGRADE.EQ.4 .AND. M.NE.N ) .OR. ( ( IGRADE.EQ.1 .OR. IGRADE.EQ.2 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ.6 ) .AND. ISYM.EQ.0 ) .OR. ( ( IGRADE.EQ.1 .OR. IGRADE.EQ.2 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ.5 ) .AND. ISYM.EQ.2 ) ) THEN
         INFO = -11
      ELSE IF( IGRADE.EQ.4 .AND. DZERO ) THEN
         INFO = -12
      ELSE IF( ( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ.5 .OR. IGRADE.EQ.6 ) .AND. ( MODEL.LT.-6 .OR. MODEL.GT.6 ) ) THEN
         INFO = -13
      ELSE IF( ( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ.5 .OR. IGRADE.EQ.6 ) .AND. ( MODEL.NE.-6 .AND. MODEL.NE.0 .AND. MODEL.NE.6 ) .AND. CONDL.LT.ONE ) THEN
         INFO = -14
      ELSE IF( ( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) .AND. ( MODER.LT.-6 .OR. MODER.GT.6 ) ) THEN
         INFO = -16
      ELSE IF( ( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) .AND. ( MODER.NE.-6 .AND. MODER.NE.0 .AND. MODER.NE.6 ) .AND. CONDR.LT.ONE ) THEN
         INFO = -17
      ELSE IF( IPVTNG.EQ.-1 .OR. ( IPVTNG.EQ.3 .AND. M.NE.N ) .OR. ( ( IPVTNG.EQ.1 .OR. IPVTNG.EQ.2 ) .AND. ( ISYM.EQ.0 .OR. ISYM.EQ.2 ) ) ) THEN
         INFO = -18
      ELSE IF( IPVTNG.NE.0 .AND. BADPVT ) THEN
         INFO = -19
      ELSE IF( KL.LT.0 ) THEN
         INFO = -20
      ELSE IF( KU.LT.0 .OR. ( ( ISYM.EQ.0 .OR. ISYM.EQ.2 ) .AND. KL.NE. KU ) ) THEN
         INFO = -21
      ELSE IF( SPARSE.LT.ZERO .OR. SPARSE.GT.ONE ) THEN
         INFO = -22
      ELSE IF( IPACK.EQ.-1 .OR. ( ( IPACK.EQ.1 .OR. IPACK.EQ.2 .OR. IPACK.EQ.5 .OR. IPACK.EQ.6 ) .AND. ISYM.EQ.1 ) .OR. ( IPACK.EQ.3 .AND. ISYM.EQ.1 .AND. ( KL.NE.0 .OR. M.NE. N ) ) .OR. ( IPACK.EQ.4 .AND. ISYM.EQ.1 .AND. ( KU.NE. 0 .OR. M.NE.N ) ) ) THEN
         INFO = -24
      ELSE IF( ( ( IPACK.EQ.0 .OR. IPACK.EQ.1 .OR. IPACK.EQ.2 ) .AND. LDA.LT.MAX( 1, M ) ) .OR. ( ( IPACK.EQ.3 .OR. IPACK.EQ. 4 ) .AND. LDA.LT.1 ) .OR. ( ( IPACK.EQ.5 .OR. IPACK.EQ. 6 ) .AND. LDA.LT.KUU+1 ) .OR. ( IPACK.EQ.7 .AND. LDA.LT.KLL+KUU+1 ) ) THEN
         INFO = -26
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLATMR', -INFO )
         RETURN
      END IF
*
*     Decide if we can pivot consistently
*
      FULBND = .FALSE.
      IF( KUU.EQ.N-1 .AND. KLL.EQ.M-1 ) FULBND = .TRUE.
*
*     Initialize random number generator
*
      DO 30 I = 1, 4
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
   30 CONTINUE
*
      ISEED( 4 ) = 2*( ISEED( 4 ) / 2 ) + 1
*
*     2)      Set up D, DL, and DR, if indicated.
*
*             Compute D according to COND and MODE
*
      CALL ZLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
      IF( MODE.NE.0 .AND. MODE.NE.-6 .AND. MODE.NE.6 ) THEN
*
*        Scale by DMAX
*
         TEMP = ABS( D( 1 ) )
         DO 40 I = 2, MNMIN
            TEMP = MAX( TEMP, ABS( D( I ) ) )
   40    CONTINUE
         IF( TEMP.EQ.ZERO .AND. DMAX.NE.CZERO ) THEN
            INFO = 2
            RETURN
         END IF
         IF( TEMP.NE.ZERO ) THEN
            CALPHA = DMAX / TEMP
         ELSE
            CALPHA = CONE
         END IF
         DO 50 I = 1, MNMIN
            D( I ) = CALPHA*D( I )
   50    CONTINUE
*
      END IF
*
*     If matrix Hermitian, make D real
*
      IF( ISYM.EQ.0 ) THEN
         DO 60 I = 1, MNMIN
            D( I ) = DBLE( D( I ) )
   60    CONTINUE
      END IF
*
*     Compute DL if grading set
*
      IF( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ. 5 .OR. IGRADE.EQ.6 ) THEN
         CALL ZLATM1( MODEL, CONDL, 0, IDIST, ISEED, DL, M, INFO )
         IF( INFO.NE.0 ) THEN
            INFO = 3
            RETURN
         END IF
      END IF
*
*     Compute DR if grading set
*
      IF( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) THEN
         CALL ZLATM1( MODER, CONDR, 0, IDIST, ISEED, DR, N, INFO )
         IF( INFO.NE.0 ) THEN
            INFO = 4
            RETURN
         END IF
      END IF
*
*     3)     Generate IWORK if pivoting
*
      IF( IPVTNG.GT.0 ) THEN
         DO 70 I = 1, NPVTS
            IWORK( I ) = I
   70    CONTINUE
         IF( FULBND ) THEN
            DO 80 I = 1, NPVTS
               K = IPIVOT( I )
               J = IWORK( I )
               IWORK( I ) = IWORK( K )
               IWORK( K ) = J
   80       CONTINUE
         ELSE
            DO 90 I = NPVTS, 1, -1
               K = IPIVOT( I )
               J = IWORK( I )
               IWORK( I ) = IWORK( K )
               IWORK( K ) = J
   90       CONTINUE
         END IF
      END IF
*
*     4)      Generate matrices for each kind of PACKing
*             Always sweep matrix columnwise (if symmetric, upper
*             half only) so that matrix generated does not depend
*             on PACK
*
      IF( FULBND ) THEN
*
*        Use ZLATM3 so matrices generated with differing PIVOTing only
*        differ only in the order of their rows and/or columns.
*
         IF( IPACK.EQ.0 ) THEN
            IF( ISYM.EQ.0 ) THEN
               DO 110 J = 1, N
                  DO 100 I = 1, J
                     CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = CTEMP
                     A( JSUB, ISUB ) = DCONJG( CTEMP )
  100             CONTINUE
  110          CONTINUE
            ELSE IF( ISYM.EQ.1 ) THEN
               DO 130 J = 1, N
                  DO 120 I = 1, M
                     CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = CTEMP
  120             CONTINUE
  130          CONTINUE
            ELSE IF( ISYM.EQ.2 ) THEN
               DO 150 J = 1, N
                  DO 140 I = 1, J
                     CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = CTEMP
                     A( JSUB, ISUB ) = CTEMP
  140             CONTINUE
  150          CONTINUE
            END IF
*
         ELSE IF( IPACK.EQ.1 ) THEN
*
            DO 170 J = 1, N
               DO 160 I = 1, J
                  CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  IF( MXSUB.EQ.ISUB .AND. ISYM.EQ.0 ) THEN
                     A( MNSUB, MXSUB ) = DCONJG( CTEMP )
                  ELSE
                     A( MNSUB, MXSUB ) = CTEMP
                  END IF
                  IF( MNSUB.NE.MXSUB ) A( MXSUB, MNSUB ) = CZERO
  160          CONTINUE
  170       CONTINUE
*
         ELSE IF( IPACK.EQ.2 ) THEN
*
            DO 190 J = 1, N
               DO 180 I = 1, J
                  CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  IF( MXSUB.EQ.JSUB .AND. ISYM.EQ.0 ) THEN
                     A( MXSUB, MNSUB ) = DCONJG( CTEMP )
                  ELSE
                     A( MXSUB, MNSUB ) = CTEMP
                  END IF
                  IF( MNSUB.NE.MXSUB ) A( MNSUB, MXSUB ) = CZERO
  180          CONTINUE
  190       CONTINUE
*
         ELSE IF( IPACK.EQ.3 ) THEN
*
            DO 210 J = 1, N
               DO 200 I = 1, J
                  CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
*
*                 Compute K = location of (ISUB,JSUB) entry in packed
*                 array
*
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  K = MXSUB*( MXSUB-1 ) / 2 + MNSUB
*
*                 Convert K to (IISUB,JJSUB) location
*
                  JJSUB = ( K-1 ) / LDA + 1
                  IISUB = K - LDA*( JJSUB-1 )
*
                  IF( MXSUB.EQ.ISUB .AND. ISYM.EQ.0 ) THEN
                     A( IISUB, JJSUB ) = DCONJG( CTEMP )
                  ELSE
                     A( IISUB, JJSUB ) = CTEMP
                  END IF
  200          CONTINUE
  210       CONTINUE
*
         ELSE IF( IPACK.EQ.4 ) THEN
*
            DO 230 J = 1, N
               DO 220 I = 1, J
                  CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
*
*                 Compute K = location of (I,J) entry in packed array
*
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  IF( MNSUB.EQ.1 ) THEN
                     K = MXSUB
                  ELSE
                     K = N*( N+1 ) / 2 - ( N-MNSUB+1 )*( N-MNSUB+2 ) / 2 + MXSUB - MNSUB + 1
                  END IF
*
*                 Convert K to (IISUB,JJSUB) location
*
                  JJSUB = ( K-1 ) / LDA + 1
                  IISUB = K - LDA*( JJSUB-1 )
*
                  IF( MXSUB.EQ.JSUB .AND. ISYM.EQ.0 ) THEN
                     A( IISUB, JJSUB ) = DCONJG( CTEMP )
                  ELSE
                     A( IISUB, JJSUB ) = CTEMP
                  END IF
  220          CONTINUE
  230       CONTINUE
*
         ELSE IF( IPACK.EQ.5 ) THEN
*
            DO 250 J = 1, N
               DO 240 I = J - KUU, J
                  IF( I.LT.1 ) THEN
                     A( J-I+1, I+N ) = CZERO
                  ELSE
                     CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     MNSUB = MIN( ISUB, JSUB )
                     MXSUB = MAX( ISUB, JSUB )
                     IF( MXSUB.EQ.JSUB .AND. ISYM.EQ.0 ) THEN
                        A( MXSUB-MNSUB+1, MNSUB ) = DCONJG( CTEMP )
                     ELSE
                        A( MXSUB-MNSUB+1, MNSUB ) = CTEMP
                     END IF
                  END IF
  240          CONTINUE
  250       CONTINUE
*
         ELSE IF( IPACK.EQ.6 ) THEN
*
            DO 270 J = 1, N
               DO 260 I = J - KUU, J
                  CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  IF( MXSUB.EQ.ISUB .AND. ISYM.EQ.0 ) THEN
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = DCONJG( CTEMP )
                  ELSE
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = CTEMP
                  END IF
  260          CONTINUE
  270       CONTINUE
*
         ELSE IF( IPACK.EQ.7 ) THEN
*
            IF( ISYM.NE.1 ) THEN
               DO 290 J = 1, N
                  DO 280 I = J - KUU, J
                     CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     MNSUB = MIN( ISUB, JSUB )
                     MXSUB = MAX( ISUB, JSUB )
                     IF( I.LT.1 ) A( J-I+1+KUU, I+N ) = CZERO
                     IF( MXSUB.EQ.ISUB .AND. ISYM.EQ.0 ) THEN
                        A( MNSUB-MXSUB+KUU+1, MXSUB ) = DCONJG( CTEMP )
                     ELSE
                        A( MNSUB-MXSUB+KUU+1, MXSUB ) = CTEMP
                     END IF
                     IF( I.GE.1 .AND. MNSUB.NE.MXSUB ) THEN
                        IF( MNSUB.EQ.ISUB .AND. ISYM.EQ.0 ) THEN
                           A( MXSUB-MNSUB+1+KUU, MNSUB ) = DCONJG( CTEMP )
                        ELSE
                           A( MXSUB-MNSUB+1+KUU, MNSUB ) = CTEMP
                        END IF
                     END IF
  280             CONTINUE
  290          CONTINUE
            ELSE IF( ISYM.EQ.1 ) THEN
               DO 310 J = 1, N
                  DO 300 I = J - KUU, J + KLL
                     CTEMP = ZLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB-JSUB+KUU+1, JSUB ) = CTEMP
  300             CONTINUE
  310          CONTINUE
            END IF
*
         END IF
*
      ELSE
*
*        Use ZLATM2
*
         IF( IPACK.EQ.0 ) THEN
            IF( ISYM.EQ.0 ) THEN
               DO 330 J = 1, N
                  DO 320 I = 1, J
                     A( I, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( J, I ) = DCONJG( A( I, J ) )
  320             CONTINUE
  330          CONTINUE
            ELSE IF( ISYM.EQ.1 ) THEN
               DO 350 J = 1, N
                  DO 340 I = 1, M
                     A( I, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  340             CONTINUE
  350          CONTINUE
            ELSE IF( ISYM.EQ.2 ) THEN
               DO 370 J = 1, N
                  DO 360 I = 1, J
                     A( I, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( J, I ) = A( I, J )
  360             CONTINUE
  370          CONTINUE
            END IF
*
         ELSE IF( IPACK.EQ.1 ) THEN
*
            DO 390 J = 1, N
               DO 380 I = 1, J
                  A( I, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )                   IF( I.NE.J ) A( J, I ) = CZERO
  380          CONTINUE
  390       CONTINUE
*
         ELSE IF( IPACK.EQ.2 ) THEN
*
            DO 410 J = 1, N
               DO 400 I = 1, J
                  IF( ISYM.EQ.0 ) THEN
                     A( J, I ) = DCONJG( ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE ) )
                  ELSE
                     A( J, I ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  END IF
                  IF( I.NE.J ) A( I, J ) = CZERO
  400          CONTINUE
  410       CONTINUE
*
         ELSE IF( IPACK.EQ.3 ) THEN
*
            ISUB = 0
            JSUB = 1
            DO 430 J = 1, N
               DO 420 I = 1, J
                  ISUB = ISUB + 1
                  IF( ISUB.GT.LDA ) THEN
                     ISUB = 1
                     JSUB = JSUB + 1
                  END IF
                  A( ISUB, JSUB ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  420          CONTINUE
  430       CONTINUE
*
         ELSE IF( IPACK.EQ.4 ) THEN
*
            IF( ISYM.EQ.0 .OR. ISYM.EQ.2 ) THEN
               DO 450 J = 1, N
                  DO 440 I = 1, J
*
*                    Compute K = location of (I,J) entry in packed array
*
                     IF( I.EQ.1 ) THEN
                        K = J
                     ELSE
                        K = N*( N+1 ) / 2 - ( N-I+1 )*( N-I+2 ) / 2 + J - I + 1
                     END IF
*
*                    Convert K to (ISUB,JSUB) location
*
                     JSUB = ( K-1 ) / LDA + 1
                     ISUB = K - LDA*( JSUB-1 )
*
                     A( ISUB, JSUB ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     IF( ISYM.EQ.0 ) A( ISUB, JSUB ) = DCONJG( A( ISUB, JSUB ) )
  440             CONTINUE
  450          CONTINUE
            ELSE
               ISUB = 0
               JSUB = 1
               DO 470 J = 1, N
                  DO 460 I = J, M
                     ISUB = ISUB + 1
                     IF( ISUB.GT.LDA ) THEN
                        ISUB = 1
                        JSUB = JSUB + 1
                     END IF
                     A( ISUB, JSUB ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  460             CONTINUE
  470          CONTINUE
            END IF
*
         ELSE IF( IPACK.EQ.5 ) THEN
*
            DO 490 J = 1, N
               DO 480 I = J - KUU, J
                  IF( I.LT.1 ) THEN
                     A( J-I+1, I+N ) = CZERO
                  ELSE
                     IF( ISYM.EQ.0 ) THEN
                        A( J-I+1, I ) = DCONJG( ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE ) )
                     ELSE
                        A( J-I+1, I ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     END IF
                  END IF
  480          CONTINUE
  490       CONTINUE
*
         ELSE IF( IPACK.EQ.6 ) THEN
*
            DO 510 J = 1, N
               DO 500 I = J - KUU, J
                  A( I-J+KUU+1, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  500          CONTINUE
  510       CONTINUE
*
         ELSE IF( IPACK.EQ.7 ) THEN
*
            IF( ISYM.NE.1 ) THEN
               DO 530 J = 1, N
                  DO 520 I = J - KUU, J
                     A( I-J+KUU+1, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     IF( I.LT.1 ) A( J-I+1+KUU, I+N ) = CZERO
                     IF( I.GE.1 .AND. I.NE.J ) THEN
                        IF( ISYM.EQ.0 ) THEN
                           A( J-I+1+KUU, I ) = DCONJG( A( I-J+KUU+1, J ) )
                        ELSE
                           A( J-I+1+KUU, I ) = A( I-J+KUU+1, J )
                        END IF
                     END IF
  520             CONTINUE
  530          CONTINUE
            ELSE IF( ISYM.EQ.1 ) THEN
               DO 550 J = 1, N
                  DO 540 I = J - KUU, J + KLL
                     A( I-J+KUU+1, J ) = ZLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  540             CONTINUE
  550          CONTINUE
            END IF
*
         END IF
*
      END IF
*
*     5)      Scaling the norm
*
      IF( IPACK.EQ.0 ) THEN
         ONORM = ZLANGE( 'M', M, N, A, LDA, TEMPA )
      ELSE IF( IPACK.EQ.1 ) THEN
         ONORM = ZLANSY( 'M', 'U', N, A, LDA, TEMPA )
      ELSE IF( IPACK.EQ.2 ) THEN
         ONORM = ZLANSY( 'M', 'L', N, A, LDA, TEMPA )
      ELSE IF( IPACK.EQ.3 ) THEN
         ONORM = ZLANSP( 'M', 'U', N, A, TEMPA )
      ELSE IF( IPACK.EQ.4 ) THEN
         ONORM = ZLANSP( 'M', 'L', N, A, TEMPA )
      ELSE IF( IPACK.EQ.5 ) THEN
         ONORM = ZLANSB( 'M', 'L', N, KLL, A, LDA, TEMPA )
      ELSE IF( IPACK.EQ.6 ) THEN
         ONORM = ZLANSB( 'M', 'U', N, KUU, A, LDA, TEMPA )
      ELSE IF( IPACK.EQ.7 ) THEN
         ONORM = ZLANGB( 'M', N, KLL, KUU, A, LDA, TEMPA )
      END IF
*
      IF( ANORM.GE.ZERO ) THEN
*
         IF( ANORM.GT.ZERO .AND. ONORM.EQ.ZERO ) THEN
*
*           Desired scaling impossible
*
            INFO = 5
            RETURN
*
         ELSE IF( ( ANORM.GT.ONE .AND. ONORM.LT.ONE ) .OR. ( ANORM.LT.ONE .AND. ONORM.GT.ONE ) ) THEN
*
*           Scale carefully to avoid over / underflow
*
            IF( IPACK.LE.2 ) THEN
               DO 560 J = 1, N
                  CALL ZDSCAL( M, ONE / ONORM, A( 1, J ), 1 )
                  CALL ZDSCAL( M, ANORM, A( 1, J ), 1 )
  560          CONTINUE
*
            ELSE IF( IPACK.EQ.3 .OR. IPACK.EQ.4 ) THEN
*
               CALL ZDSCAL( N*( N+1 ) / 2, ONE / ONORM, A, 1 )
               CALL ZDSCAL( N*( N+1 ) / 2, ANORM, A, 1 )
*
            ELSE IF( IPACK.GE.5 ) THEN
*
               DO 570 J = 1, N
                  CALL ZDSCAL( KLL+KUU+1, ONE / ONORM, A( 1, J ), 1 )
                  CALL ZDSCAL( KLL+KUU+1, ANORM, A( 1, J ), 1 )
  570          CONTINUE
*
            END IF
*
         ELSE
*
*           Scale straightforwardly
*
            IF( IPACK.LE.2 ) THEN
               DO 580 J = 1, N
                  CALL ZDSCAL( M, ANORM / ONORM, A( 1, J ), 1 )
  580          CONTINUE
*
            ELSE IF( IPACK.EQ.3 .OR. IPACK.EQ.4 ) THEN
*
               CALL ZDSCAL( N*( N+1 ) / 2, ANORM / ONORM, A, 1 )
*
            ELSE IF( IPACK.GE.5 ) THEN
*
               DO 590 J = 1, N
                  CALL ZDSCAL( KLL+KUU+1, ANORM / ONORM, A( 1, J ), 1 )
  590          CONTINUE
            END IF
*
         END IF
*
      END IF
*
*     End of ZLATMR
*
      END
