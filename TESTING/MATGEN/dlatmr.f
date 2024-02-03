      SUBROUTINE DLATMR( M, N, DIST, ISEED, SYM, D, MODE, COND, DMAX, RSIGN, GRADE, DL, MODEL, CONDL, DR, MODER, CONDR, PIVTNG, IPIVOT, KL, KU, SPARSE, ANORM, PACK, A, LDA, IWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             DIST, GRADE, PACK, PIVTNG, RSIGN, SYM;
      int                INFO, KL, KU, LDA, M, MODE, MODEL, MODER, N;
      double             ANORM, COND, CONDL, CONDR, DMAX, SPARSE;
*     ..
*     .. Array Arguments ..
      int                IPIVOT( * ), ISEED( 4 ), IWORK( * );
      double             A( LDA, * ), D( * ), DL( * ), DR( * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO;
      PARAMETER          ( ZERO = 0.0D0 )
      double             ONE;
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      bool               BADPVT, DZERO, FULBND;
      int                I, IDIST, IGRADE, IISUB, IPACK, IPVTNG, IRSIGN, ISUB, ISYM, J, JJSUB, JSUB, K, KLL, KUU, MNMIN, MNSUB, MXSUB, NPVTS;
      double             ALPHA, ONORM, TEMP;
*     ..
*     .. Local Arrays ..
      double             TEMPA( 1 );
*     ..
*     .. External Functions ..
      bool               LSAME;
      double             DLANGB, DLANGE, DLANSB, DLANSP, DLANSY, DLATM2, DLATM3       EXTERNAL           LSAME, DLANGB, DLANGE, DLANSB, DLANSP, DLANSY, DLATM2, DLATM3;
*     ..
*     .. External Subroutines ..
      // EXTERNAL DLATM1, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, MOD
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
      IF( LSAME( SYM, 'S' ) ) THEN
         ISYM = 0
      ELSE IF( LSAME( SYM, 'N' ) ) THEN
         ISYM = 1
      ELSE IF( LSAME( SYM, 'H' ) ) THEN
         ISYM = 0
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
      ELSE IF( LSAME( GRADE, 'H' ) .OR. LSAME( GRADE, 'S' ) ) THEN
         IGRADE = 5
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
            IF( DL( I ).EQ.ZERO ) DZERO = .TRUE.
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
      ELSE IF( M.NE.N .AND. ISYM.EQ.0 ) THEN
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
      ELSE IF( IGRADE.EQ.-1 .OR. ( IGRADE.EQ.4 .AND. M.NE.N ) .OR. ( ( IGRADE.GE.1 .AND. IGRADE.LE.4 ) .AND. ISYM.EQ.0 ) ) THEN
         INFO = -11
      ELSE IF( IGRADE.EQ.4 .AND. DZERO ) THEN
         INFO = -12
      ELSE IF( ( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ.5 ) .AND. ( MODEL.LT.-6 .OR. MODEL.GT.6 ) ) THEN
         INFO = -13
      ELSE IF( ( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ.5 ) .AND. ( MODEL.NE.-6 .AND. MODEL.NE.0 .AND. MODEL.NE.6 ) .AND. CONDL.LT.ONE ) THEN
         INFO = -14
      ELSE IF( ( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) .AND. ( MODER.LT.-6 .OR. MODER.GT.6 ) ) THEN
         INFO = -16
      ELSE IF( ( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) .AND. ( MODER.NE.-6 .AND. MODER.NE.0 .AND. MODER.NE.6 ) .AND. CONDR.LT.ONE ) THEN
         INFO = -17
      ELSE IF( IPVTNG.EQ.-1 .OR. ( IPVTNG.EQ.3 .AND. M.NE.N ) .OR. ( ( IPVTNG.EQ.1 .OR. IPVTNG.EQ.2 ) .AND. ISYM.EQ.0 ) ) THEN
         INFO = -18
      ELSE IF( IPVTNG.NE.0 .AND. BADPVT ) THEN
         INFO = -19
      ELSE IF( KL.LT.0 ) THEN
         INFO = -20
      ELSE IF( KU.LT.0 .OR. ( ISYM.EQ.0 .AND. KL.NE.KU ) ) THEN
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
         CALL XERBLA( 'DLATMR', -INFO )
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
      CALL DLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, MNMIN, INFO )
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
         IF( TEMP.EQ.ZERO .AND. DMAX.NE.ZERO ) THEN
            INFO = 2
            RETURN
         END IF
         IF( TEMP.NE.ZERO ) THEN
            ALPHA = DMAX / TEMP
         ELSE
            ALPHA = ONE
         END IF
         DO 50 I = 1, MNMIN
            D( I ) = ALPHA*D( I )
   50    CONTINUE
*
      END IF
*
*     Compute DL if grading set
*
      IF( IGRADE.EQ.1 .OR. IGRADE.EQ.3 .OR. IGRADE.EQ.4 .OR. IGRADE.EQ. 5 ) THEN
         CALL DLATM1( MODEL, CONDL, 0, IDIST, ISEED, DL, M, INFO )
         IF( INFO.NE.0 ) THEN
            INFO = 3
            RETURN
         END IF
      END IF
*
*     Compute DR if grading set
*
      IF( IGRADE.EQ.2 .OR. IGRADE.EQ.3 ) THEN
         CALL DLATM1( MODER, CONDR, 0, IDIST, ISEED, DR, N, INFO )
         IF( INFO.NE.0 ) THEN
            INFO = 4
            RETURN
         END IF
      END IF
*
*     3)     Generate IWORK if pivoting
*
      IF( IPVTNG.GT.0 ) THEN
         DO 60 I = 1, NPVTS
            IWORK( I ) = I
   60    CONTINUE
         IF( FULBND ) THEN
            DO 70 I = 1, NPVTS
               K = IPIVOT( I )
               J = IWORK( I )
               IWORK( I ) = IWORK( K )
               IWORK( K ) = J
   70       CONTINUE
         ELSE
            DO 80 I = NPVTS, 1, -1
               K = IPIVOT( I )
               J = IWORK( I )
               IWORK( I ) = IWORK( K )
               IWORK( K ) = J
   80       CONTINUE
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
*        Use DLATM3 so matrices generated with differing PIVOTing only
*        differ only in the order of their rows and/or columns.
*
         IF( IPACK.EQ.0 ) THEN
            IF( ISYM.EQ.0 ) THEN
               DO 100 J = 1, N
                  DO 90 I = 1, J
                     TEMP = DLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = TEMP
                     A( JSUB, ISUB ) = TEMP
   90             CONTINUE
  100          CONTINUE
            ELSE IF( ISYM.EQ.1 ) THEN
               DO 120 J = 1, N
                  DO 110 I = 1, M
                     TEMP = DLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB, JSUB ) = TEMP
  110             CONTINUE
  120          CONTINUE
            END IF
*
         ELSE IF( IPACK.EQ.1 ) THEN
*
            DO 140 J = 1, N
               DO 130 I = 1, J
                  TEMP = DLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  A( MNSUB, MXSUB ) = TEMP
                  IF( MNSUB.NE.MXSUB ) A( MXSUB, MNSUB ) = ZERO
  130          CONTINUE
  140       CONTINUE
*
         ELSE IF( IPACK.EQ.2 ) THEN
*
            DO 160 J = 1, N
               DO 150 I = 1, J
                  TEMP = DLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  A( MXSUB, MNSUB ) = TEMP
                  IF( MNSUB.NE.MXSUB ) A( MNSUB, MXSUB ) = ZERO
  150          CONTINUE
  160       CONTINUE
*
         ELSE IF( IPACK.EQ.3 ) THEN
*
            DO 180 J = 1, N
               DO 170 I = 1, J
                  TEMP = DLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
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
                  A( IISUB, JJSUB ) = TEMP
  170          CONTINUE
  180       CONTINUE
*
         ELSE IF( IPACK.EQ.4 ) THEN
*
            DO 200 J = 1, N
               DO 190 I = 1, J
                  TEMP = DLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
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
                  A( IISUB, JJSUB ) = TEMP
  190          CONTINUE
  200       CONTINUE
*
         ELSE IF( IPACK.EQ.5 ) THEN
*
            DO 220 J = 1, N
               DO 210 I = J - KUU, J
                  IF( I.LT.1 ) THEN
                     A( J-I+1, I+N ) = ZERO
                  ELSE
                     TEMP = DLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     MNSUB = MIN( ISUB, JSUB )
                     MXSUB = MAX( ISUB, JSUB )
                     A( MXSUB-MNSUB+1, MNSUB ) = TEMP
                  END IF
  210          CONTINUE
  220       CONTINUE
*
         ELSE IF( IPACK.EQ.6 ) THEN
*
            DO 240 J = 1, N
               DO 230 I = J - KUU, J
                  TEMP = DLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  MNSUB = MIN( ISUB, JSUB )
                  MXSUB = MAX( ISUB, JSUB )
                  A( MNSUB-MXSUB+KUU+1, MXSUB ) = TEMP
  230          CONTINUE
  240       CONTINUE
*
         ELSE IF( IPACK.EQ.7 ) THEN
*
            IF( ISYM.EQ.0 ) THEN
               DO 260 J = 1, N
                  DO 250 I = J - KUU, J
                     TEMP = DLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     MNSUB = MIN( ISUB, JSUB )
                     MXSUB = MAX( ISUB, JSUB )
                     A( MNSUB-MXSUB+KUU+1, MXSUB ) = TEMP
                     IF( I.LT.1 ) A( J-I+1+KUU, I+N ) = ZERO                      IF( I.GE.1 .AND. MNSUB.NE.MXSUB ) A( MXSUB-MNSUB+1+KUU, MNSUB ) = TEMP
  250             CONTINUE
  260          CONTINUE
            ELSE IF( ISYM.EQ.1 ) THEN
               DO 280 J = 1, N
                  DO 270 I = J - KUU, J + KLL
                     TEMP = DLATM3( M, N, I, J, ISUB, JSUB, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( ISUB-JSUB+KUU+1, JSUB ) = TEMP
  270             CONTINUE
  280          CONTINUE
            END IF
*
         END IF
*
      ELSE
*
*        Use DLATM2
*
         IF( IPACK.EQ.0 ) THEN
            IF( ISYM.EQ.0 ) THEN
               DO 300 J = 1, N
                  DO 290 I = 1, J
                     A( I, J ) = DLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     A( J, I ) = A( I, J )
  290             CONTINUE
  300          CONTINUE
            ELSE IF( ISYM.EQ.1 ) THEN
               DO 320 J = 1, N
                  DO 310 I = 1, M
                     A( I, J ) = DLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  310             CONTINUE
  320          CONTINUE
            END IF
*
         ELSE IF( IPACK.EQ.1 ) THEN
*
            DO 340 J = 1, N
               DO 330 I = 1, J
                  A( I, J ) = DLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )                   IF( I.NE.J ) A( J, I ) = ZERO
  330          CONTINUE
  340       CONTINUE
*
         ELSE IF( IPACK.EQ.2 ) THEN
*
            DO 360 J = 1, N
               DO 350 I = 1, J
                  A( J, I ) = DLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )                   IF( I.NE.J ) A( I, J ) = ZERO
  350          CONTINUE
  360       CONTINUE
*
         ELSE IF( IPACK.EQ.3 ) THEN
*
            ISUB = 0
            JSUB = 1
            DO 380 J = 1, N
               DO 370 I = 1, J
                  ISUB = ISUB + 1
                  IF( ISUB.GT.LDA ) THEN
                     ISUB = 1
                     JSUB = JSUB + 1
                  END IF
                  A( ISUB, JSUB ) = DLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  370          CONTINUE
  380       CONTINUE
*
         ELSE IF( IPACK.EQ.4 ) THEN
*
            IF( ISYM.EQ.0 ) THEN
               DO 400 J = 1, N
                  DO 390 I = 1, J
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
                     A( ISUB, JSUB ) = DLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  390             CONTINUE
  400          CONTINUE
            ELSE
               ISUB = 0
               JSUB = 1
               DO 420 J = 1, N
                  DO 410 I = J, M
                     ISUB = ISUB + 1
                     IF( ISUB.GT.LDA ) THEN
                        ISUB = 1
                        JSUB = JSUB + 1
                     END IF
                     A( ISUB, JSUB ) = DLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  410             CONTINUE
  420          CONTINUE
            END IF
*
         ELSE IF( IPACK.EQ.5 ) THEN
*
            DO 440 J = 1, N
               DO 430 I = J - KUU, J
                  IF( I.LT.1 ) THEN
                     A( J-I+1, I+N ) = ZERO
                  ELSE
                     A( J-I+1, I ) = DLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                  END IF
  430          CONTINUE
  440       CONTINUE
*
         ELSE IF( IPACK.EQ.6 ) THEN
*
            DO 460 J = 1, N
               DO 450 I = J - KUU, J
                  A( I-J+KUU+1, J ) = DLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  450          CONTINUE
  460       CONTINUE
*
         ELSE IF( IPACK.EQ.7 ) THEN
*
            IF( ISYM.EQ.0 ) THEN
               DO 480 J = 1, N
                  DO 470 I = J - KUU, J
                     A( I-J+KUU+1, J ) = DLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
                     IF( I.LT.1 ) A( J-I+1+KUU, I+N ) = ZERO                      IF( I.GE.1 .AND. I.NE.J ) A( J-I+1+KUU, I ) = A( I-J+KUU+1, J )
  470             CONTINUE
  480          CONTINUE
            ELSE IF( ISYM.EQ.1 ) THEN
               DO 500 J = 1, N
                  DO 490 I = J - KUU, J + KLL
                     A( I-J+KUU+1, J ) = DLATM2( M, N, I, J, KL, KU, IDIST, ISEED, D, IGRADE, DL, DR, IPVTNG, IWORK, SPARSE )
  490             CONTINUE
  500          CONTINUE
            END IF
*
         END IF
*
      END IF
*
*     5)      Scaling the norm
*
      IF( IPACK.EQ.0 ) THEN
         ONORM = DLANGE( 'M', M, N, A, LDA, TEMPA )
      ELSE IF( IPACK.EQ.1 ) THEN
         ONORM = DLANSY( 'M', 'U', N, A, LDA, TEMPA )
      ELSE IF( IPACK.EQ.2 ) THEN
         ONORM = DLANSY( 'M', 'L', N, A, LDA, TEMPA )
      ELSE IF( IPACK.EQ.3 ) THEN
         ONORM = DLANSP( 'M', 'U', N, A, TEMPA )
      ELSE IF( IPACK.EQ.4 ) THEN
         ONORM = DLANSP( 'M', 'L', N, A, TEMPA )
      ELSE IF( IPACK.EQ.5 ) THEN
         ONORM = DLANSB( 'M', 'L', N, KLL, A, LDA, TEMPA )
      ELSE IF( IPACK.EQ.6 ) THEN
         ONORM = DLANSB( 'M', 'U', N, KUU, A, LDA, TEMPA )
      ELSE IF( IPACK.EQ.7 ) THEN
         ONORM = DLANGB( 'M', N, KLL, KUU, A, LDA, TEMPA )
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
               DO 510 J = 1, N
                  CALL DSCAL( M, ONE / ONORM, A( 1, J ), 1 )
                  CALL DSCAL( M, ANORM, A( 1, J ), 1 )
  510          CONTINUE
*
            ELSE IF( IPACK.EQ.3 .OR. IPACK.EQ.4 ) THEN
*
               CALL DSCAL( N*( N+1 ) / 2, ONE / ONORM, A, 1 )
               CALL DSCAL( N*( N+1 ) / 2, ANORM, A, 1 )
*
            ELSE IF( IPACK.GE.5 ) THEN
*
               DO 520 J = 1, N
                  CALL DSCAL( KLL+KUU+1, ONE / ONORM, A( 1, J ), 1 )
                  CALL DSCAL( KLL+KUU+1, ANORM, A( 1, J ), 1 )
  520          CONTINUE
*
            END IF
*
         ELSE
*
*           Scale straightforwardly
*
            IF( IPACK.LE.2 ) THEN
               DO 530 J = 1, N
                  CALL DSCAL( M, ANORM / ONORM, A( 1, J ), 1 )
  530          CONTINUE
*
            ELSE IF( IPACK.EQ.3 .OR. IPACK.EQ.4 ) THEN
*
               CALL DSCAL( N*( N+1 ) / 2, ANORM / ONORM, A, 1 )
*
            ELSE IF( IPACK.GE.5 ) THEN
*
               DO 540 J = 1, N
                  CALL DSCAL( KLL+KUU+1, ANORM / ONORM, A( 1, J ), 1 )
  540          CONTINUE
            END IF
*
         END IF
*
      END IF
*
*     End of DLATMR
*
      END
