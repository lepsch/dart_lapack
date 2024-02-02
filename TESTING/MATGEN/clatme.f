      SUBROUTINE CLATME( N, DIST, ISEED, D, MODE, COND, DMAX,
     $  RSIGN,
     $                   UPPER, SIM, DS, MODES, CONDS, KL, KU, ANORM,
     $  A,
     $                   LDA, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIST, RSIGN, SIM, UPPER
      INTEGER            INFO, KL, KU, LDA, MODE, MODES, N
      REAL               ANORM, COND, CONDS
      COMPLEX            DMAX
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      REAL               DS( * )
      COMPLEX            A( LDA, * ), D( * ), WORK( * )
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
      COMPLEX            CONE
      PARAMETER          ( CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            BADS
      INTEGER            I, IC, ICOLS, IDIST, IINFO, IR, IROWS, IRSIGN,
     $                   ISIM, IUPPER, J, JC, JCR
      REAL               RALPHA, TEMP
      COMPLEX            ALPHA, TAU, XNORMS
*     ..
*     .. Local Arrays ..
      REAL               TEMPA( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               CLANGE
      COMPLEX            CLARND
      EXTERNAL           LSAME, CLANGE, CLARND
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CGEMV, CGERC, CLACGV, CLARFG, CLARGE,
     $                   CLARNV, CLATM1, CLASET, CSCAL, CSSCAL, SLATM1,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG, MAX, MOD
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
      IF( N.EQ.0 )
     $   RETURN
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
*     Decode RSIGN
*
      IF( LSAME( RSIGN, 'T' ) ) THEN
         IRSIGN = 1
      ELSE IF( LSAME( RSIGN, 'F' ) ) THEN
         IRSIGN = 0
      ELSE
         IRSIGN = -1
      END IF
*
*     Decode UPPER
*
      IF( LSAME( UPPER, 'T' ) ) THEN
         IUPPER = 1
      ELSE IF( LSAME( UPPER, 'F' ) ) THEN
         IUPPER = 0
      ELSE
         IUPPER = -1
      END IF
*
*     Decode SIM
*
      IF( LSAME( SIM, 'T' ) ) THEN
         ISIM = 1
      ELSE IF( LSAME( SIM, 'F' ) ) THEN
         ISIM = 0
      ELSE
         ISIM = -1
      END IF
*
*     Check DS, if MODES=0 and ISIM=1
*
      BADS = .FALSE.
      IF( MODES.EQ.0 .AND. ISIM.EQ.1 ) THEN
         DO 10 J = 1, N
            IF( DS( J ).EQ.ZERO )
     $         BADS = .TRUE.
   10    CONTINUE
      END IF
*
*     Set INFO if an error
*
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( IDIST.EQ.-1 ) THEN
         INFO = -2
      ELSE IF( ABS( MODE ).GT.6 ) THEN
         INFO = -5
      ELSE IF( ( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) .AND. COND.LT.ONE )
     $          THEN
         INFO = -6
      ELSE IF( IRSIGN.EQ.-1 ) THEN
         INFO = -9
      ELSE IF( IUPPER.EQ.-1 ) THEN
         INFO = -10
      ELSE IF( ISIM.EQ.-1 ) THEN
         INFO = -11
      ELSE IF( BADS ) THEN
         INFO = -12
      ELSE IF( ISIM.EQ.1 .AND. ABS( MODES ).GT.5 ) THEN
         INFO = -13
      ELSE IF( ISIM.EQ.1 .AND. MODES.NE.0 .AND. CONDS.LT.ONE ) THEN
         INFO = -14
      ELSE IF( KL.LT.1 ) THEN
         INFO = -15
      ELSE IF( KU.LT.1 .OR. ( KU.LT.N-1 .AND. KL.LT.N-1 ) ) THEN
         INFO = -16
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -19
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CLATME', -INFO )
         RETURN
      END IF
*
*     Initialize random number generator
*
      DO 20 I = 1, 4
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
   20 CONTINUE
*
      IF( MOD( ISEED( 4 ), 2 ).NE.1 )
     $   ISEED( 4 ) = ISEED( 4 ) + 1
*
*     2)      Set up diagonal of A
*
*             Compute D according to COND and MODE
*
      CALL CLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, N, IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
      IF( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) THEN
*
*        Scale by DMAX
*
         TEMP = ABS( D( 1 ) )
         DO 30 I = 2, N
            TEMP = MAX( TEMP, ABS( D( I ) ) )
   30    CONTINUE
*
         IF( TEMP.GT.ZERO ) THEN
            ALPHA = DMAX / TEMP
         ELSE
            INFO = 2
            RETURN
         END IF
*
         CALL CSCAL( N, ALPHA, D, 1 )
*
      END IF
*
      CALL CLASET( 'Full', N, N, CZERO, CZERO, A, LDA )
      CALL CCOPY( N, D, 1, A, LDA+1 )
*
*     3)      If UPPER='T', set upper triangle of A to random numbers.
*
      IF( IUPPER.NE.0 ) THEN
         DO 40 JC = 2, N
            CALL CLARNV( IDIST, ISEED, JC-1, A( 1, JC ) )
   40    CONTINUE
      END IF
*
*     4)      If SIM='T', apply similarity transformation.
*
*                                -1
*             Transform is  X A X  , where X = U S V, thus
*
*             it is  U S V A V' (1/S) U'
*
      IF( ISIM.NE.0 ) THEN
*
*        Compute S (singular values of the eigenvector matrix)
*        according to CONDS and MODES
*
         CALL SLATM1( MODES, CONDS, 0, 0, ISEED, DS, N, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 3
            RETURN
         END IF
*
*        Multiply by V and V'
*
         CALL CLARGE( N, A, LDA, ISEED, WORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 4
            RETURN
         END IF
*
*        Multiply by S and (1/S)
*
         DO 50 J = 1, N
            CALL CSSCAL( N, DS( J ), A( J, 1 ), LDA )
            IF( DS( J ).NE.ZERO ) THEN
               CALL CSSCAL( N, ONE / DS( J ), A( 1, J ), 1 )
            ELSE
               INFO = 5
               RETURN
            END IF
   50    CONTINUE
*
*        Multiply by U and U'
*
         CALL CLARGE( N, A, LDA, ISEED, WORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 4
            RETURN
         END IF
      END IF
*
*     5)      Reduce the bandwidth.
*
      IF( KL.LT.N-1 ) THEN
*
*        Reduce bandwidth -- kill column
*
         DO 60 JCR = KL + 1, N - 1
            IC = JCR - KL
            IROWS = N + 1 - JCR
            ICOLS = N + KL - JCR
*
            CALL CCOPY( IROWS, A( JCR, IC ), 1, WORK, 1 )
            XNORMS = WORK( 1 )
            CALL CLARFG( IROWS, XNORMS, WORK( 2 ), 1, TAU )
            TAU = CONJG( TAU )
            WORK( 1 ) = CONE
            ALPHA = CLARND( 5, ISEED )
*
            CALL CGEMV( 'C', IROWS, ICOLS, CONE, A( JCR, IC+1 ), LDA,
     $                  WORK, 1, CZERO, WORK( IROWS+1 ), 1 )
            CALL CGERC( IROWS, ICOLS, -TAU, WORK, 1, WORK( IROWS+1 ), 1,
     $                  A( JCR, IC+1 ), LDA )
*
            CALL CGEMV( 'N', N, IROWS, CONE, A( 1, JCR ), LDA, WORK, 1,
     $                  CZERO, WORK( IROWS+1 ), 1 )
            CALL CGERC( N, IROWS, -CONJG( TAU ), WORK( IROWS+1 ), 1,
     $                  WORK, 1, A( 1, JCR ), LDA )
*
            A( JCR, IC ) = XNORMS
            CALL CLASET( 'Full', IROWS-1, 1, CZERO, CZERO,
     $                   A( JCR+1, IC ), LDA )
*
            CALL CSCAL( ICOLS+1, ALPHA, A( JCR, IC ), LDA )
            CALL CSCAL( N, CONJG( ALPHA ), A( 1, JCR ), 1 )
   60    CONTINUE
      ELSE IF( KU.LT.N-1 ) THEN
*
*        Reduce upper bandwidth -- kill a row at a time.
*
         DO 70 JCR = KU + 1, N - 1
            IR = JCR - KU
            IROWS = N + KU - JCR
            ICOLS = N + 1 - JCR
*
            CALL CCOPY( ICOLS, A( IR, JCR ), LDA, WORK, 1 )
            XNORMS = WORK( 1 )
            CALL CLARFG( ICOLS, XNORMS, WORK( 2 ), 1, TAU )
            TAU = CONJG( TAU )
            WORK( 1 ) = CONE
            CALL CLACGV( ICOLS-1, WORK( 2 ), 1 )
            ALPHA = CLARND( 5, ISEED )
*
            CALL CGEMV( 'N', IROWS, ICOLS, CONE, A( IR+1, JCR ), LDA,
     $                  WORK, 1, CZERO, WORK( ICOLS+1 ), 1 )
            CALL CGERC( IROWS, ICOLS, -TAU, WORK( ICOLS+1 ), 1, WORK, 1,
     $                  A( IR+1, JCR ), LDA )
*
            CALL CGEMV( 'C', ICOLS, N, CONE, A( JCR, 1 ), LDA, WORK, 1,
     $                  CZERO, WORK( ICOLS+1 ), 1 )
            CALL CGERC( ICOLS, N, -CONJG( TAU ), WORK, 1,
     $                  WORK( ICOLS+1 ), 1, A( JCR, 1 ), LDA )
*
            A( IR, JCR ) = XNORMS
            CALL CLASET( 'Full', 1, ICOLS-1, CZERO, CZERO,
     $                   A( IR, JCR+1 ), LDA )
*
            CALL CSCAL( IROWS+1, ALPHA, A( IR, JCR ), 1 )
            CALL CSCAL( N, CONJG( ALPHA ), A( JCR, 1 ), LDA )
   70    CONTINUE
      END IF
*
*     Scale the matrix to have norm ANORM
*
      IF( ANORM.GE.ZERO ) THEN
         TEMP = CLANGE( 'M', N, N, A, LDA, TEMPA )
         IF( TEMP.GT.ZERO ) THEN
            RALPHA = ANORM / TEMP
            DO 80 J = 1, N
               CALL CSSCAL( N, RALPHA, A( 1, J ), 1 )
   80       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of CLATME
*
      END