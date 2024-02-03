      SUBROUTINE DLATME( N, DIST, ISEED, D, MODE, COND, DMAX, EI, RSIGN, UPPER, SIM, DS, MODES, CONDS, KL, KU, ANORM, A, LDA, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             DIST, RSIGN, SIM, UPPER;
      int                INFO, KL, KU, LDA, MODE, MODES, N
      DOUBLE PRECISION   ANORM, COND, CONDS, DMAX
*     ..
*     .. Array Arguments ..
      String             EI( * );
      int                ISEED( 4 )
      DOUBLE PRECISION   A( LDA, * ), D( * ), DS( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 1.0D0 / 2.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BADEI, BADS, USEEI
      int                I, IC, ICOLS, IDIST, IINFO, IR, IROWS, IRSIGN, ISIM, IUPPER, J, JC, JCR, JR
      DOUBLE PRECISION   ALPHA, TAU, TEMP, XNORMS
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   TEMPA( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLANGE, DLARAN
      EXTERNAL           LSAME, DLANGE, DLARAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMV, DGER, DLARFG, DLARGE, DLARNV, DLASET, DLATM1, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MOD
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
      IF( N.EQ.0 ) RETURN
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
*     Check EI
*
      USEEI = .TRUE.
      BADEI = .FALSE.
      IF( LSAME( EI( 1 ), ' ' ) .OR. MODE.NE.0 ) THEN
         USEEI = .FALSE.
      ELSE
         IF( LSAME( EI( 1 ), 'R' ) ) THEN
            DO 10 J = 2, N
               IF( LSAME( EI( J ), 'I' ) ) THEN
                  IF( LSAME( EI( J-1 ), 'I' ) ) BADEI = .TRUE.
               ELSE
                  IF( .NOT.LSAME( EI( J ), 'R' ) ) BADEI = .TRUE.
               END IF
   10       CONTINUE
         ELSE
            BADEI = .TRUE.
         END IF
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
         DO 20 J = 1, N
            IF( DS( J ).EQ.ZERO ) BADS = .TRUE.
   20    CONTINUE
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
      ELSE IF( ( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) .AND. COND.LT.ONE ) THEN
         INFO = -6
      ELSE IF( BADEI ) THEN
         INFO = -8
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
         CALL XERBLA( 'DLATME', -INFO )
         RETURN
      END IF
*
*     Initialize random number generator
*
      DO 30 I = 1, 4
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
   30 CONTINUE
*
      IF( MOD( ISEED( 4 ), 2 ).NE.1 ) ISEED( 4 ) = ISEED( 4 ) + 1
*
*     2)      Set up diagonal of A
*
*             Compute D according to COND and MODE
*
      CALL DLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, N, IINFO )
      IF( IINFO.NE.0 ) THEN
         INFO = 1
         RETURN
      END IF
      IF( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) THEN
*
*        Scale by DMAX
*
         TEMP = ABS( D( 1 ) )
         DO 40 I = 2, N
            TEMP = MAX( TEMP, ABS( D( I ) ) )
   40    CONTINUE
*
         IF( TEMP.GT.ZERO ) THEN
            ALPHA = DMAX / TEMP
         ELSE IF( DMAX.NE.ZERO ) THEN
            INFO = 2
            RETURN
         ELSE
            ALPHA = ZERO
         END IF
*
         CALL DSCAL( N, ALPHA, D, 1 )
*
      END IF
*
      CALL DLASET( 'Full', N, N, ZERO, ZERO, A, LDA )
      CALL DCOPY( N, D, 1, A, LDA+1 )
*
*     Set up complex conjugate pairs
*
      IF( MODE.EQ.0 ) THEN
         IF( USEEI ) THEN
            DO 50 J = 2, N
               IF( LSAME( EI( J ), 'I' ) ) THEN
                  A( J-1, J ) = A( J, J )
                  A( J, J-1 ) = -A( J, J )
                  A( J, J ) = A( J-1, J-1 )
               END IF
   50       CONTINUE
         END IF
*
      ELSE IF( ABS( MODE ).EQ.5 ) THEN
*
         DO 60 J = 2, N, 2
            IF( DLARAN( ISEED ).GT.HALF ) THEN
               A( J-1, J ) = A( J, J )
               A( J, J-1 ) = -A( J, J )
               A( J, J ) = A( J-1, J-1 )
            END IF
   60    CONTINUE
      END IF
*
*     3)      If UPPER='T', set upper triangle of A to random numbers.
*             (but don't modify the corners of 2x2 blocks.)
*
      IF( IUPPER.NE.0 ) THEN
         DO 70 JC = 2, N
            IF( A( JC-1, JC ).NE.ZERO ) THEN
               JR = JC - 2
            ELSE
               JR = JC - 1
            END IF
            CALL DLARNV( IDIST, ISEED, JR, A( 1, JC ) )
   70    CONTINUE
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
         CALL DLATM1( MODES, CONDS, 0, 0, ISEED, DS, N, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 3
            RETURN
         END IF
*
*        Multiply by V and V'
*
         CALL DLARGE( N, A, LDA, ISEED, WORK, IINFO )
         IF( IINFO.NE.0 ) THEN
            INFO = 4
            RETURN
         END IF
*
*        Multiply by S and (1/S)
*
         DO 80 J = 1, N
            CALL DSCAL( N, DS( J ), A( J, 1 ), LDA )
            IF( DS( J ).NE.ZERO ) THEN
               CALL DSCAL( N, ONE / DS( J ), A( 1, J ), 1 )
            ELSE
               INFO = 5
               RETURN
            END IF
   80    CONTINUE
*
*        Multiply by U and U'
*
         CALL DLARGE( N, A, LDA, ISEED, WORK, IINFO )
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
         DO 90 JCR = KL + 1, N - 1
            IC = JCR - KL
            IROWS = N + 1 - JCR
            ICOLS = N + KL - JCR
*
            CALL DCOPY( IROWS, A( JCR, IC ), 1, WORK, 1 )
            XNORMS = WORK( 1 )
            CALL DLARFG( IROWS, XNORMS, WORK( 2 ), 1, TAU )
            WORK( 1 ) = ONE
*
            CALL DGEMV( 'T', IROWS, ICOLS, ONE, A( JCR, IC+1 ), LDA, WORK, 1, ZERO, WORK( IROWS+1 ), 1 )             CALL DGER( IROWS, ICOLS, -TAU, WORK, 1, WORK( IROWS+1 ), 1, A( JCR, IC+1 ), LDA )
*
            CALL DGEMV( 'N', N, IROWS, ONE, A( 1, JCR ), LDA, WORK, 1, ZERO, WORK( IROWS+1 ), 1 )             CALL DGER( N, IROWS, -TAU, WORK( IROWS+1 ), 1, WORK, 1, A( 1, JCR ), LDA )
*
            A( JCR, IC ) = XNORMS
            CALL DLASET( 'Full', IROWS-1, 1, ZERO, ZERO, A( JCR+1, IC ), LDA )
   90    CONTINUE
      ELSE IF( KU.LT.N-1 ) THEN
*
*        Reduce upper bandwidth -- kill a row at a time.
*
         DO 100 JCR = KU + 1, N - 1
            IR = JCR - KU
            IROWS = N + KU - JCR
            ICOLS = N + 1 - JCR
*
            CALL DCOPY( ICOLS, A( IR, JCR ), LDA, WORK, 1 )
            XNORMS = WORK( 1 )
            CALL DLARFG( ICOLS, XNORMS, WORK( 2 ), 1, TAU )
            WORK( 1 ) = ONE
*
            CALL DGEMV( 'N', IROWS, ICOLS, ONE, A( IR+1, JCR ), LDA, WORK, 1, ZERO, WORK( ICOLS+1 ), 1 )             CALL DGER( IROWS, ICOLS, -TAU, WORK( ICOLS+1 ), 1, WORK, 1, A( IR+1, JCR ), LDA )
*
            CALL DGEMV( 'C', ICOLS, N, ONE, A( JCR, 1 ), LDA, WORK, 1, ZERO, WORK( ICOLS+1 ), 1 )             CALL DGER( ICOLS, N, -TAU, WORK, 1, WORK( ICOLS+1 ), 1, A( JCR, 1 ), LDA )
*
            A( IR, JCR ) = XNORMS
            CALL DLASET( 'Full', 1, ICOLS-1, ZERO, ZERO, A( IR, JCR+1 ), LDA )
  100    CONTINUE
      END IF
*
*     Scale the matrix to have norm ANORM
*
      IF( ANORM.GE.ZERO ) THEN
         TEMP = DLANGE( 'M', N, N, A, LDA, TEMPA )
         IF( TEMP.GT.ZERO ) THEN
            ALPHA = ANORM / TEMP
            DO 110 J = 1, N
               CALL DSCAL( N, ALPHA, A( 1, J ), 1 )
  110       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DLATME
*
      END
