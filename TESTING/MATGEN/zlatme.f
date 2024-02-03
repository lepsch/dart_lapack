      SUBROUTINE ZLATME( N, DIST, ISEED, D, MODE, COND, DMAX, RSIGN, UPPER, SIM, DS, MODES, CONDS, KL, KU, ANORM, A, LDA, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIST, RSIGN, SIM, UPPER;
      int                INFO, KL, KU, LDA, MODE, MODES, N;
      double             ANORM, COND, CONDS;
      COMPLEX*16         DMAX
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             DS( * );
      COMPLEX*16         A( LDA, * ), D( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      double             ONE;
      const              ONE = 1.0D+0 ;
      COMPLEX*16         CZERO
      const              CZERO = ( 0.0D+0, 0.0D+0 ) ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               BADS;
      int                I, IC, ICOLS, IDIST, IINFO, IR, IROWS, IRSIGN, ISIM, IUPPER, J, JC, JCR;
      double             RALPHA, TEMP;
      COMPLEX*16         ALPHA, TAU, XNORMS
      // ..
      // .. Local Arrays ..
      double             TEMPA( 1 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             ZLANGE;
      COMPLEX*16         ZLARND
      // EXTERNAL LSAME, ZLANGE, ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLATM1, XERBLA, ZCOPY, ZDSCAL, ZGEMV, ZGERC, ZLACGV, ZLARFG, ZLARGE, ZLARNV, ZLASET, ZLATM1, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX, MOD
      // ..
      // .. Executable Statements ..

      // 1)      Decode and Test the input parameters.
              // Initialize flags & seed.

      INFO = 0

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Decode DIST

      if ( LSAME( DIST, 'U' ) ) {
         IDIST = 1
      } else if ( LSAME( DIST, 'S' ) ) {
         IDIST = 2
      } else if ( LSAME( DIST, 'N' ) ) {
         IDIST = 3
      } else if ( LSAME( DIST, 'D' ) ) {
         IDIST = 4
      } else {
         IDIST = -1
      }

      // Decode RSIGN

      if ( LSAME( RSIGN, 'T' ) ) {
         IRSIGN = 1
      } else if ( LSAME( RSIGN, 'F' ) ) {
         IRSIGN = 0
      } else {
         IRSIGN = -1
      }

      // Decode UPPER

      if ( LSAME( UPPER, 'T' ) ) {
         IUPPER = 1
      } else if ( LSAME( UPPER, 'F' ) ) {
         IUPPER = 0
      } else {
         IUPPER = -1
      }

      // Decode SIM

      if ( LSAME( SIM, 'T' ) ) {
         ISIM = 1
      } else if ( LSAME( SIM, 'F' ) ) {
         ISIM = 0
      } else {
         ISIM = -1
      }

      // Check DS, if MODES=0 and ISIM=1

      BADS = .FALSE.
      if ( MODES.EQ.0 .AND. ISIM.EQ.1 ) {
         DO 10 J = 1, N
            IF( DS( J ).EQ.ZERO ) BADS = .TRUE.
   10    CONTINUE
      }

      // Set INFO if an error

      if ( N.LT.0 ) {
         INFO = -1
      } else if ( IDIST.EQ.-1 ) {
         INFO = -2
      } else if ( ABS( MODE ).GT.6 ) {
         INFO = -5
      } else if ( ( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) .AND. COND.LT.ONE ) {
         INFO = -6
      } else if ( IRSIGN.EQ.-1 ) {
         INFO = -9
      } else if ( IUPPER.EQ.-1 ) {
         INFO = -10
      } else if ( ISIM.EQ.-1 ) {
         INFO = -11
      } else if ( BADS ) {
         INFO = -12
      } else if ( ISIM.EQ.1 .AND. ABS( MODES ).GT.5 ) {
         INFO = -13
      } else if ( ISIM.EQ.1 .AND. MODES.NE.0 .AND. CONDS.LT.ONE ) {
         INFO = -14
      } else if ( KL.LT.1 ) {
         INFO = -15
      } else if ( KU.LT.1 .OR. ( KU.LT.N-1 .AND. KL.LT.N-1 ) ) {
         INFO = -16
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -19
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZLATME', -INFO )
         RETURN
      }

      // Initialize random number generator

      DO 20 I = 1, 4
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
   20 CONTINUE

      IF( MOD( ISEED( 4 ), 2 ).NE.1 ) ISEED( 4 ) = ISEED( 4 ) + 1

      // 2)      Set up diagonal of A

              // Compute D according to COND and MODE

      CALL ZLATM1( MODE, COND, IRSIGN, IDIST, ISEED, D, N, IINFO )
      if ( IINFO.NE.0 ) {
         INFO = 1
         RETURN
      }
      if ( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) {

         // Scale by DMAX

         TEMP = ABS( D( 1 ) )
         DO 30 I = 2, N
            TEMP = MAX( TEMP, ABS( D( I ) ) )
   30    CONTINUE

         if ( TEMP.GT.ZERO ) {
            ALPHA = DMAX / TEMP
         } else {
            INFO = 2
            RETURN
         }

         CALL ZSCAL( N, ALPHA, D, 1 )

      }

      CALL ZLASET( 'Full', N, N, CZERO, CZERO, A, LDA )
      CALL ZCOPY( N, D, 1, A, LDA+1 )

      // 3)      If UPPER='T', set upper triangle of A to random numbers.

      if ( IUPPER.NE.0 ) {
         DO 40 JC = 2, N
            CALL ZLARNV( IDIST, ISEED, JC-1, A( 1, JC ) )
   40    CONTINUE
      }

      // 4)      If SIM='T', apply similarity transformation.

                                 // -1
              // Transform is  X A X  , where X = U S V, thus

              // it is  U S V A V' (1/S) U'

      if ( ISIM.NE.0 ) {

         // Compute S (singular values of the eigenvector matrix)
         // according to CONDS and MODES

         CALL DLATM1( MODES, CONDS, 0, 0, ISEED, DS, N, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = 3
            RETURN
         }

         // Multiply by V and V'

         CALL ZLARGE( N, A, LDA, ISEED, WORK, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = 4
            RETURN
         }

         // Multiply by S and (1/S)

         DO 50 J = 1, N
            CALL ZDSCAL( N, DS( J ), A( J, 1 ), LDA )
            if ( DS( J ).NE.ZERO ) {
               CALL ZDSCAL( N, ONE / DS( J ), A( 1, J ), 1 )
            } else {
               INFO = 5
               RETURN
            }
   50    CONTINUE

         // Multiply by U and U'

         CALL ZLARGE( N, A, LDA, ISEED, WORK, IINFO )
         if ( IINFO.NE.0 ) {
            INFO = 4
            RETURN
         }
      }

      // 5)      Reduce the bandwidth.

      if ( KL.LT.N-1 ) {

         // Reduce bandwidth -- kill column

         DO 60 JCR = KL + 1, N - 1
            IC = JCR - KL
            IROWS = N + 1 - JCR
            ICOLS = N + KL - JCR

            CALL ZCOPY( IROWS, A( JCR, IC ), 1, WORK, 1 )
            XNORMS = WORK( 1 )
            CALL ZLARFG( IROWS, XNORMS, WORK( 2 ), 1, TAU )
            TAU = DCONJG( TAU )
            WORK( 1 ) = CONE
            ALPHA = ZLARND( 5, ISEED )

            CALL ZGEMV( 'C', IROWS, ICOLS, CONE, A( JCR, IC+1 ), LDA, WORK, 1, CZERO, WORK( IROWS+1 ), 1 )             CALL ZGERC( IROWS, ICOLS, -TAU, WORK, 1, WORK( IROWS+1 ), 1, A( JCR, IC+1 ), LDA )

            CALL ZGEMV( 'N', N, IROWS, CONE, A( 1, JCR ), LDA, WORK, 1, CZERO, WORK( IROWS+1 ), 1 )             CALL ZGERC( N, IROWS, -DCONJG( TAU ), WORK( IROWS+1 ), 1, WORK, 1, A( 1, JCR ), LDA )

            A( JCR, IC ) = XNORMS
            CALL ZLASET( 'Full', IROWS-1, 1, CZERO, CZERO, A( JCR+1, IC ), LDA )

            CALL ZSCAL( ICOLS+1, ALPHA, A( JCR, IC ), LDA )
            CALL ZSCAL( N, DCONJG( ALPHA ), A( 1, JCR ), 1 )
   60    CONTINUE
      } else if ( KU.LT.N-1 ) {

         // Reduce upper bandwidth -- kill a row at a time.

         DO 70 JCR = KU + 1, N - 1
            IR = JCR - KU
            IROWS = N + KU - JCR
            ICOLS = N + 1 - JCR

            CALL ZCOPY( ICOLS, A( IR, JCR ), LDA, WORK, 1 )
            XNORMS = WORK( 1 )
            CALL ZLARFG( ICOLS, XNORMS, WORK( 2 ), 1, TAU )
            TAU = DCONJG( TAU )
            WORK( 1 ) = CONE
            CALL ZLACGV( ICOLS-1, WORK( 2 ), 1 )
            ALPHA = ZLARND( 5, ISEED )

            CALL ZGEMV( 'N', IROWS, ICOLS, CONE, A( IR+1, JCR ), LDA, WORK, 1, CZERO, WORK( ICOLS+1 ), 1 )             CALL ZGERC( IROWS, ICOLS, -TAU, WORK( ICOLS+1 ), 1, WORK, 1, A( IR+1, JCR ), LDA )

            CALL ZGEMV( 'C', ICOLS, N, CONE, A( JCR, 1 ), LDA, WORK, 1, CZERO, WORK( ICOLS+1 ), 1 )             CALL ZGERC( ICOLS, N, -DCONJG( TAU ), WORK, 1, WORK( ICOLS+1 ), 1, A( JCR, 1 ), LDA )

            A( IR, JCR ) = XNORMS
            CALL ZLASET( 'Full', 1, ICOLS-1, CZERO, CZERO, A( IR, JCR+1 ), LDA )

            CALL ZSCAL( IROWS+1, ALPHA, A( IR, JCR ), 1 )
            CALL ZSCAL( N, DCONJG( ALPHA ), A( JCR, 1 ), LDA )
   70    CONTINUE
      }

      // Scale the matrix to have norm ANORM

      if ( ANORM.GE.ZERO ) {
         TEMP = ZLANGE( 'M', N, N, A, LDA, TEMPA )
         if ( TEMP.GT.ZERO ) {
            RALPHA = ANORM / TEMP
            DO 80 J = 1, N
               CALL ZDSCAL( N, RALPHA, A( 1, J ), 1 )
   80       CONTINUE
         }
      }

      RETURN

      // End of ZLATME

      }
