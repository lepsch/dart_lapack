      SUBROUTINE SLATME( N, DIST, ISEED, D, MODE, COND, DMAX, EI, RSIGN, UPPER, SIM, DS, MODES, CONDS, KL, KU, ANORM, A, LDA, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIST, RSIGN, SIM, UPPER;
      int                INFO, KL, KU, LDA, MODE, MODES, N;
      REAL               ANORM, COND, CONDS, DMAX
      // ..
      // .. Array Arguments ..
      String             EI( * );
      int                ISEED( 4 );
      REAL               A( LDA, * ), D( * ), DS( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO
      const              ZERO = 0.0E0 ;
      REAL               ONE
      const              ONE = 1.0E0 ;
      REAL               HALF
      const              HALF = 1.0E0 / 2.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               BADEI, BADS, USEEI;
      int                I, IC, ICOLS, IDIST, IINFO, IR, IROWS, IRSIGN, ISIM, IUPPER, J, JC, JCR, JR;
      REAL               ALPHA, TAU, TEMP, XNORMS
      // ..
      // .. Local Arrays ..
      REAL               TEMPA( 1 )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLANGE, SLARAN
      // EXTERNAL LSAME, SLANGE, SLARAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMV, SGER, SLARFG, SLARGE, SLARNV, SLATM1, SLASET, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MOD
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
      } else {
         IDIST = -1
      }

      // Check EI

      USEEI = .TRUE.
      BADEI = .FALSE.
      if ( LSAME( EI( 1 ), ' ' ) .OR. MODE.NE.0 ) {
         USEEI = .FALSE.
      } else {
         if ( LSAME( EI( 1 ), 'R' ) ) {
            for (J = 2; J <= N; J++) { // 10
               if ( LSAME( EI( J ), 'I' ) ) {
                  IF( LSAME( EI( J-1 ), 'I' ) ) BADEI = .TRUE.
               } else {
                  IF( .NOT.LSAME( EI( J ), 'R' ) ) BADEI = .TRUE.
               }
            } // 10
         } else {
            BADEI = .TRUE.
         }
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
         for (J = 1; J <= N; J++) { // 20
            IF( DS( J ).EQ.ZERO ) BADS = .TRUE.
         } // 20
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
      } else if ( BADEI ) {
         INFO = -8
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
         xerbla('SLATME', -INFO );
         RETURN
      }

      // Initialize random number generator

      for (I = 1; I <= 4; I++) { // 30
         ISEED( I ) = MOD( ABS( ISEED( I ) ), 4096 )
      } // 30

      IF( MOD( ISEED( 4 ), 2 ).NE.1 ) ISEED( 4 ) = ISEED( 4 ) + 1

      // 2)      Set up diagonal of A

              // Compute D according to COND and MODE

      slatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, N, IINFO );
      if ( IINFO.NE.0 ) {
         INFO = 1
         RETURN
      }
      if ( MODE.NE.0 .AND. ABS( MODE ).NE.6 ) {

         // Scale by DMAX

         TEMP = ABS( D( 1 ) )
         for (I = 2; I <= N; I++) { // 40
            TEMP = MAX( TEMP, ABS( D( I ) ) )
         } // 40

         if ( TEMP.GT.ZERO ) {
            ALPHA = DMAX / TEMP
         } else if ( DMAX.NE.ZERO ) {
            INFO = 2
            RETURN
         } else {
            ALPHA = ZERO
         }

         sscal(N, ALPHA, D, 1 );

      }

      slaset('Full', N, N, ZERO, ZERO, A, LDA );
      scopy(N, D, 1, A, LDA+1 );

      // Set up complex conjugate pairs

      if ( MODE.EQ.0 ) {
         if ( USEEI ) {
            for (J = 2; J <= N; J++) { // 50
               if ( LSAME( EI( J ), 'I' ) ) {
                  A( J-1, J ) = A( J, J )
                  A( J, J-1 ) = -A( J, J )
                  A( J, J ) = A( J-1, J-1 )
               }
            } // 50
         }

      } else if ( ABS( MODE ).EQ.5 ) {

         DO 60 J = 2, N, 2
            if ( SLARAN( ISEED ).GT.HALF ) {
               A( J-1, J ) = A( J, J )
               A( J, J-1 ) = -A( J, J )
               A( J, J ) = A( J-1, J-1 )
            }
         } // 60
      }

      // 3)      If UPPER='T', set upper triangle of A to random numbers.
              // (but don't modify the corners of 2x2 blocks.)

      if ( IUPPER.NE.0 ) {
         for (JC = 2; JC <= N; JC++) { // 70
            if ( A( JC-1, JC ).NE.ZERO ) {
               JR = JC - 2
            } else {
               JR = JC - 1
            }
            slarnv(IDIST, ISEED, JR, A( 1, JC ) );
         } // 70
      }

      // 4)      If SIM='T', apply similarity transformation.

                                 // -1
              // Transform is  X A X  , where X = U S V, thus

              // it is  U S V A V' (1/S) U'

      if ( ISIM.NE.0 ) {

         // Compute S (singular values of the eigenvector matrix)
         // according to CONDS and MODES

         slatm1(MODES, CONDS, 0, 0, ISEED, DS, N, IINFO );
         if ( IINFO.NE.0 ) {
            INFO = 3
            RETURN
         }

         // Multiply by V and V'

         slarge(N, A, LDA, ISEED, WORK, IINFO );
         if ( IINFO.NE.0 ) {
            INFO = 4
            RETURN
         }

         // Multiply by S and (1/S)

         for (J = 1; J <= N; J++) { // 80
            sscal(N, DS( J ), A( J, 1 ), LDA );
            if ( DS( J ).NE.ZERO ) {
               sscal(N, ONE / DS( J ), A( 1, J ), 1 );
            } else {
               INFO = 5
               RETURN
            }
         } // 80

         // Multiply by U and U'

         slarge(N, A, LDA, ISEED, WORK, IINFO );
         if ( IINFO.NE.0 ) {
            INFO = 4
            RETURN
         }
      }

      // 5)      Reduce the bandwidth.

      if ( KL.LT.N-1 ) {

         // Reduce bandwidth -- kill column

         DO 90 JCR = KL + 1, N - 1
            IC = JCR - KL
            IROWS = N + 1 - JCR
            ICOLS = N + KL - JCR

            scopy(IROWS, A( JCR, IC ), 1, WORK, 1 );
            XNORMS = WORK( 1 )
            slarfg(IROWS, XNORMS, WORK( 2 ), 1, TAU );
            WORK( 1 ) = ONE

            sgemv('T', IROWS, ICOLS, ONE, A( JCR, IC+1 ), LDA, WORK, 1, ZERO, WORK( IROWS+1 ), 1 )             CALL SGER( IROWS, ICOLS, -TAU, WORK, 1, WORK( IROWS+1 ), 1, A( JCR, IC+1 ), LDA );

            sgemv('N', N, IROWS, ONE, A( 1, JCR ), LDA, WORK, 1, ZERO, WORK( IROWS+1 ), 1 )             CALL SGER( N, IROWS, -TAU, WORK( IROWS+1 ), 1, WORK, 1, A( 1, JCR ), LDA );

            A( JCR, IC ) = XNORMS
            slaset('Full', IROWS-1, 1, ZERO, ZERO, A( JCR+1, IC ), LDA );
         } // 90
      } else if ( KU.LT.N-1 ) {

         // Reduce upper bandwidth -- kill a row at a time.

         DO 100 JCR = KU + 1, N - 1
            IR = JCR - KU
            IROWS = N + KU - JCR
            ICOLS = N + 1 - JCR

            scopy(ICOLS, A( IR, JCR ), LDA, WORK, 1 );
            XNORMS = WORK( 1 )
            slarfg(ICOLS, XNORMS, WORK( 2 ), 1, TAU );
            WORK( 1 ) = ONE

            sgemv('N', IROWS, ICOLS, ONE, A( IR+1, JCR ), LDA, WORK, 1, ZERO, WORK( ICOLS+1 ), 1 )             CALL SGER( IROWS, ICOLS, -TAU, WORK( ICOLS+1 ), 1, WORK, 1, A( IR+1, JCR ), LDA );

            sgemv('C', ICOLS, N, ONE, A( JCR, 1 ), LDA, WORK, 1, ZERO, WORK( ICOLS+1 ), 1 )             CALL SGER( ICOLS, N, -TAU, WORK, 1, WORK( ICOLS+1 ), 1, A( JCR, 1 ), LDA );

            A( IR, JCR ) = XNORMS
            slaset('Full', 1, ICOLS-1, ZERO, ZERO, A( IR, JCR+1 ), LDA );
         } // 100
      }

      // Scale the matrix to have norm ANORM

      if ( ANORM.GE.ZERO ) {
         TEMP = SLANGE( 'M', N, N, A, LDA, TEMPA )
         if ( TEMP.GT.ZERO ) {
            ALPHA = ANORM / TEMP
            for (J = 1; J <= N; J++) { // 110
               sscal(N, ALPHA, A( 1, J ), 1 );
            } // 110
         }
      }

      RETURN

      // End of SLATME

      }
