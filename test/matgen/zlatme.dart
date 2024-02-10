      void zlatme(N, DIST, final Array<int> ISEED, D, MODE, COND, DMAX, RSIGN, UPPER, SIM, DS, MODES, CONDS, KL, KU, ANORM, final Matrix<double> A, final int LDA, final Array<double> _WORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             DIST, RSIGN, SIM, UPPER;
      int                INFO, KL, KU, LDA, MODE, MODES, N;
      double             ANORM, COND, CONDS;
      Complex         DMAX;
      int                ISEED( 4 );
      double             DS( * );
      Complex         A( LDA, * ), D( * ), WORK( * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      double             ONE;
      const              ONE = 1.0 ;
      Complex         CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      bool               BADS;
      int                I, IC, ICOLS, IDIST, IINFO, IR, IROWS, IRSIGN, ISIM, IUPPER, J, JC, JCR;
      double             RALPHA, TEMP;
      Complex         ALPHA, TAU, XNORMS;
      double             TEMPA( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             ZLANGE;
      //- Complex         ZLARND;
      // EXTERNAL lsame, ZLANGE, ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLATM1, XERBLA, ZCOPY, ZDSCAL, ZGEMV, ZGERC, ZLACGV, ZLARFG, ZLARGE, ZLARNV, ZLASET, ZLATM1, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCONJG, MAX, MOD

      // 1)      Decode and Test the input parameters.
      //         Initialize flags & seed.

      INFO = 0;

      // Quick return if possible

      if (N == 0) return;

      // Decode DIST

      if ( lsame( DIST, 'U' ) ) {
         IDIST = 1;
      } else if ( lsame( DIST, 'S' ) ) {
         IDIST = 2;
      } else if ( lsame( DIST, 'N' ) ) {
         IDIST = 3;
      } else if ( lsame( DIST, 'D' ) ) {
         IDIST = 4;
      } else {
         IDIST = -1;
      }

      // Decode RSIGN

      if ( lsame( RSIGN, 'T' ) ) {
         IRSIGN = 1;
      } else if ( lsame( RSIGN, 'F' ) ) {
         IRSIGN = 0;
      } else {
         IRSIGN = -1;
      }

      // Decode UPPER

      if ( lsame( UPPER, 'T' ) ) {
         IUPPER = 1;
      } else if ( lsame( UPPER, 'F' ) ) {
         IUPPER = 0;
      } else {
         IUPPER = -1;
      }

      // Decode SIM

      if ( lsame( SIM, 'T' ) ) {
         ISIM = 1;
      } else if ( lsame( SIM, 'F' ) ) {
         ISIM = 0;
      } else {
         ISIM = -1;
      }

      // Check DS, if MODES=0 and ISIM=1

      BADS = false;
      if ( MODES == 0 && ISIM == 1 ) {
         for (J = 1; J <= N; J++) { // 10
            if( DS( J ) == ZERO ) BADS = true;
         } // 10
      }

      // Set INFO if an error

      if ( N < 0 ) {
         INFO = -1;
      } else if ( IDIST == -1 ) {
         INFO = -2;
      } else if ( ( MODE ).abs() > 6 ) {
         INFO = -5;
      } else if ( ( MODE != 0 && ( MODE ).abs() != 6 ) && COND < ONE ) {
         INFO = -6;
      } else if ( IRSIGN == -1 ) {
         INFO = -9;
      } else if ( IUPPER == -1 ) {
         INFO = -10;
      } else if ( ISIM == -1 ) {
         INFO = -11;
      } else if ( BADS ) {
         INFO = -12;
      } else if ( ISIM == 1 && ( MODES ).abs() > 5 ) {
         INFO = -13;
      } else if ( ISIM == 1 && MODES != 0 && CONDS < ONE ) {
         INFO = -14;
      } else if ( KL < 1 ) {
         INFO = -15;
      } else if ( KU < 1 || ( KU < N-1 && KL < N-1 ) ) {
         INFO = -16;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -19;
      }

      if ( INFO != 0 ) {
         xerbla('ZLATME', -INFO );
         return;
      }

      // Initialize random number generator

      for (I = 1; I <= 4; I++) { // 20
         ISEED[I] = (( ISEED( I ) ).abs() % 4096);
      } // 20

      if( (ISEED( 4 ) % 2) != 1 ) ISEED( 4 ) = ISEED( 4 ) + 1;

      // 2)      Set up diagonal of A

              // Compute D according to COND and MODE

      zlatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, N, IINFO );
      if ( IINFO != 0 ) {
         INFO = 1;
         return;
      }
      if ( MODE != 0 && ( MODE ).abs() != 6 ) {

         // Scale by DMAX

         TEMP = ( D( 1 ) ).abs();
         for (I = 2; I <= N; I++) { // 30
            TEMP = max( TEMP, ( D( I ) ).abs() );
         } // 30

         if ( TEMP > ZERO ) {
            ALPHA = DMAX / TEMP;
         } else {
            INFO = 2;
            return;
         }

         zscal(N, ALPHA, D, 1 );

      }

      zlaset('Full', N, N, CZERO, CZERO, A, LDA );
      zcopy(N, D, 1, A, LDA+1 );

      // 3)      If UPPER='T', set upper triangle of A to random numbers.

      if ( IUPPER != 0 ) {
         for (JC = 2; JC <= N; JC++) { // 40
            zlarnv(IDIST, ISEED, JC-1, A( 1, JC ) );
         } // 40
      }

      // 4)      If SIM='T', apply similarity transformation.

                                 // -1
              // Transform is  X A X  , where X = U S V, thus

              // it is  U S V A V' (1/S) U'

      if ( ISIM != 0 ) {

         // Compute S (singular values of the eigenvector matrix)
         // according to CONDS and MODES

         dlatm1(MODES, CONDS, 0, 0, ISEED, DS, N, IINFO );
         if ( IINFO != 0 ) {
            INFO = 3;
            return;
         }

         // Multiply by V and V'

         zlarge(N, A, LDA, ISEED, WORK, IINFO );
         if ( IINFO != 0 ) {
            INFO = 4;
            return;
         }

         // Multiply by S and (1/S)

         for (J = 1; J <= N; J++) { // 50
            zdscal(N, DS( J ), A( J, 1 ), LDA );
            if ( DS( J ) != ZERO ) {
               zdscal(N, ONE / DS( J ), A( 1, J ), 1 );
            } else {
               INFO = 5;
               return;
            }
         } // 50

         // Multiply by U and U'

         zlarge(N, A, LDA, ISEED, WORK, IINFO );
         if ( IINFO != 0 ) {
            INFO = 4;
            return;
         }
      }

      // 5)      Reduce the bandwidth.

      if ( KL < N-1 ) {

         // Reduce bandwidth -- kill column

         for (JCR = KL + 1; JCR <= N - 1; JCR++) { // 60
            IC = JCR - KL;
            IROWS = N + 1 - JCR;
            ICOLS = N + KL - JCR;

            zcopy(IROWS, A( JCR, IC ), 1, WORK, 1 );
            XNORMS = WORK( 1 );
            zlarfg(IROWS, XNORMS, WORK( 2 ), 1, TAU );
            TAU = DCONJG( TAU );
            WORK[1] = CONE;
            ALPHA = ZLARND( 5, ISEED );

            zgemv('C', IROWS, ICOLS, CONE, A( JCR, IC+1 ), LDA, WORK, 1, CZERO, WORK( IROWS+1 ), 1 );
            zgerc(IROWS, ICOLS, -TAU, WORK, 1, WORK( IROWS+1 ), 1, A( JCR, IC+1 ), LDA );

            zgemv('N', N, IROWS, CONE, A( 1, JCR ), LDA, WORK, 1, CZERO, WORK( IROWS+1 ), 1 );
            zgerc(N, IROWS, -DCONJG( TAU ), WORK( IROWS+1 ), 1, WORK, 1, A( 1, JCR ), LDA );

            A[JCR][IC] = XNORMS;
            zlaset('Full', IROWS-1, 1, CZERO, CZERO, A( JCR+1, IC ), LDA );

            zscal(ICOLS+1, ALPHA, A( JCR, IC ), LDA );
            zscal(N, DCONJG( ALPHA ), A( 1, JCR ), 1 );
         } // 60
      } else if ( KU < N-1 ) {

         // Reduce upper bandwidth -- kill a row at a time.

         for (JCR = KU + 1; JCR <= N - 1; JCR++) { // 70
            IR = JCR - KU;
            IROWS = N + KU - JCR;
            ICOLS = N + 1 - JCR;

            zcopy(ICOLS, A( IR, JCR ), LDA, WORK, 1 );
            XNORMS = WORK( 1 );
            zlarfg(ICOLS, XNORMS, WORK( 2 ), 1, TAU );
            TAU = DCONJG( TAU );
            WORK[1] = CONE;
            zlacgv(ICOLS-1, WORK( 2 ), 1 );
            ALPHA = ZLARND( 5, ISEED );

            zgemv('N', IROWS, ICOLS, CONE, A( IR+1, JCR ), LDA, WORK, 1, CZERO, WORK( ICOLS+1 ), 1 );
            zgerc(IROWS, ICOLS, -TAU, WORK( ICOLS+1 ), 1, WORK, 1, A( IR+1, JCR ), LDA );

            zgemv('C', ICOLS, N, CONE, A( JCR, 1 ), LDA, WORK, 1, CZERO, WORK( ICOLS+1 ), 1 );
            zgerc(ICOLS, N, -DCONJG( TAU ), WORK, 1, WORK( ICOLS+1 ), 1, A( JCR, 1 ), LDA );

            A[IR][JCR] = XNORMS;
            zlaset('Full', 1, ICOLS-1, CZERO, CZERO, A( IR, JCR+1 ), LDA );

            zscal(IROWS+1, ALPHA, A( IR, JCR ), 1 );
            zscal(N, DCONJG( ALPHA ), A( JCR, 1 ), LDA );
         } // 70
      }

      // Scale the matrix to have norm ANORM

      if ( ANORM >= ZERO ) {
         TEMP = ZLANGE( 'M', N, N, A, LDA, TEMPA );
         if ( TEMP > ZERO ) {
            RALPHA = ANORM / TEMP;
            for (J = 1; J <= N; J++) { // 80
               zdscal(N, RALPHA, A( 1, J ), 1 );
            } // 80
         }
      }

      }
