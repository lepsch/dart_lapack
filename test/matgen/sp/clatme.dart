      void clatme(N, DIST, ISEED, D, MODE, COND, DMAX, RSIGN, UPPER, SIM, DS, MODES, CONDS, KL, KU, ANORM, A, LDA, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIST, RSIGN, SIM, UPPER;
      int                INFO, KL, KU, LDA, MODE, MODES, N;
      double               ANORM, COND, CONDS;
      Complex            DMAX;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double               DS( * );
      Complex            A( LDA, * ), D( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO;
      const              ZERO = 0.0 ;
      double               ONE;
      const              ONE = 1.0 ;
      Complex            CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               BADS;
      int                I, IC, ICOLS, IDIST, IINFO, IR, IROWS, IRSIGN, ISIM, IUPPER, J, JC, JCR;
      double               RALPHA, TEMP;
      Complex            ALPHA, TAU, XNORMS;
      // ..
      // .. Local Arrays ..
      double               TEMPA( 1 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANGE;
      //- COMPLEX            CLARND;
      // EXTERNAL lsame, CLANGE, CLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMV, CGERC, CLACGV, CLARFG, CLARGE, CLARNV, CLATM1, CLASET, CSCAL, CSSCAL, SLATM1, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, MOD
      // ..
      // .. Executable Statements ..

      // 1)      Decode and Test the input parameters.
              // Initialize flags & seed.

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
         xerbla('CLATME', -INFO );
         return;
      }

      // Initialize random number generator

      for (I = 1; I <= 4; I++) { // 20
         ISEED[I] = (( ISEED( I ) ).abs() % 4096);
      } // 20

      if( (ISEED( 4 ) % 2) != 1 ) ISEED( 4 ) = ISEED( 4 ) + 1;

      // 2)      Set up diagonal of A

              // Compute D according to COND and MODE

      clatm1(MODE, COND, IRSIGN, IDIST, ISEED, D, N, IINFO );
      if ( IINFO != 0 ) {
         INFO = 1;
         return;
      }
      if ( MODE != 0 && ( MODE ).abs() != 6 ) {

         // Scale by DMAX

         TEMP = ( D( 1 ) ).abs();
         for (I = 2; I <= N; I++) { // 30
            TEMP = max( TEMP, ( D( I ) ) ).abs();
         } // 30

         if ( TEMP > ZERO ) {
            ALPHA = DMAX / TEMP;
         } else {
            INFO = 2;
            return;
         }

         cscal(N, ALPHA, D, 1 );

      }

      claset('Full', N, N, CZERO, CZERO, A, LDA );
      ccopy(N, D, 1, A, LDA+1 );

      // 3)      If UPPER='T', set upper triangle of A to random numbers.

      if ( IUPPER != 0 ) {
         for (JC = 2; JC <= N; JC++) { // 40
            clarnv(IDIST, ISEED, JC-1, A( 1, JC ) );
         } // 40
      }

      // 4)      If SIM='T', apply similarity transformation.

                                 // -1
              // Transform is  X A X  , where X = U S V, thus

              // it is  U S V A V' (1/S) U'

      if ( ISIM != 0 ) {

         // Compute S (singular values of the eigenvector matrix)
         // according to CONDS and MODES

         slatm1(MODES, CONDS, 0, 0, ISEED, DS, N, IINFO );
         if ( IINFO != 0 ) {
            INFO = 3;
            return;
         }

         // Multiply by V and V'

         clarge(N, A, LDA, ISEED, WORK, IINFO );
         if ( IINFO != 0 ) {
            INFO = 4;
            return;
         }

         // Multiply by S and (1/S)

         for (J = 1; J <= N; J++) { // 50
            csscal(N, DS( J ), A( J, 1 ), LDA );
            if ( DS( J ) != ZERO ) {
               csscal(N, ONE / DS( J ), A( 1, J ), 1 );
            } else {
               INFO = 5;
               return;
            }
         } // 50

         // Multiply by U and U'

         clarge(N, A, LDA, ISEED, WORK, IINFO );
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

            ccopy(IROWS, A( JCR, IC ), 1, WORK, 1 );
            XNORMS = WORK( 1 );
            clarfg(IROWS, XNORMS, WORK( 2 ), 1, TAU );
            TAU = CONJG( TAU );
            WORK[1] = CONE;
            ALPHA = CLARND( 5, ISEED );

            cgemv('C', IROWS, ICOLS, CONE, A( JCR, IC+1 ), LDA, WORK, 1, CZERO, WORK( IROWS+1 ), 1 );
            cgerc(IROWS, ICOLS, -TAU, WORK, 1, WORK( IROWS+1 ), 1, A( JCR, IC+1 ), LDA );

            cgemv('N', N, IROWS, CONE, A( 1, JCR ), LDA, WORK, 1, CZERO, WORK( IROWS+1 ), 1 );
            cgerc(N, IROWS, -CONJG( TAU ), WORK( IROWS+1 ), 1, WORK, 1, A( 1, JCR ), LDA );

            A[JCR][IC] = XNORMS;
            claset('Full', IROWS-1, 1, CZERO, CZERO, A( JCR+1, IC ), LDA );

            cscal(ICOLS+1, ALPHA, A( JCR, IC ), LDA );
            cscal(N, CONJG( ALPHA ), A( 1, JCR ), 1 );
         } // 60
      } else if ( KU < N-1 ) {

         // Reduce upper bandwidth -- kill a row at a time.

         for (JCR = KU + 1; JCR <= N - 1; JCR++) { // 70
            IR = JCR - KU;
            IROWS = N + KU - JCR;
            ICOLS = N + 1 - JCR;

            ccopy(ICOLS, A( IR, JCR ), LDA, WORK, 1 );
            XNORMS = WORK( 1 );
            clarfg(ICOLS, XNORMS, WORK( 2 ), 1, TAU );
            TAU = CONJG( TAU );
            WORK[1] = CONE;
            clacgv(ICOLS-1, WORK( 2 ), 1 );
            ALPHA = CLARND( 5, ISEED );

            cgemv('N', IROWS, ICOLS, CONE, A( IR+1, JCR ), LDA, WORK, 1, CZERO, WORK( ICOLS+1 ), 1 );
            cgerc(IROWS, ICOLS, -TAU, WORK( ICOLS+1 ), 1, WORK, 1, A( IR+1, JCR ), LDA );

            cgemv('C', ICOLS, N, CONE, A( JCR, 1 ), LDA, WORK, 1, CZERO, WORK( ICOLS+1 ), 1 );
            cgerc(ICOLS, N, -CONJG( TAU ), WORK, 1, WORK( ICOLS+1 ), 1, A( JCR, 1 ), LDA );

            A[IR][JCR] = XNORMS;
            claset('Full', 1, ICOLS-1, CZERO, CZERO, A( IR, JCR+1 ), LDA );

            cscal(IROWS+1, ALPHA, A( IR, JCR ), 1 );
            cscal(N, CONJG( ALPHA ), A( JCR, 1 ), LDA );
         } // 70
      }

      // Scale the matrix to have norm ANORM

      if ( ANORM >= ZERO ) {
         TEMP = CLANGE( 'M', N, N, A, LDA, TEMPA );
         if ( TEMP > ZERO ) {
            RALPHA = ANORM / TEMP;
            for (J = 1; J <= N; J++) { // 80
               csscal(N, RALPHA, A( 1, J ), 1 );
            } // 80
         }
      }

      return;
      }
