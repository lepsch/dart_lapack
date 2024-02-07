      void ssb2st_kernels(UPLO, WANTZ, TTYPE, ST, ED, SWEEP, N, NB, IB, A, LDA, V, TAU, LDVT, WORK) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      bool               WANTZ;
      int                TTYPE, ST, ED, SWEEP, N, NB, IB, LDA, LDVT;
      double               A( LDA, * ), V( * ), TAU( * ), WORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      bool               UPPER;
      int                I, J1, J2, LM, LN, VPOS, TAUPOS, DPOS, OFDPOS, AJETER;
      double               CTMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARFG, SLARFX, SLARFY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..

      AJETER = IB + LDVT;
      UPPER = lsame( UPLO, 'U' );

      if ( UPPER ) {
          DPOS    = 2 * NB + 1;
          OFDPOS  = 2 * NB;
      } else {
          DPOS    = 1;
          OFDPOS  = 2;
      }


      // Upper case

      if ( UPPER ) {

          if ( WANTZ ) {
              VPOS   = (SWEEP-1 % 2) * N + ST;
              TAUPOS = (SWEEP-1 % 2) * N + ST;
          } else {
              VPOS   = (SWEEP-1 % 2) * N + ST;
              TAUPOS = (SWEEP-1 % 2) * N + ST;
          }

          if ( TTYPE == 1 ) {
              LM = ED - ST + 1;

              V[VPOS] = ONE;
              for (I = 1; I <= LM-1; I++) { // 10
                  V[VPOS+I] = ( A( OFDPOS-I, ST+I ) );
                  A[OFDPOS-I][ST+I] = ZERO;
              } // 10
              CTMP = ( A( OFDPOS, ST ) );
              slarfg(LM, CTMP, V( VPOS+1 ), 1, TAU( TAUPOS ) );
              A[OFDPOS][ST] = CTMP;

              LM = ED - ST + 1;
              slarfy(UPLO, LM, V( VPOS ), 1, ( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK);
          }

          if ( TTYPE == 3 ) {

              LM = ED - ST + 1;
              slarfy(UPLO, LM, V( VPOS ), 1, ( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK);
          }

          if ( TTYPE == 2 ) {
              J1 = ED+1;
              J2 = min( ED+NB, N );
              LN = ED-ST+1;
              LM = J2-J1+1;
              if ( LM > 0) {
                  slarfx('Left', LN, LM, V( VPOS ), ( TAU( TAUPOS ) ), A( DPOS-NB, J1 ), LDA-1, WORK);

                  if ( WANTZ ) {
                      VPOS   = (SWEEP-1 % 2) * N + J1;
                      TAUPOS = (SWEEP-1 % 2) * N + J1;
                  } else {
                      VPOS   = (SWEEP-1 % 2) * N + J1;
                      TAUPOS = (SWEEP-1 % 2) * N + J1;
                  }

                  V[VPOS] = ONE;
                  for (I = 1; I <= LM-1; I++) { // 30
                      V[VPOS+I] = ( A( DPOS-NB-I, J1+I ) );
                      A[DPOS-NB-I][J1+I] = ZERO;
                  } // 30
                  CTMP = ( A( DPOS-NB, J1 ) );
                  slarfg(LM, CTMP, V( VPOS+1 ), 1, TAU( TAUPOS ) );
                  A[DPOS-NB][J1] = CTMP;

                  slarfx('Right', LN-1, LM, V( VPOS ), TAU( TAUPOS ), A( DPOS-NB+1, J1 ), LDA-1, WORK);
              }
          }

      // Lower case

      } else {

          if ( WANTZ ) {
              VPOS   = (SWEEP-1 % 2) * N + ST;
              TAUPOS = (SWEEP-1 % 2) * N + ST;
          } else {
              VPOS   = (SWEEP-1 % 2) * N + ST;
              TAUPOS = (SWEEP-1 % 2) * N + ST;
          }

          if ( TTYPE == 1 ) {
              LM = ED - ST + 1;

              V[VPOS] = ONE;
              for (I = 1; I <= LM-1; I++) { // 20
                  V[VPOS+I] = A( OFDPOS+I, ST-1 );
                  A[OFDPOS+I][ST-1] = ZERO;
              } // 20
              slarfg(LM, A( OFDPOS, ST-1 ), V( VPOS+1 ), 1, TAU( TAUPOS ) );

              LM = ED - ST + 1;

              slarfy(UPLO, LM, V( VPOS ), 1, ( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK);

          }

          if ( TTYPE == 3 ) {
              LM = ED - ST + 1;

              slarfy(UPLO, LM, V( VPOS ), 1, ( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK);

          }

          if ( TTYPE == 2 ) {
              J1 = ED+1;
              J2 = min( ED+NB, N );
              LN = ED-ST+1;
              LM = J2-J1+1;

              if ( LM > 0) {
                  slarfx('Right', LM, LN, V( VPOS ), TAU( TAUPOS ), A( DPOS+NB, ST ), LDA-1, WORK);

                  if ( WANTZ ) {
                      VPOS   = (SWEEP-1 % 2) * N + J1;
                      TAUPOS = (SWEEP-1 % 2) * N + J1;
                  } else {
                      VPOS   = (SWEEP-1 % 2) * N + J1;
                      TAUPOS = (SWEEP-1 % 2) * N + J1;
                  }

                  V[VPOS] = ONE;
                  for (I = 1; I <= LM-1; I++) { // 40
                      V[VPOS+I] = A( DPOS+NB+I, ST );
                      A[DPOS+NB+I][ST] = ZERO;
                  } // 40
                  slarfg(LM, A( DPOS+NB, ST ), V( VPOS+1 ), 1, TAU( TAUPOS ) );

                  slarfx('Left', LM, LN-1, V( VPOS ), ( TAU( TAUPOS ) ), A( DPOS+NB-1, ST+1 ), LDA-1, WORK);

              }
          }
      }

      return;
      }
