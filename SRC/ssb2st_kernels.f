      SUBROUTINE  SSB2ST_KERNELS( UPLO, WANTZ, TTYPE, ST, ED, SWEEP, N, NB, IB, A, LDA, V, TAU, LDVT, WORK)

      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      bool               WANTZ;
      int                TTYPE, ST, ED, SWEEP, N, NB, IB, LDA, LDVT;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), V( * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, J1, J2, LM, LN, VPOS, TAUPOS, DPOS, OFDPOS, AJETER;
      REAL               CTMP
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARFG, SLARFX, SLARFY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MOD
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // ..
      // .. Executable Statements ..

      AJETER = IB + LDVT
      UPPER = LSAME( UPLO, 'U' )

      if ( UPPER ) {
          DPOS    = 2 * NB + 1
          OFDPOS  = 2 * NB
      } else {
          DPOS    = 1
          OFDPOS  = 2
      }


      // Upper case

      if ( UPPER ) {

          if ( WANTZ ) {
              VPOS   = MOD( SWEEP-1, 2 ) * N + ST
              TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
          } else {
              VPOS   = MOD( SWEEP-1, 2 ) * N + ST
              TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
          }

          if ( TTYPE == 1 ) {
              LM = ED - ST + 1

              V( VPOS ) = ONE
              for (I = 1; I <= LM-1; I++) { // 10
                  V( VPOS+I )         = ( A( OFDPOS-I, ST+I ) )
                  A( OFDPOS-I, ST+I ) = ZERO
              } // 10
              CTMP = ( A( OFDPOS, ST ) )
              slarfg(LM, CTMP, V( VPOS+1 ), 1, TAU( TAUPOS ) );
              A( OFDPOS, ST ) = CTMP

              LM = ED - ST + 1
              slarfy(UPLO, LM, V( VPOS ), 1, ( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK);
          }

          if ( TTYPE == 3 ) {

              LM = ED - ST + 1
              slarfy(UPLO, LM, V( VPOS ), 1, ( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK);
          }

          if ( TTYPE == 2 ) {
              J1 = ED+1
              J2 = MIN( ED+NB, N )
              LN = ED-ST+1
              LM = J2-J1+1
              if ( LM.GT.0) {
                  slarfx('Left', LN, LM, V( VPOS ), ( TAU( TAUPOS ) ), A( DPOS-NB, J1 ), LDA-1, WORK);

                  if ( WANTZ ) {
                      VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                      TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
                  } else {
                      VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                      TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
                  }

                  V( VPOS ) = ONE
                  for (I = 1; I <= LM-1; I++) { // 30
                      V( VPOS+I )          = ( A( DPOS-NB-I, J1+I ) )
                      A( DPOS-NB-I, J1+I ) = ZERO
                  } // 30
                  CTMP = ( A( DPOS-NB, J1 ) )
                  slarfg(LM, CTMP, V( VPOS+1 ), 1, TAU( TAUPOS ) );
                  A( DPOS-NB, J1 ) = CTMP

                  slarfx('Right', LN-1, LM, V( VPOS ), TAU( TAUPOS ), A( DPOS-NB+1, J1 ), LDA-1, WORK);
              }
          }

      // Lower case

      } else {

          if ( WANTZ ) {
              VPOS   = MOD( SWEEP-1, 2 ) * N + ST
              TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
          } else {
              VPOS   = MOD( SWEEP-1, 2 ) * N + ST
              TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
          }

          if ( TTYPE == 1 ) {
              LM = ED - ST + 1

              V( VPOS ) = ONE
              for (I = 1; I <= LM-1; I++) { // 20
                  V( VPOS+I )         = A( OFDPOS+I, ST-1 )
                  A( OFDPOS+I, ST-1 ) = ZERO
              } // 20
              slarfg(LM, A( OFDPOS, ST-1 ), V( VPOS+1 ), 1, TAU( TAUPOS ) );

              LM = ED - ST + 1

              slarfy(UPLO, LM, V( VPOS ), 1, ( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK);

          }

          if ( TTYPE == 3 ) {
              LM = ED - ST + 1

              slarfy(UPLO, LM, V( VPOS ), 1, ( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK);

          }

          if ( TTYPE == 2 ) {
              J1 = ED+1
              J2 = MIN( ED+NB, N )
              LN = ED-ST+1
              LM = J2-J1+1

              if ( LM.GT.0) {
                  slarfx('Right', LM, LN, V( VPOS ), TAU( TAUPOS ), A( DPOS+NB, ST ), LDA-1, WORK);

                  if ( WANTZ ) {
                      VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                      TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
                  } else {
                      VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                      TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
                  }

                  V( VPOS ) = ONE
                  for (I = 1; I <= LM-1; I++) { // 40
                      V( VPOS+I )        = A( DPOS+NB+I, ST )
                      A( DPOS+NB+I, ST ) = ZERO
                  } // 40
                  slarfg(LM, A( DPOS+NB, ST ), V( VPOS+1 ), 1, TAU( TAUPOS ) );

                  slarfx('Left', LM, LN-1, V( VPOS ), ( TAU( TAUPOS ) ), A( DPOS+NB-1, ST+1 ), LDA-1, WORK);

              }
          }
      }

      RETURN

      // End of SSB2ST_KERNELS

      }
