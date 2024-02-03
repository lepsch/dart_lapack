      SUBROUTINE  CHB2ST_KERNELS( UPLO, WANTZ, TTYPE, ST, ED, SWEEP, N, NB, IB, A, LDA, V, TAU, LDVT, WORK)

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
      COMPLEX            A( LDA, * ), V( * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0E+0, 0.0E+0 ), ONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, J1, J2, LM, LN, VPOS, TAUPOS, DPOS, OFDPOS, AJETER;
      COMPLEX            CTMP
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARFG, CLARFX, CLARFY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MOD
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
      ENDIF


      // Upper case

      if ( UPPER ) {

          if ( WANTZ ) {
              VPOS   = MOD( SWEEP-1, 2 ) * N + ST
              TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
          } else {
              VPOS   = MOD( SWEEP-1, 2 ) * N + ST
              TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
          ENDIF

          if ( TTYPE.EQ.1 ) {
              LM = ED - ST + 1

              V( VPOS ) = ONE
              DO 10 I = 1, LM-1
                  V( VPOS+I )         = CONJG( A( OFDPOS-I, ST+I ) )
                  A( OFDPOS-I, ST+I ) = ZERO
              } // 10
              CTMP = CONJG( A( OFDPOS, ST ) )
              clarfg(LM, CTMP, V( VPOS+1 ), 1, TAU( TAUPOS ) );
              A( OFDPOS, ST ) = CTMP

              LM = ED - ST + 1
              clarfy(UPLO, LM, V( VPOS ), 1, CONJG( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK);
          ENDIF

          if ( TTYPE.EQ.3 ) {

              LM = ED - ST + 1
              clarfy(UPLO, LM, V( VPOS ), 1, CONJG( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK);
          ENDIF

          if ( TTYPE.EQ.2 ) {
              J1 = ED+1
              J2 = MIN( ED+NB, N )
              LN = ED-ST+1
              LM = J2-J1+1
              if ( LM.GT.0) {
                  clarfx('Left', LN, LM, V( VPOS ), CONJG( TAU( TAUPOS ) ), A( DPOS-NB, J1 ), LDA-1, WORK);

                  if ( WANTZ ) {
                      VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                      TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
                  } else {
                      VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                      TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
                  ENDIF

                  V( VPOS ) = ONE
                  DO 30 I = 1, LM-1
                      V( VPOS+I )          = CONJG( A( DPOS-NB-I, J1+I ) )
                      A( DPOS-NB-I, J1+I ) = ZERO
                  } // 30
                  CTMP = CONJG( A( DPOS-NB, J1 ) )
                  clarfg(LM, CTMP, V( VPOS+1 ), 1, TAU( TAUPOS ) );
                  A( DPOS-NB, J1 ) = CTMP

                  clarfx('Right', LN-1, LM, V( VPOS ), TAU( TAUPOS ), A( DPOS-NB+1, J1 ), LDA-1, WORK);
              ENDIF
          ENDIF

      // Lower case

      } else {

          if ( WANTZ ) {
              VPOS   = MOD( SWEEP-1, 2 ) * N + ST
              TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
          } else {
              VPOS   = MOD( SWEEP-1, 2 ) * N + ST
              TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
          ENDIF

          if ( TTYPE.EQ.1 ) {
              LM = ED - ST + 1

              V( VPOS ) = ONE
              DO 20 I = 1, LM-1
                  V( VPOS+I )         = A( OFDPOS+I, ST-1 )
                  A( OFDPOS+I, ST-1 ) = ZERO
              } // 20
              clarfg(LM, A( OFDPOS, ST-1 ), V( VPOS+1 ), 1, TAU( TAUPOS ) );

              LM = ED - ST + 1

              clarfy(UPLO, LM, V( VPOS ), 1, CONJG( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK);

          ENDIF

          if ( TTYPE.EQ.3 ) {
              LM = ED - ST + 1

              clarfy(UPLO, LM, V( VPOS ), 1, CONJG( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK);

          ENDIF

          if ( TTYPE.EQ.2 ) {
              J1 = ED+1
              J2 = MIN( ED+NB, N )
              LN = ED-ST+1
              LM = J2-J1+1

              if ( LM.GT.0) {
                  clarfx('Right', LM, LN, V( VPOS ), TAU( TAUPOS ), A( DPOS+NB, ST ), LDA-1, WORK);

                  if ( WANTZ ) {
                      VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                      TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
                  } else {
                      VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                      TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
                  ENDIF

                  V( VPOS ) = ONE
                  DO 40 I = 1, LM-1
                      V( VPOS+I )        = A( DPOS+NB+I, ST )
                      A( DPOS+NB+I, ST ) = ZERO
                  } // 40
                  clarfg(LM, A( DPOS+NB, ST ), V( VPOS+1 ), 1, TAU( TAUPOS ) );

                  clarfx('Left', LM, LN-1, V( VPOS ), CONJG( TAU( TAUPOS ) ), A( DPOS+NB-1, ST+1 ), LDA-1, WORK);

              ENDIF
          ENDIF
      ENDIF

      RETURN

      // End of CHB2ST_KERNELS

      }
