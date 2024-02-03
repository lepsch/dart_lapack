      SUBROUTINE  DSB2ST_KERNELS( UPLO, WANTZ, TTYPE, ST, ED, SWEEP, N, NB, IB, A, LDA, V, TAU, LDVT, WORK)

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
      double             A( LDA, * ), V( * ), TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I, J1, J2, LM, LN, VPOS, TAUPOS, DPOS, OFDPOS, AJETER;
      double             CTMP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARFG, DLARFX, DLARFY
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

      IF( UPPER ) THEN
          DPOS    = 2 * NB + 1
          OFDPOS  = 2 * NB
      ELSE
          DPOS    = 1
          OFDPOS  = 2
      ENDIF


      // Upper case

      IF( UPPER ) THEN

          IF( WANTZ ) THEN
              VPOS   = MOD( SWEEP-1, 2 ) * N + ST
              TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
          ELSE
              VPOS   = MOD( SWEEP-1, 2 ) * N + ST
              TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
          ENDIF

          IF( TTYPE.EQ.1 ) THEN
              LM = ED - ST + 1

              V( VPOS ) = ONE
              DO 10 I = 1, LM-1
                  V( VPOS+I )         = ( A( OFDPOS-I, ST+I ) )
                  A( OFDPOS-I, ST+I ) = ZERO
   10         CONTINUE
              CTMP = ( A( OFDPOS, ST ) )
              CALL DLARFG( LM, CTMP, V( VPOS+1 ), 1, TAU( TAUPOS ) )
              A( OFDPOS, ST ) = CTMP

              LM = ED - ST + 1
              CALL DLARFY( UPLO, LM, V( VPOS ), 1, ( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK)
          ENDIF

          IF( TTYPE.EQ.3 ) THEN

              LM = ED - ST + 1
              CALL DLARFY( UPLO, LM, V( VPOS ), 1, ( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK)
          ENDIF

          IF( TTYPE.EQ.2 ) THEN
              J1 = ED+1
              J2 = MIN( ED+NB, N )
              LN = ED-ST+1
              LM = J2-J1+1
              IF( LM.GT.0) THEN
                  CALL DLARFX( 'Left', LN, LM, V( VPOS ), ( TAU( TAUPOS ) ), A( DPOS-NB, J1 ), LDA-1, WORK)

                  IF( WANTZ ) THEN
                      VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                      TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
                  ELSE
                      VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                      TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
                  ENDIF

                  V( VPOS ) = ONE
                  DO 30 I = 1, LM-1
                      V( VPOS+I )          = ( A( DPOS-NB-I, J1+I ) )
                      A( DPOS-NB-I, J1+I ) = ZERO
   30             CONTINUE
                  CTMP = ( A( DPOS-NB, J1 ) )
                  CALL DLARFG( LM, CTMP, V( VPOS+1 ), 1, TAU( TAUPOS ) )
                  A( DPOS-NB, J1 ) = CTMP

                  CALL DLARFX( 'Right', LN-1, LM, V( VPOS ), TAU( TAUPOS ), A( DPOS-NB+1, J1 ), LDA-1, WORK)
              ENDIF
          ENDIF

      // Lower case

      ELSE

          IF( WANTZ ) THEN
              VPOS   = MOD( SWEEP-1, 2 ) * N + ST
              TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
          ELSE
              VPOS   = MOD( SWEEP-1, 2 ) * N + ST
              TAUPOS = MOD( SWEEP-1, 2 ) * N + ST
          ENDIF

          IF( TTYPE.EQ.1 ) THEN
              LM = ED - ST + 1

              V( VPOS ) = ONE
              DO 20 I = 1, LM-1
                  V( VPOS+I )         = A( OFDPOS+I, ST-1 )
                  A( OFDPOS+I, ST-1 ) = ZERO
   20         CONTINUE
              CALL DLARFG( LM, A( OFDPOS, ST-1 ), V( VPOS+1 ), 1, TAU( TAUPOS ) )

              LM = ED - ST + 1

              CALL DLARFY( UPLO, LM, V( VPOS ), 1, ( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK)

          ENDIF

          IF( TTYPE.EQ.3 ) THEN
              LM = ED - ST + 1

              CALL DLARFY( UPLO, LM, V( VPOS ), 1, ( TAU( TAUPOS ) ), A( DPOS, ST ), LDA-1, WORK)

          ENDIF

          IF( TTYPE.EQ.2 ) THEN
              J1 = ED+1
              J2 = MIN( ED+NB, N )
              LN = ED-ST+1
              LM = J2-J1+1

              IF( LM.GT.0) THEN
                  CALL DLARFX( 'Right', LM, LN, V( VPOS ), TAU( TAUPOS ), A( DPOS+NB, ST ), LDA-1, WORK)

                  IF( WANTZ ) THEN
                      VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                      TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
                  ELSE
                      VPOS   = MOD( SWEEP-1, 2 ) * N + J1
                      TAUPOS = MOD( SWEEP-1, 2 ) * N + J1
                  ENDIF

                  V( VPOS ) = ONE
                  DO 40 I = 1, LM-1
                      V( VPOS+I )        = A( DPOS+NB+I, ST )
                      A( DPOS+NB+I, ST ) = ZERO
   40             CONTINUE
                  CALL DLARFG( LM, A( DPOS+NB, ST ), V( VPOS+1 ), 1, TAU( TAUPOS ) )

                  CALL DLARFX( 'Left', LM, LN-1, V( VPOS ), ( TAU( TAUPOS ) ), A( DPOS+NB-1, ST+1 ), LDA-1, WORK)

              ENDIF
          ENDIF
      ENDIF

      RETURN

      // End of DSB2ST_KERNELS

      }
