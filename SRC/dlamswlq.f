      SUBROUTINE DLAMSWLQ( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T, LDT, C, LDC, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), WORK( * ), C( LDC, * ), T( LDT, * );
      // ..

* =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LEFT, RIGHT, TRAN, NOTRAN, LQUERY;
      int                I, II, KK, CTR, LW, MINMNK, LWMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // .. External Subroutines ..
      // EXTERNAL DTPMLQT, DGEMLQT, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      LQUERY  = ( LWORK.EQ.-1 )
      NOTRAN  = LSAME( TRANS, 'N' )
      TRAN    = LSAME( TRANS, 'T' )
      LEFT    = LSAME( SIDE, 'L' )
      RIGHT   = LSAME( SIDE, 'R' )
      if ( LEFT ) {
        LW = N * MB
      } else {
        LW = M * MB
      }

      MINMNK = MIN( M, N, K )
      if ( MINMNK.EQ.0 ) {
        LWMIN = 1
      } else {
        LWMIN = MAX( 1, LW )
      }

      INFO = 0
      if ( .NOT.LEFT .AND. .NOT.RIGHT ) {
        INFO = -1
      } else if ( .NOT.TRAN .AND. .NOT.NOTRAN ) {
        INFO = -2
      } else if ( K.LT.0 ) {
        INFO = -5
      } else if ( M.LT.K ) {
        INFO = -3
      } else if ( N.LT.0 ) {
        INFO = -4
      } else if ( K.LT.MB .OR. MB.LT.1 ) {
        INFO = -6
      } else if ( LDA.LT.MAX( 1, K ) ) {
        INFO = -9
      } else if ( LDT.LT.MAX( 1, MB ) ) {
        INFO = -11
      } else if ( LDC.LT.MAX( 1, M ) ) {
        INFO = -13
      } else if ( LWORK.LT.LWMIN .AND. (.NOT.LQUERY) ) {
        INFO = -15
      }

      if ( INFO.EQ.0 ) {
        WORK( 1 ) = LWMIN
      }
      if ( INFO.NE.0 ) {
        CALL XERBLA( 'DLAMSWLQ', -INFO )
        RETURN
      } else if ( LQUERY ) {
        RETURN
      }

      // Quick return if possible

      if ( MINMNK.EQ.0 ) {
        RETURN
      }

      if ((NB.LE.K).OR.(NB.GE.MAX(M,N,K))) {
        CALL DGEMLQT( SIDE, TRANS, M, N, K, MB, A, LDA, T, LDT, C, LDC, WORK, INFO)
        RETURN
      }

      if (LEFT.AND.TRAN) {

          // Multiply Q to the last block of C

          KK = MOD((M-K),(NB-K))
          CTR = (M-K)/(NB-K)
          if (KK.GT.0) {
            II=M-KK+1
            CALL DTPMLQT('L','T',KK , N, K, 0, MB, A(1,II), LDA, T(1,CTR*K+1), LDT, C(1,1), LDC, C(II,1), LDC, WORK, INFO )
          } else {
            II=M+1
          }

          DO I=II-(NB-K),NB+1,-(NB-K)

          // Multiply Q to the current block of C (1:M,I:I+NB)

            CTR = CTR - 1
            CALL DTPMLQT('L','T',NB-K , N, K, 0,MB, A(1,I), LDA, T(1, CTR*K+1),LDT, C(1,1), LDC, C(I,1), LDC, WORK, INFO )

          END DO

          // Multiply Q to the first block of C (1:M,1:NB)

          CALL DGEMLQT('L','T',NB , N, K, MB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO )

      } else if (LEFT.AND.NOTRAN) {

          // Multiply Q to the first block of C

         KK = MOD((M-K),(NB-K))
         II=M-KK+1
         CTR = 1
         CALL DGEMLQT('L','N',NB , N, K, MB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO )

         DO I=NB+1,II-NB+K,(NB-K)

          // Multiply Q to the current block of C (I:I+NB,1:N)

          CALL DTPMLQT('L','N',NB-K , N, K, 0,MB, A(1,I), LDA, T(1,CTR*K+1), LDT, C(1,1), LDC, C(I,1), LDC, WORK, INFO )
          CTR = CTR + 1

         END DO
         if (II.LE.M) {

          // Multiply Q to the last block of C

          CALL DTPMLQT('L','N',KK , N, K, 0, MB, A(1,II), LDA, T(1,CTR*K+1), LDT, C(1,1), LDC, C(II,1), LDC, WORK, INFO )

         }

      } else if (RIGHT.AND.NOTRAN) {

          // Multiply Q to the last block of C

          KK = MOD((N-K),(NB-K))
          CTR = (N-K)/(NB-K)
          if (KK.GT.0) {
            II=N-KK+1
            CALL DTPMLQT('R','N',M , KK, K, 0, MB, A(1, II), LDA, T(1,CTR *K+1), LDT, C(1,1), LDC, C(1,II), LDC, WORK, INFO )
          } else {
            II=N+1
          }

          DO I=II-(NB-K),NB+1,-(NB-K)

          // Multiply Q to the current block of C (1:M,I:I+MB)

             CTR = CTR - 1
             CALL DTPMLQT('R','N', M, NB-K, K, 0, MB, A(1, I), LDA, T(1,CTR*K+1), LDT, C(1,1), LDC, C(1,I), LDC, WORK, INFO )

          END DO

          // Multiply Q to the first block of C (1:M,1:MB)

          CALL DGEMLQT('R','N',M , NB, K, MB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO )

      } else if (RIGHT.AND.TRAN) {

        // Multiply Q to the first block of C

         KK = MOD((N-K),(NB-K))
         CTR = 1
         II=N-KK+1
         CALL DGEMLQT('R','T',M , NB, K, MB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO )

         DO I=NB+1,II-NB+K,(NB-K)

          // Multiply Q to the current block of C (1:M,I:I+MB)

          CALL DTPMLQT('R','T',M , NB-K, K, 0,MB, A(1,I), LDA, T(1,CTR*K+1), LDT, C(1,1), LDC, C(1,I), LDC, WORK, INFO )
          CTR = CTR + 1

         END DO
         if (II.LE.N) {

        // Multiply Q to the last block of C

          CALL DTPMLQT('R','T',M , KK, K, 0,MB, A(1,II), LDA, T(1,CTR*K+1),LDT, C(1,1), LDC, C(1,II), LDC, WORK, INFO )

         }

      }

      WORK( 1 ) = LWMIN

      RETURN

      // End of DLAMSWLQ

      }
