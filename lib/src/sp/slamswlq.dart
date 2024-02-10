      void slamswlq(final int SIDE, final int TRANS, final int M, final int N, final int K, final int MB, final int NB, final Matrix<double> A, final int LDA, final Matrix<double> T, final int LDT, final Matrix<double> C, final int LDC, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             SIDE, TRANS;
      int                INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC;
      double               A( LDA, * ), WORK( * ), C( LDC, * ), T( LDT, * );
      // ..

// =====================================================================

      bool               LEFT, RIGHT, TRAN, NOTRAN, LQUERY;
      int                I, II, KK, LW, CTR, MINMNK, LWMIN;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      double               SROUNDUP_LWORK;
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL STPMLQT, SGEMLQT, XERBLA

      // Test the input arguments

      LQUERY  = ( LWORK == -1 );
      NOTRAN  = lsame( TRANS, 'N' );
      TRAN    = lsame( TRANS, 'T' );
      LEFT    = lsame( SIDE, 'L' );
      RIGHT   = lsame( SIDE, 'R' );
      if ( LEFT ) {
        LW = N * MB;
      } else {
        LW = M * MB;
      }

      MINMNK = min( M, N, K );
      if ( MINMNK == 0 ) {
        LWMIN = 1;
      } else {
        LWMIN = max( 1, LW );
      }

      INFO = 0;
      if ( !LEFT && !RIGHT ) {
        INFO = -1;
      } else if ( !TRAN && !NOTRAN ) {
        INFO = -2;
      } else if ( K < 0 ) {
        INFO = -5;
      } else if ( M < K ) {
        INFO = -3;
      } else if ( N < 0 ) {
        INFO = -4;
      } else if ( K < MB || MB < 1 ) {
        INFO = -6;
      } else if ( LDA < max( 1, K ) ) {
        INFO = -9;
      } else if ( LDT < max( 1, MB ) ) {
        INFO = -11;
      } else if ( LDC < max( 1, M ) ) {
        INFO = -13;
      } else if ( LWORK < LWMIN && ( !LQUERY) ) {
        INFO = -15;
      }

      if ( INFO == 0 ) {
        WORK[1] = SROUNDUP_LWORK( LWMIN );
      }
      if ( INFO != 0 ) {
        xerbla('SLAMSWLQ', -INFO );
        return;
      } else if ( LQUERY ) {
        return;
      }

      // Quick return if possible

      if ( MINMNK == 0 ) {
        return;
      }

      if ((NB <= K) || (NB >= max(M,N,K))) {
        sgemlqt(SIDE, TRANS, M, N, K, MB, A, LDA, T, LDT, C, LDC, WORK, INFO);
        return;
      }

      if (LEFT && TRAN) {

          // Multiply Q to the last block of C

          KK = (M-K % NB-K);
          CTR = (M-K)/(NB-K);

          if (KK > 0) {
            II=M-KK+1;
            stpmlqt('L','T',KK , N, K, 0, MB, A(1,II), LDA, T(1,CTR*K+1), LDT, C(1,1), LDC, C(II,1), LDC, WORK, INFO );
          } else {
            II=M+1;
          }

          DO I=II-(NB-K),NB+1,-(NB-K);

          // Multiply Q to the current block of C (1:M,I:I+NB)

            CTR = CTR - 1;
            stpmlqt('L','T',NB-K , N, K, 0,MB, A(1,I), LDA, T(1,CTR*K+1),LDT, C(1,1), LDC, C(I,1), LDC, WORK, INFO );
          }

          // Multiply Q to the first block of C (1:M,1:NB)

          sgemlqt('L','T',NB , N, K, MB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO );

      } else if (LEFT && NOTRAN) {

          // Multiply Q to the first block of C

         KK = (M-K % NB-K);
         II=M-KK+1;
         CTR = 1;
         sgemlqt('L','N',NB , N, K, MB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO );

         for (I = NB+1; (NB-K) < 0 ? I >= II-NB+K : I <= II-NB+K; I += (NB-K)) {

          // Multiply Q to the current block of C (I:I+NB,1:N)

          stpmlqt('L','N',NB-K , N, K, 0,MB, A(1,I), LDA, T(1,CTR * K+1), LDT, C(1,1), LDC, C(I,1), LDC, WORK, INFO );
          CTR = CTR + 1;

         }
         if (II <= M) {

          // Multiply Q to the last block of C

          stpmlqt('L','N',KK , N, K, 0, MB, A(1,II), LDA, T(1,CTR*K+1), LDT, C(1,1), LDC, C(II,1), LDC, WORK, INFO );

         }

      } else if (RIGHT && NOTRAN) {

          // Multiply Q to the last block of C

          KK = (N-K % NB-K);
          CTR = (N-K)/(NB-K);
          if (KK > 0) {
            II=N-KK+1;
            stpmlqt('R','N',M , KK, K, 0, MB, A(1, II), LDA, T(1,CTR*K+1), LDT, C(1,1), LDC, C(1,II), LDC, WORK, INFO );
          } else {
            II=N+1;
          }

          DO I=II-(NB-K),NB+1,-(NB-K);

          // Multiply Q to the current block of C (1:M,I:I+MB)

             CTR = CTR - 1;
             stpmlqt('R','N', M, NB-K, K, 0, MB, A(1, I), LDA, T(1,CTR*K+1), LDT, C(1,1), LDC, C(1,I), LDC, WORK, INFO );

          }

          // Multiply Q to the first block of C (1:M,1:MB)

          sgemlqt('R','N',M , NB, K, MB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO );

      } else if (RIGHT && TRAN) {

        // Multiply Q to the first block of C

         KK = (N-K % NB-K);
         II=N-KK+1;
         CTR = 1;
         sgemlqt('R','T',M , NB, K, MB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO );

         for (I = NB+1; (NB-K) < 0 ? I >= II-NB+K : I <= II-NB+K; I += (NB-K)) {

          // Multiply Q to the current block of C (1:M,I:I+MB)

          stpmlqt('R','T',M , NB-K, K, 0,MB, A(1,I), LDA, T(1, CTR*K+1), LDT, C(1,1), LDC, C(1,I), LDC, WORK, INFO );
          CTR = CTR + 1;

         }
         if (II <= N) {

        // Multiply Q to the last block of C

          stpmlqt('R','T',M , KK, K, 0,MB, A(1,II), LDA, T(1,CTR*K+1),LDT, C(1,1), LDC, C(1,II), LDC, WORK, INFO );

         }

      }

      WORK[1] = SROUNDUP_LWORK( LWMIN );
      }
