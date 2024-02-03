      void zlaswlq(M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N, MB, NB, LWORK, LDT;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), WORK( * ), T( LDT, * );
      // ..

// =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, II, KK, CTR, MINMN, LWMIN;
      // ..
      // .. EXTERNAL FUNCTIONS ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. EXTERNAL SUBROUTINES ..
      // EXTERNAL ZGELQT, ZTPLQT, XERBLA
      // ..
      // .. INTRINSIC FUNCTIONS ..
      // INTRINSIC MAX, MIN, MOD
      // ..
      // .. EXECUTABLE STATEMENTS ..

      // TEST THE INPUT ARGUMENTS

      INFO = 0;

      LQUERY = ( LWORK == -1 );

      MINMN = min( M, N );
      if ( MINMN == 0 ) {
        LWMIN = 1;
      } else {
        LWMIN = M*MB;
      }

      if ( M < 0 ) {
        INFO = -1;
      } else if ( N < 0 || N < M ) {
        INFO = -2;
      } else if ( MB < 1 || ( MB > M && M > 0 ) ) {
        INFO = -3;
      } else if ( NB <= 0 ) {
        INFO = -4;
      } else if ( LDA < max( 1, M ) ) {
        INFO = -6;
      } else if ( LDT < MB ) {
        INFO = -8;
      } else if ( LWORK < LWMIN && ( !LQUERY) ) {
        INFO = -10;
      }

      if ( INFO == 0 ) {
        WORK( 1 ) = LWMIN;
      }

      if ( INFO != 0 ) {
        xerbla('ZLASWLQ', -INFO );
        return;
      } else if ( LQUERY ) {
        return;
      }

      // Quick return if possible

      if ( MINMN == 0 ) {
        return;
      }

      // The LQ Decomposition

      if ( (M >= N) || (NB <= M) || (NB >= N) ) {
        zgelqt(M, N, MB, A, LDA, T, LDT, WORK, INFO );
        return;
      }

      KK = MOD((N-M),(NB-M));
      II = N-KK+1;

      // Compute the LQ factorization of the first block A(1:M,1:NB)

      zgelqt(M, NB, MB, A(1,1), LDA, T, LDT, WORK, INFO );
      CTR = 1;

      DO I = NB+1, II-NB+M, (NB-M);

        // Compute the QR factorization of the current block A(1:M,I:I+NB-M)

        ztplqt(M, NB-M, 0, MB, A(1,1), LDA, A( 1, I ), LDA, T(1, CTR * M + 1), LDT, WORK, INFO );
        CTR = CTR + 1;
      }

      // Compute the QR factorization of the last block A(1:M,II:N)

      if ( II <= N ) {
        ztplqt(M, KK, 0, MB, A(1,1), LDA, A( 1, II ), LDA, T(1, CTR * M + 1), LDT, WORK, INFO );
      }

      WORK( 1 ) = LWMIN;
      return;

      // End of ZLASWLQ

      }
