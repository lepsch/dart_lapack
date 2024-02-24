      void slatsqr(final int M, final int N, final int MB, final int NB, final Matrix<double> A_, final int LDA, final Matrix<double> T_, final int LDT, final Array<double> WORK_, final int LWORK, final Box<int> INFO,) {
  final A = A_.dim();
  final T = T_.dim();
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
      int                INFO, LDA, M, N, MB, NB, LDT, LWORK;
      double               A( LDA, * ), WORK( * ), T( LDT, * );
      // ..

// =====================================================================

      bool               LQUERY;
      int                I, II, KK, CTR, MINMN, LWMIN;
      // ..
      // .. EXTERNAL FUNCTIONS ..
      //- bool               lsame;
      // EXTERNAL lsame
      double               SROUNDUP_LWORK;
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. EXTERNAL SUBROUTINES ..
      // EXTERNAL SGEQRT, STPQRT, XERBLA
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
        LWMIN = N*NB;
      }

      if ( M < 0 ) {
        INFO = -1;
      } else if ( N < 0 || M < N ) {
        INFO = -2;
      } else if ( MB < 1 ) {
        INFO = -3;
      } else if ( NB < 1 || ( NB > N && N > 0 ) ) {
        INFO = -4;
      } else if ( LDA < max( 1, M ) ) {
        INFO = -6;
      } else if ( LDT < NB ) {
        INFO = -8;
      } else if ( LWORK < LWMIN && ( !LQUERY) ) {
        INFO = -10;
      }

      if ( INFO == 0 ) {
        WORK[1] = SROUNDUP_LWORK( LWMIN );
      }
      if ( INFO != 0 ) {
        xerbla('SLATSQR', -INFO );
        return;
      } else if ( LQUERY ) {
        return;
      }

      // Quick return if possible

      if ( MINMN == 0 ) {
        return;
      }

      // The QR Decomposition

      if ( (MB <= N) || (MB >= M) ) {
        sgeqrt(M, N, NB, A, LDA, T, LDT, WORK, INFO );
        return;
      }
      KK = ((M-N) % (MB-N));
      II = M-KK+1;

      // Compute the QR factorization of the first block A(1:MB,1:N)

      sgeqrt(MB, N, NB, A(1,1), LDA, T, LDT, WORK, INFO );

      CTR = 1;
      for (I = MB+1; (MB-N) < 0 ? I >= II-MB+N : I <= II-MB+N; I += (MB-N)) {

        // Compute the QR factorization of the current block A(I:I+MB-N,1:N)

        stpqrt(MB-N, N, 0, NB, A(1,1), LDA, A( I, 1 ), LDA, T(1, CTR * N + 1), LDT, WORK, INFO );
        CTR = CTR + 1;
      }

      // Compute the QR factorization of the last block A(II:M,1:N)

      if ( II <= M ) {
        stpqrt(KK, N, 0, NB, A(1,1), LDA, A( II, 1 ), LDA, T(1, CTR * N + 1), LDT, WORK, INFO );
      }

      WORK[1] = SROUNDUP_LWORK( LWMIN );
      }
