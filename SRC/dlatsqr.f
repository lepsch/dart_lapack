      SUBROUTINE DLATSQR( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N, MB, NB, LDT, LWORK;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), WORK( * ), T( LDT, * );
      // ..

*  =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, II, KK, CTR, MINMN, LWMIN;
      // ..
      // .. EXTERNAL FUNCTIONS ..
      bool               LSAME;
      // EXTERNAL LSAME
      // .. EXTERNAL SUBROUTINES ..
      // EXTERNAL DGEQRT, DTPQRT, XERBLA
      // .. INTRINSIC FUNCTIONS ..
      // INTRINSIC MAX, MIN, MOD
      // ..
      // .. EXECUTABLE STATEMENTS ..

      // TEST THE INPUT ARGUMENTS

      INFO = 0

      LQUERY = ( LWORK == -1 )

      MINMN = MIN( M, N )
      if ( MINMN == 0 ) {
        LWMIN = 1
      } else {
        LWMIN = N*NB
      }

      if ( M.LT.0 ) {
        INFO = -1
      } else if ( N.LT.0 .OR. M.LT.N ) {
        INFO = -2
      } else if ( MB.LT.1 ) {
        INFO = -3
      } else if ( NB.LT.1 .OR. ( NB.GT.N .AND. N.GT.0 ) ) {
        INFO = -4
      } else if ( LDA.LT.MAX( 1, M ) ) {
        INFO = -6
      } else if ( LDT.LT.NB ) {
        INFO = -8
      } else if ( LWORK.LT.LWMIN .AND. (.NOT.LQUERY) ) {
        INFO = -10
      }

      if ( INFO == 0 ) {
        WORK( 1 ) = LWMIN
      }

      if ( INFO != 0 ) {
        xerbla('DLATSQR', -INFO );
        RETURN
      } else if ( LQUERY ) {
        RETURN
      }

      // Quick return if possible

      if ( MINMN == 0 ) {
        RETURN
      }

      // The QR Decomposition

      if ( (MB.LE.N) .OR. (MB.GE.M) ) {
        dgeqrt(M, N, NB, A, LDA, T, LDT, WORK, INFO );
        RETURN
      }

      KK = MOD((M-N),(MB-N))
      II = M-KK+1

      // Compute the QR factorization of the first block A(1:MB,1:N)

      dgeqrt(MB, N, NB, A(1,1), LDA, T, LDT, WORK, INFO );

      CTR = 1
      DO I = MB+1, II-MB+N, (MB-N)

        // Compute the QR factorization of the current block A(I:I+MB-N,1:N)

        dtpqrt(MB-N, N, 0, NB, A(1,1), LDA, A( I, 1 ), LDA, T(1, CTR * N + 1), LDT, WORK, INFO );
        CTR = CTR + 1
      }

      // Compute the QR factorization of the last block A(II:M,1:N)

      if ( II.LE.M ) {
        dtpqrt(KK, N, 0, NB, A(1,1), LDA, A( II, 1 ), LDA, T(1, CTR * N + 1), LDT, WORK, INFO );
      }

      WORK( 1 ) = LWMIN
      RETURN

      // End of DLATSQR

      }
