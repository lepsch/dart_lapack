      SUBROUTINE DLASWLQ( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N, MB, NB, LWORK, LDT;
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
      // ..
      // .. EXTERNAL SUBROUTINES ..
      // EXTERNAL DGELQT, DTPLQT, XERBLA
      // ..
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
         LWMIN = M*MB
      }

      if ( M.LT.0 ) {
        INFO = -1
      } else if ( N.LT.0 .OR. N.LT.M ) {
        INFO = -2
      } else if ( MB.LT.1 .OR. ( MB.GT.M .AND. M.GT.0 ) ) {
        INFO = -3
      } else if ( NB.LT.0 ) {
        INFO = -4
      } else if ( LDA.LT.MAX( 1, M ) ) {
        INFO = -6
      } else if ( LDT.LT.MB ) {
        INFO = -8
      } else if ( LWORK.LT.LWMIN .AND. (.NOT.LQUERY) ) {
        INFO = -10
      }

      if ( INFO == 0 ) {
        WORK( 1 ) = LWMIN
      }

      if ( INFO.NE.0 ) {
        xerbla('DLASWLQ', -INFO );
        RETURN
      } else if ( LQUERY ) {
        RETURN
      }

      // Quick return if possible

      if ( MINMN == 0 ) {
        RETURN
      }

      // The LQ Decomposition

      if ( (M.GE.N) .OR. (NB.LE.M) .OR. (NB.GE.N) ) {
        dgelqt(M, N, MB, A, LDA, T, LDT, WORK, INFO );
        RETURN
      }

      KK = MOD((N-M),(NB-M))
      II = N-KK+1

      // Compute the LQ factorization of the first block A(1:M,1:NB)

      dgelqt(M, NB, MB, A(1,1), LDA, T, LDT, WORK, INFO );
      CTR = 1

      DO I = NB+1, II-NB+M, (NB-M)

        // Compute the QR factorization of the current block A(1:M,I:I+NB-M)

        dtplqt(M, NB-M, 0, MB, A(1,1), LDA, A( 1, I ), LDA, T(1, CTR * M + 1), LDT, WORK, INFO );
        CTR = CTR + 1
      }

      // Compute the QR factorization of the last block A(1:M,II:N)

      if ( II.LE.N ) {
        dtplqt(M, KK, 0, MB, A(1,1), LDA, A( 1, II ), LDA, T(1, CTR * M + 1), LDT, WORK, INFO );
      }

      WORK( 1 ) = LWMIN

      RETURN

      // End of DLASWLQ

      }
