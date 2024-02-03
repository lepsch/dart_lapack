      SUBROUTINE ZLATSQR( M, N, MB, NB, A, LDA, T, LDT, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N, MB, NB, LDT, LWORK;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), WORK( * ), T( LDT, * )
      // ..

*  =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, II, KK, CTR, LWMIN, MINMN;
      // ..
      // .. EXTERNAL FUNCTIONS ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. EXTERNAL SUBROUTINES ..
      // EXTERNAL ZGEQRT, ZTPQRT, XERBLA
      // ..
      // .. INTRINSIC FUNCTIONS ..
      // INTRINSIC MAX, MIN, MOD
      // ..
      // .. EXECUTABLE STATEMENTS ..

      // TEST THE INPUT ARGUMENTS

      INFO = 0

      LQUERY = ( LWORK.EQ.-1 )

      MINMN = MIN( M, N )
      if ( MINMN.EQ.0 ) {
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

      if ( INFO.EQ.0 ) {
        WORK( 1 ) = LWMIN
      }
      if ( INFO.NE.0 ) {
        CALL XERBLA( 'ZLATSQR', -INFO )
        RETURN
      } else if ( LQUERY ) {
        RETURN
      }

      // Quick return if possible

      if ( MINMN.EQ.0 ) {
        RETURN
      }

      // The QR Decomposition

      if ( (MB.LE.N) .OR. (MB.GE.M) ) {
        CALL ZGEQRT( M, N, NB, A, LDA, T, LDT, WORK, INFO )
        RETURN
      }
      KK = MOD((M-N),(MB-N))
      II = M-KK+1

      // Compute the QR factorization of the first block A(1:MB,1:N)

      CALL ZGEQRT( MB, N, NB, A(1,1), LDA, T, LDT, WORK, INFO )
      CTR = 1

      DO I = MB+1, II-MB+N, (MB-N)

        // Compute the QR factorization of the current block A(I:I+MB-N,1:N)

        CALL ZTPQRT( MB-N, N, 0, NB, A(1,1), LDA, A( I, 1 ), LDA, T(1, CTR * N + 1), LDT, WORK, INFO )
        CTR = CTR + 1
      END DO

      // Compute the QR factorization of the last block A(II:M,1:N)

      if ( II.LE.M ) {
        CALL ZTPQRT( KK, N, 0, NB, A(1,1), LDA, A( II, 1 ), LDA, T(1,CTR * N + 1), LDT, WORK, INFO )
      }

      WORK( 1 ) = LWMIN
      RETURN

      // End of ZLATSQR

      }
