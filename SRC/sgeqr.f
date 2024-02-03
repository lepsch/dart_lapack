      SUBROUTINE SGEQR( M, N, A, LDA, T, TSIZE, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N, TSIZE, LWORK;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), T( * ), WORK( * )
      // ..

*  =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LQUERY, LMINWS, MINT, MINW;
      int                MB, NB, MINTSZ, NBLCKS, LWMIN, LWREQ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      REAL               SROUNDUP_LWORK
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLATSQR, SGEQRT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, MOD
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Executable statements ..

      // Test the input arguments

      INFO = 0

      LQUERY = ( TSIZE.EQ.-1 .OR. TSIZE.EQ.-2 .OR. LWORK.EQ.-1 .OR. LWORK.EQ.-2 )

      MINT = .FALSE.
      MINW = .FALSE.
      IF( TSIZE.EQ.-2 .OR. LWORK.EQ.-2 ) THEN
        IF( TSIZE.NE.-1 ) MINT = .TRUE.
        IF( LWORK.NE.-1 ) MINW = .TRUE.
      END IF

      // Determine the block size

      IF( MIN( M, N ).GT.0 ) THEN
        MB = ILAENV( 1, 'SGEQR ', ' ', M, N, 1, -1 )
        NB = ILAENV( 1, 'SGEQR ', ' ', M, N, 2, -1 )
      ELSE
        MB = M
        NB = 1
      END IF
      IF( MB.GT.M .OR. MB.LE.N ) MB = M
      IF( NB.GT.MIN( M, N ) .OR. NB.LT.1 ) NB = 1
      MINTSZ = N + 5
      IF ( MB.GT.N .AND. M.GT.N ) THEN
        IF( MOD( M - N, MB - N ).EQ.0 ) THEN
          NBLCKS = ( M - N ) / ( MB - N )
        ELSE
          NBLCKS = ( M - N ) / ( MB - N ) + 1
        END IF
      ELSE
        NBLCKS = 1
      END IF

      // Determine if the workspace size satisfies minimal size

      LWMIN = MAX( 1, N )
      LWREQ = MAX( 1, N*NB )
      LMINWS = .FALSE.
      IF( ( TSIZE.LT.MAX( 1, NB*N*NBLCKS + 5 ) .OR. LWORK.LT.LWREQ ) .AND. ( LWORK.GE.N ) .AND. ( TSIZE.GE.MINTSZ ) .AND. ( .NOT.LQUERY ) ) THEN
        IF( TSIZE.LT.MAX( 1, NB*N*NBLCKS + 5 ) ) THEN
          LMINWS = .TRUE.
          NB = 1
          MB = M
        END IF
        IF( LWORK.LT.LWREQ ) THEN
          LMINWS = .TRUE.
          NB = 1
        END IF
      END IF

      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -4
      ELSE IF( TSIZE.LT.MAX( 1, NB*N*NBLCKS + 5 ) .AND. ( .NOT.LQUERY ) .AND. ( .NOT.LMINWS ) ) THEN
        INFO = -6
      ELSE IF( ( LWORK.LT.LWREQ ) .AND. ( .NOT.LQUERY ) .AND. ( .NOT.LMINWS ) ) THEN
        INFO = -8
      END IF

      IF( INFO.EQ.0 ) THEN
        IF( MINT ) THEN
          T( 1 ) = MINTSZ
        ELSE
          T( 1 ) = NB*N*NBLCKS + 5
        END IF
        T( 2 ) = MB
        T( 3 ) = NB
        IF( MINW ) THEN
          WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
        ELSE
          WORK( 1 ) = SROUNDUP_LWORK( LWREQ )
        END IF
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'SGEQR', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF

      // Quick return if possible

      IF( MIN( M, N ).EQ.0 ) THEN
        RETURN
      END IF

      // The QR Decomposition

      IF( ( M.LE.N ) .OR. ( MB.LE.N ) .OR. ( MB.GE.M ) ) THEN
        CALL SGEQRT( M, N, NB, A, LDA, T( 6 ), NB, WORK, INFO )
      ELSE
        CALL SLATSQR( M, N, MB, NB, A, LDA, T( 6 ), NB, WORK, LWORK, INFO )
      END IF

      WORK( 1 ) = SROUNDUP_LWORK( LWREQ )

      RETURN

      // End of SGEQR

      END
