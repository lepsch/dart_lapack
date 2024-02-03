      SUBROUTINE DGELQ( M, N, A, LDA, T, TSIZE, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N, TSIZE, LWORK;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), T( * ), WORK( * );
      // ..

*  =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LQUERY, LMINWS, MINT, MINW;
      int                MB, NB, MINTSZ, NBLCKS, LWMIN, LWOPT, LWREQ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGELQT, DLASWLQ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, MOD
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Executable Statements ..

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
        MB = ILAENV( 1, 'DGELQ ', ' ', M, N, 1, -1 )
        NB = ILAENV( 1, 'DGELQ ', ' ', M, N, 2, -1 )
      ELSE
        MB = 1
        NB = N
      END IF
      IF( MB.GT.MIN( M, N ) .OR. MB.LT.1 ) MB = 1
      IF( NB.GT.N .OR. NB.LE.M ) NB = N
      MINTSZ = M + 5
      IF ( NB.GT.M .AND. N.GT.M ) THEN
        IF( MOD( N - M, NB - M ).EQ.0 ) THEN
          NBLCKS = ( N - M ) / ( NB - M )
        ELSE
          NBLCKS = ( N - M ) / ( NB - M ) + 1
        END IF
      ELSE
        NBLCKS = 1
      END IF

      // Determine if the workspace size satisfies minimal size

      IF( ( N.LE.M ) .OR. ( NB.LE.M ) .OR. ( NB.GE.N ) ) THEN
         LWMIN = MAX( 1, N )
         LWOPT = MAX( 1, MB*N )
      ELSE
         LWMIN = MAX( 1, M )
         LWOPT = MAX( 1, MB*M )
      END IF
      LMINWS = .FALSE.
      IF( ( TSIZE.LT.MAX( 1, MB*M*NBLCKS + 5 ) .OR. LWORK.LT.LWOPT ) .AND. ( LWORK.GE.LWMIN ) .AND. ( TSIZE.GE.MINTSZ ) .AND. ( .NOT.LQUERY ) ) THEN
        IF( TSIZE.LT.MAX( 1, MB*M*NBLCKS + 5 ) ) THEN
            LMINWS = .TRUE.
            MB = 1
            NB = N
        END IF
        IF( LWORK.LT.LWOPT ) THEN
            LMINWS = .TRUE.
            MB = 1
        END IF
      END IF
      IF( ( N.LE.M ) .OR. ( NB.LE.M ) .OR. ( NB.GE.N ) ) THEN
         LWREQ = MAX( 1, MB*N )
      ELSE
         LWREQ = MAX( 1, MB*M )
      END IF

      IF( M.LT.0 ) THEN
        INFO = -1
      ELSE IF( N.LT.0 ) THEN
        INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
        INFO = -4
      ELSE IF( TSIZE.LT.MAX( 1, MB*M*NBLCKS + 5 ) .AND. ( .NOT.LQUERY ) .AND. ( .NOT.LMINWS ) ) THEN
        INFO = -6
      ELSE IF( ( LWORK.LT.LWREQ ) .AND .( .NOT.LQUERY ) .AND. ( .NOT.LMINWS ) ) THEN
        INFO = -8
      END IF

      IF( INFO.EQ.0 ) THEN
        IF( MINT ) THEN
          T( 1 ) = MINTSZ
        ELSE
          T( 1 ) = MB*M*NBLCKS + 5
        END IF
        T( 2 ) = MB
        T( 3 ) = NB
        IF( MINW ) THEN
          WORK( 1 ) = LWMIN
        ELSE
          WORK( 1 ) = LWREQ
        END IF
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DGELQ', -INFO )
        RETURN
      ELSE IF( LQUERY ) THEN
        RETURN
      END IF

      // Quick return if possible

      IF( MIN( M, N ).EQ.0 ) THEN
        RETURN
      END IF

      // The LQ Decomposition

      IF( ( N.LE.M ) .OR. ( NB.LE.M ) .OR. ( NB.GE.N ) ) THEN
        CALL DGELQT( M, N, MB, A, LDA, T( 6 ), MB, WORK, INFO )
      ELSE
        CALL DLASWLQ( M, N, MB, NB, A, LDA, T( 6 ), MB, WORK, LWORK, INFO )
      END IF

      WORK( 1 ) = LWREQ

      RETURN

      // End of DGELQ

      END
