      SUBROUTINE CGELQ( M, N, A, LDA, T, TSIZE, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N, TSIZE, LWORK;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), T( * ), WORK( * )
      // ..

*  =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LQUERY, LMINWS, MINT, MINW;
      int                MB, NB, MINTSZ, NBLCKS, LWMIN, LWOPT, LWREQ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGELQT, CLASWLQ, XERBLA
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
      if ( TSIZE.EQ.-2 .OR. LWORK.EQ.-2 ) {
        IF( TSIZE.NE.-1 ) MINT = .TRUE.
        IF( LWORK.NE.-1 ) MINW = .TRUE.
      }

      // Determine the block size

      if ( MIN( M, N ).GT.0 ) {
        MB = ILAENV( 1, 'CGELQ ', ' ', M, N, 1, -1 )
        NB = ILAENV( 1, 'CGELQ ', ' ', M, N, 2, -1 )
      } else {
        MB = 1
        NB = N
      }
      IF( MB.GT.MIN( M, N ) .OR. MB.LT.1 ) MB = 1
      IF( NB.GT.N .OR. NB.LE.M ) NB = N
      MINTSZ = M + 5
      if ( NB.GT.M .AND. N.GT.M ) {
        if ( MOD( N - M, NB - M ).EQ.0 ) {
          NBLCKS = ( N - M ) / ( NB - M )
        } else {
          NBLCKS = ( N - M ) / ( NB - M ) + 1
        }
      } else {
        NBLCKS = 1
      }

      // Determine if the workspace size satisfies minimal size

      if ( ( N.LE.M ) .OR. ( NB.LE.M ) .OR. ( NB.GE.N ) ) {
         LWMIN = MAX( 1, N )
         LWOPT = MAX( 1, MB*N )
      } else {
         LWMIN = MAX( 1, M )
         LWOPT = MAX( 1, MB*M )
      }
      LMINWS = .FALSE.
      if ( ( TSIZE.LT.MAX( 1, MB*M*NBLCKS + 5 ) .OR. LWORK.LT.LWOPT ) .AND. ( LWORK.GE.LWMIN ) .AND. ( TSIZE.GE.MINTSZ ) .AND. ( .NOT.LQUERY ) ) {
        if ( TSIZE.LT.MAX( 1, MB*M*NBLCKS + 5 ) ) {
          LMINWS = .TRUE.
          MB = 1
          NB = N
        }
        if ( LWORK.LT.LWOPT ) {
          LMINWS = .TRUE.
          MB = 1
        }
      }
      if ( ( N.LE.M ) .OR. ( NB.LE.M ) .OR. ( NB.GE.N ) ) {
         LWREQ = MAX( 1, MB*N )
      } else {
         LWREQ = MAX( 1, MB*M )
      }

      if ( M.LT.0 ) {
        INFO = -1
      } else if ( N.LT.0 ) {
        INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
        INFO = -4
      } else if ( TSIZE.LT.MAX( 1, MB*M*NBLCKS + 5 ) .AND. ( .NOT.LQUERY ) .AND. ( .NOT.LMINWS ) ) {
        INFO = -6
      } else if ( ( LWORK.LT.LWREQ ) .AND .( .NOT.LQUERY ) .AND. ( .NOT.LMINWS ) ) {
        INFO = -8
      }

      if ( INFO.EQ.0 ) {
        if ( MINT ) {
          T( 1 ) = MINTSZ
        } else {
          T( 1 ) = MB*M*NBLCKS + 5
        }
        T( 2 ) = MB
        T( 3 ) = NB
        if ( MINW ) {
          WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
        } else {
          WORK( 1 ) = SROUNDUP_LWORK( LWREQ )
        }
      }
      if ( INFO.NE.0 ) {
        xerbla('CGELQ', -INFO );
        RETURN
      } else if ( LQUERY ) {
        RETURN
      }

      // Quick return if possible

      if ( MIN( M, N ).EQ.0 ) {
        RETURN
      }

      // The LQ Decomposition

      if ( ( N.LE.M ) .OR. ( NB.LE.M ) .OR. ( NB.GE.N ) ) {
        cgelqt(M, N, MB, A, LDA, T( 6 ), MB, WORK, INFO );
      } else {
        claswlq(M, N, MB, NB, A, LDA, T( 6 ), MB, WORK, LWORK, INFO );
      }

      WORK( 1 ) = SROUNDUP_LWORK( LWREQ )

      RETURN

      // End of CGELQ

      }
