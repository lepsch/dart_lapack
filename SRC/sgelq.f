      SUBROUTINE SGELQ( M, N, A, LDA, T, TSIZE, WORK, LWORK, INFO )

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
      int                MB, NB, MINTSZ, NBLCKS, LWMIN, LWOPT, LWREQ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGELQT, SLASWLQ, XERBLA
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

      LQUERY = ( TSIZE == -1 || TSIZE == -2 || LWORK == -1 || LWORK == -2 )

      MINT = false;
      MINW = false;
      if ( TSIZE == -2 || LWORK == -2 ) {
        if (TSIZE != -1) MINT = true ;
        if (LWORK != -1) MINW = true ;
      }

      // Determine the block size

      if ( MIN( M, N ) > 0 ) {
        MB = ILAENV( 1, 'SGELQ ', ' ', M, N, 1, -1 )
        NB = ILAENV( 1, 'SGELQ ', ' ', M, N, 2, -1 )
      } else {
        MB = 1
        NB = N
      }
      IF( MB > MIN( M, N ) || MB < 1 ) MB = 1
      if (NB > N || NB <= M) NB = N;
      MINTSZ = M + 5
      if ( NB > M && N > M ) {
        if ( MOD( N - M, NB - M ) == 0 ) {
          NBLCKS = ( N - M ) / ( NB - M )
        } else {
          NBLCKS = ( N - M ) / ( NB - M ) + 1
        }
      } else {
        NBLCKS = 1
      }

      // Determine if the workspace size satisfies minimal size

      if ( ( N <= M ) || ( NB <= M ) || ( NB >= N ) ) {
         LWMIN = MAX( 1, N )
         LWOPT = MAX( 1, MB*N )
      } else {
         LWMIN = MAX( 1, M )
         LWOPT = MAX( 1, MB*M )
      }
      LMINWS = false;
      if ( ( TSIZE < MAX( 1, MB*M*NBLCKS + 5 ) || LWORK < LWOPT ) && ( LWORK >= LWMIN ) && ( TSIZE >= MINTSZ ) && ( .NOT.LQUERY ) ) {
        if ( TSIZE < MAX( 1, MB*M*NBLCKS + 5 ) ) {
          LMINWS = true;
          MB = 1
          NB = N
        }
        if ( LWORK < LWOPT ) {
          LMINWS = true;
          MB = 1
        }
      }
      if ( ( N <= M ) || ( NB <= M ) || ( NB >= N ) ) {
         LWREQ = MAX( 1, MB*N )
      } else {
         LWREQ = MAX( 1, MB*M )
      }

      if ( M < 0 ) {
        INFO = -1
      } else if ( N < 0 ) {
        INFO = -2
      } else if ( LDA < MAX( 1, M ) ) {
        INFO = -4
      } else if ( TSIZE < MAX( 1, MB*M*NBLCKS + 5 ) && ( .NOT.LQUERY ) && ( .NOT.LMINWS ) ) {
        INFO = -6
      } else if ( ( LWORK < LWREQ ) .AND .( .NOT.LQUERY ) && ( .NOT.LMINWS ) ) {
        INFO = -8
      }

      if ( INFO == 0 ) {
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
      if ( INFO != 0 ) {
        xerbla('SGELQ', -INFO );
        RETURN
      } else if ( LQUERY ) {
        RETURN
      }

      // Quick return if possible

      if ( MIN( M, N ) == 0 ) {
        RETURN
      }

      // The LQ Decomposition

      if ( ( N <= M ) || ( NB <= M ) || ( NB >= N ) ) {
        sgelqt(M, N, MB, A, LDA, T( 6 ), MB, WORK, INFO );
      } else {
        slaswlq(M, N, MB, NB, A, LDA, T( 6 ), MB, WORK, LWORK, INFO );
      }

      WORK( 1 ) = SROUNDUP_LWORK( LWREQ )
      RETURN

      // End of SGELQ

      }
