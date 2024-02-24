      void zgelq(final int M, final int N, final Matrix<double> A_, final int LDA, final int T, final int TSIZE, final Array<double> WORK_, final int LWORK, final Box<int> INFO,) {
  final A = A_.dim();
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
      int                INFO, LDA, M, N, TSIZE, LWORK;
      Complex         A( LDA, * ), T( * ), WORK( * );
      // ..

// =====================================================================

      bool               LQUERY, LMINWS, MINT, MINW;
      int                MB, NB, MINTSZ, NBLCKS, LWMIN, LWOPT, LWREQ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGELQT, ZLASWLQ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, MOD
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV

      // Test the input arguments

      INFO = 0;

      LQUERY = ( TSIZE == -1 || TSIZE == -2 || LWORK == -1 || LWORK == -2 );

      MINT = false;
      MINW = false;
      if ( TSIZE == -2 || LWORK == -2 ) {
        if (TSIZE != -1) MINT = true ;
        if (LWORK != -1) MINW = true ;
      }

      // Determine the block size

      if ( min( M, N ) > 0 ) {
        MB = ilaenv( 1, 'ZGELQ ', ' ', M, N, 1, -1 );
        NB = ilaenv( 1, 'ZGELQ ', ' ', M, N, 2, -1 );
      } else {
        MB = 1;
        NB = N;
      }
      if( MB > min( M, N ) || MB < 1 ) MB = 1;
      if (NB > N || NB <= M) NB = N;
      MINTSZ = M + 5;
      if ( NB > M && N > M ) {
        if ( ((N - M) % (NB - M)) == 0 ) {
          NBLCKS = ( N - M ) / ( NB - M );
        } else {
          NBLCKS = ( N - M ) / ( NB - M ) + 1;
        }
      } else {
        NBLCKS = 1;
      }

      // Determine if the workspace size satisfies minimal size

      if ( ( N <= M ) || ( NB <= M ) || ( NB >= N ) ) {
         LWMIN = max( 1, N );
         LWOPT = max( 1, MB*N );
      } else {
         LWMIN = max( 1, M );
         LWOPT = max( 1, MB*M );
      }
      LMINWS = false;
      if ( ( TSIZE < max( 1, MB*M*NBLCKS + 5 ) || LWORK < LWOPT ) && ( LWORK >= LWMIN ) && ( TSIZE >= MINTSZ ) && ( !LQUERY ) ) {
        if ( TSIZE < max( 1, MB*M*NBLCKS + 5 ) ) {
            LMINWS = true;
            MB = 1;
            NB = N;
        }
        if ( LWORK < LWOPT ) {
            LMINWS = true;
            MB = 1;
        }
      }
      if ( ( N <= M ) || ( NB <= M ) || ( NB >= N ) ) {
         LWREQ = max( 1, MB*N );
      } else {
         LWREQ = max( 1, MB*M );
      }

      if ( M < 0 ) {
        INFO = -1;
      } else if ( N < 0 ) {
        INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
        INFO = -4;
      } else if ( TSIZE < max( 1, MB*M*NBLCKS + 5 ) && ( !LQUERY ) && ( !LMINWS ) ) {
        INFO = -6;
      } else if ( ( LWORK < LWREQ ) .AND .( !LQUERY ) && ( !LMINWS ) ) {
        INFO = -8;
      }

      if ( INFO == 0 ) {
        if ( MINT ) {
          T[1] = MINTSZ;
        } else {
          T[1] = MB*M*NBLCKS + 5;
        }
        T[2] = MB;
        T[3] = NB;
        if ( MINW ) {
          WORK[1] = LWMIN;
        } else {
          WORK[1] = LWREQ;
        }
      }
      if ( INFO != 0 ) {
        xerbla('ZGELQ', -INFO );
        return;
      } else if ( LQUERY ) {
        return;
      }

      // Quick return if possible

      if ( min( M, N ) == 0 ) {
        return;
      }

      // The LQ Decomposition

      if ( ( N <= M ) || ( NB <= M ) || ( NB >= N ) ) {
        zgelqt(M, N, MB, A, LDA, T( 6 ), MB, WORK, INFO );
      } else {
        zlaswlq(M, N, MB, NB, A, LDA, T( 6 ), MB, WORK, LWORK, INFO );
      }

      WORK[1] = LWREQ;

      }
