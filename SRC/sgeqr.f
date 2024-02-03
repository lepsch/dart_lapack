      void sgeqr(M, N, A, LDA, T, TSIZE, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N, TSIZE, LWORK;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), T( * ), WORK( * );
      // ..

// =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LQUERY, LMINWS, MINT, MINW;
      int                MB, NB, MINTSZ, NBLCKS, LWMIN, LWREQ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      REAL               SROUNDUP_LWORK;
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
        MB = ILAENV( 1, 'SGEQR ', ' ', M, N, 1, -1 );
        NB = ILAENV( 1, 'SGEQR ', ' ', M, N, 2, -1 );
      } else {
        MB = M;
        NB = 1;
      }
      if (MB > M || MB <= N) MB = M;
      if( NB > min( M, N ) || NB < 1 ) NB = 1;
      MINTSZ = N + 5;
      if ( MB > N && M > N ) {
        if ( MOD( M - N, MB - N ) == 0 ) {
          NBLCKS = ( M - N ) / ( MB - N );
        } else {
          NBLCKS = ( M - N ) / ( MB - N ) + 1;
        }
      } else {
        NBLCKS = 1;
      }

      // Determine if the workspace size satisfies minimal size

      LWMIN = max( 1, N );
      LWREQ = max( 1, N*NB );
      LMINWS = false;
      if ( ( TSIZE < max( 1, NB*N*NBLCKS + 5 ) || LWORK < LWREQ ) && ( LWORK >= N ) && ( TSIZE >= MINTSZ ) && ( !LQUERY ) ) {
        if ( TSIZE < max( 1, NB*N*NBLCKS + 5 ) ) {
          LMINWS = true;
          NB = 1;
          MB = M;
        }
        if ( LWORK < LWREQ ) {
          LMINWS = true;
          NB = 1;
        }
      }

      if ( M < 0 ) {
        INFO = -1;
      } else if ( N < 0 ) {
        INFO = -2;
      } else if ( LDA < max( 1, M ) ) {
        INFO = -4;
      } else if ( TSIZE < max( 1, NB*N*NBLCKS + 5 ) && ( !LQUERY ) && ( !LMINWS ) ) {
        INFO = -6;
      } else if ( ( LWORK < LWREQ ) && ( !LQUERY ) && ( !LMINWS ) ) {
        INFO = -8;
      }

      if ( INFO == 0 ) {
        if ( MINT ) {
          T( 1 ) = MINTSZ;
        } else {
          T( 1 ) = NB*N*NBLCKS + 5;
        }
        T( 2 ) = MB;
        T( 3 ) = NB;
        if ( MINW ) {
          WORK( 1 ) = SROUNDUP_LWORK( LWMIN );
        } else {
          WORK( 1 ) = SROUNDUP_LWORK( LWREQ );
        }
      }
      if ( INFO != 0 ) {
        xerbla('SGEQR', -INFO );
        return;
      } else if ( LQUERY ) {
        return;
      }

      // Quick return if possible

      if ( min( M, N ) == 0 ) {
        return;
      }

      // The QR Decomposition

      if ( ( M <= N ) || ( MB <= N ) || ( MB >= M ) ) {
        sgeqrt(M, N, NB, A, LDA, T( 6 ), NB, WORK, INFO );
      } else {
        slatsqr(M, N, MB, NB, A, LDA, T( 6 ), NB, WORK, LWORK, INFO );
      }

      WORK( 1 ) = SROUNDUP_LWORK( LWREQ );

      return;
      }
