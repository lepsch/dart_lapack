      SUBROUTINE DGEQR( M, N, A, LDA, T, TSIZE, WORK, LWORK, INFO )

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
      int                MB, NB, MINTSZ, NBLCKS, LWMIN, LWREQ;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLATSQR, DGEQRT, XERBLA
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

      MINT = false;
      MINW = false;
      if ( TSIZE.EQ.-2 .OR. LWORK.EQ.-2 ) {
        if (TSIZE.NE.-1) MINT = true ;
        if (LWORK.NE.-1) MINW = true ;
      }

      // Determine the block size

      if ( MIN( M, N ).GT.0 ) {
        MB = ILAENV( 1, 'DGEQR ', ' ', M, N, 1, -1 )
        NB = ILAENV( 1, 'DGEQR ', ' ', M, N, 2, -1 )
      } else {
        MB = M
        NB = 1
      }
      if (MB.GT.M .OR. MB.LE.N) MB = M;
      IF( NB.GT.MIN( M, N ) .OR. NB.LT.1 ) NB = 1
      MINTSZ = N + 5
      if ( MB.GT.N .AND. M.GT.N ) {
        if ( MOD( M - N, MB - N ).EQ.0 ) {
          NBLCKS = ( M - N ) / ( MB - N )
        } else {
          NBLCKS = ( M - N ) / ( MB - N ) + 1
        }
      } else {
        NBLCKS = 1
      }

      // Determine if the workspace size satisfies minimal size

      LWMIN = MAX( 1, N )
      LWREQ = MAX( 1, N*NB )
      LMINWS = false;
      if ( ( TSIZE.LT.MAX( 1, NB*N*NBLCKS + 5 ) .OR. LWORK.LT.LWREQ ) .AND. ( LWORK.GE.N ) .AND. ( TSIZE.GE.MINTSZ ) .AND. ( .NOT.LQUERY ) ) {
        if ( TSIZE.LT.MAX( 1, NB*N*NBLCKS + 5 ) ) {
          LMINWS = true;
          NB = 1
          MB = M
        }
        if ( LWORK.LT.LWREQ ) {
          LMINWS = true;
          NB = 1
        }
      }

      if ( M.LT.0 ) {
        INFO = -1
      } else if ( N.LT.0 ) {
        INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
        INFO = -4
      } else if ( TSIZE.LT.MAX( 1, NB*N*NBLCKS + 5 ) .AND. ( .NOT.LQUERY ) .AND. ( .NOT.LMINWS ) ) {
        INFO = -6
      } else if ( ( LWORK.LT.LWREQ ) .AND. ( .NOT.LQUERY ) .AND. ( .NOT.LMINWS ) ) {
        INFO = -8
      }

      if ( INFO.EQ.0 ) {
        if ( MINT ) {
          T( 1 ) = MINTSZ
        } else {
          T( 1 ) = NB*N*NBLCKS + 5
        }
        T( 2 ) = MB
        T( 3 ) = NB
        if ( MINW ) {
          WORK( 1 ) = LWMIN
        } else {
          WORK( 1 ) = LWREQ
        }
      }
      if ( INFO.NE.0 ) {
        xerbla('DGEQR', -INFO );
        RETURN
      } else if ( LQUERY ) {
        RETURN
      }

      // Quick return if possible

      if ( MIN( M, N ).EQ.0 ) {
        RETURN
      }

      // The QR Decomposition

      if ( ( M.LE.N ) .OR. ( MB.LE.N ) .OR. ( MB.GE.M ) ) {
        dgeqrt(M, N, NB, A, LDA, T( 6 ), NB, WORK, INFO );
      } else {
        dlatsqr(M, N, MB, NB, A, LDA, T( 6 ), NB, WORK, LWORK, INFO );
      }

      WORK( 1 ) = LWREQ

      RETURN

      // End of DGEQR

      }
