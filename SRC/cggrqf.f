      SUBROUTINE CGGRQF( M, P, N, A, LDA, TAUA, B, LDB, TAUB, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LOPT, LWKOPT, NB, NB1, NB2, NB3;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQRF, CGERQF, CUNMRQ, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      NB1 = ILAENV( 1, 'CGERQF', ' ', M, N, -1, -1 )
      NB2 = ILAENV( 1, 'CGEQRF', ' ', P, N, -1, -1 )
      NB3 = ILAENV( 1, 'CUNMRQ', ' ', M, N, P, -1 )
      NB = MAX( NB1, NB2, NB3 )
      LWKOPT = MAX( 1, MAX( N, M, P )*NB )
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      LQUERY = ( LWORK == -1 )
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( P.LT.0 ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, P ) ) {
         INFO = -8
      } else if ( LWORK.LT.MAX( 1, M, P, N ) .AND. .NOT.LQUERY ) {
         INFO = -11
      }
      if ( INFO != 0 ) {
         xerbla('CGGRQF', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // RQ factorization of M-by-N matrix A: A = R*Q

      cgerqf(M, N, A, LDA, TAUA, WORK, LWORK, INFO );
      LOPT = INT( WORK( 1 ) )

      // Update B := B*Q**H

      cunmrq('Right', 'Conjugate Transpose', P, N, MIN( M, N ), A( MAX( 1, M-N+1 ), 1 ), LDA, TAUA, B, LDB, WORK, LWORK, INFO );
      LOPT = MAX( LOPT, INT( WORK( 1 ) ) )

      // QR factorization of P-by-N matrix B: B = Z*T

      cgeqrf(P, N, B, LDB, TAUB, WORK, LWORK, INFO );
      WORK( 1 ) = SROUNDUP_LWORK( MAX( LOPT, INT( WORK( 1 ) ) ) )

      RETURN

      // End of CGGRQF

      }
