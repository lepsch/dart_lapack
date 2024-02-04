      void cggrqf(M, P, N, A, LDA, TAUA, B, LDB, TAUB, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      Complex            A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LOPT, LWKOPT, NB, NB1, NB2, NB3;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQRF, CGERQF, CUNMRQ, XERBLA
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0;
      NB1 = ILAENV( 1, 'CGERQF', ' ', M, N, -1, -1 );
      NB2 = ILAENV( 1, 'CGEQRF', ' ', P, N, -1, -1 );
      NB3 = ILAENV( 1, 'CUNMRQ', ' ', M, N, P, -1 );
      NB = max( NB1, NB2, NB3 );
      LWKOPT = max( 1, max( N, M, P )*NB );
      WORK[1] = SROUNDUP_LWORK( LWKOPT );
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( P < 0 ) {
         INFO = -2;
      } else if ( N < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, P ) ) {
         INFO = -8;
      } else if ( LWORK < max( 1, M, P, N ) && !LQUERY ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('CGGRQF', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // RQ factorization of M-by-N matrix A: A = R*Q

      cgerqf(M, N, A, LDA, TAUA, WORK, LWORK, INFO );
      LOPT = INT( WORK( 1 ) );

      // Update B := B*Q**H

      cunmrq('Right', 'Conjugate Transpose', P, N, min( M, N ), A( max( 1, M-N+1 ), 1 ), LDA, TAUA, B, LDB, WORK, LWORK, INFO );
      LOPT = max( LOPT, INT( WORK( 1 ) ) );

      // QR factorization of P-by-N matrix B: B = Z*T

      cgeqrf(P, N, B, LDB, TAUB, WORK, LWORK, INFO );
      WORK[1] = SROUNDUP_LWORK( max( LOPT, INT( WORK( 1 ) ) ) );

      return;
      }