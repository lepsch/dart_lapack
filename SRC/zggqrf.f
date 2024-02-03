      SUBROUTINE ZGGQRF( N, M, P, A, LDA, TAUA, B, LDB, TAUB, WORK, LWORK, INFO );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LOPT, LWKOPT, NB, NB1, NB2, NB3;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQRF, ZGERQF, ZUNMQR
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0;
      NB1 = ILAENV( 1, 'ZGEQRF', ' ', N, M, -1, -1 );
      NB2 = ILAENV( 1, 'ZGERQF', ' ', N, P, -1, -1 );
      NB3 = ILAENV( 1, 'ZUNMQR', ' ', N, M, P, -1 );
      NB = MAX( NB1, NB2, NB3 );
      LWKOPT = MAX( 1, MAX( N, M, P )*NB );
      WORK( 1 ) = LWKOPT;
      LQUERY = ( LWORK == -1 );
      if ( N < 0 ) {
         INFO = -1;
      } else if ( M < 0 ) {
         INFO = -2;
      } else if ( P < 0 ) {
         INFO = -3;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -8;
      } else if ( LWORK < MAX( 1, N, M, P ) && !LQUERY ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('ZGGQRF', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // QR factorization of N-by-M matrix A: A = Q*R

      zgeqrf(N, M, A, LDA, TAUA, WORK, LWORK, INFO );
      LOPT = INT( WORK( 1 ) );

      // Update B := Q**H*B.

      zunmqr('Left', 'Conjugate Transpose', N, P, MIN( N, M ), A, LDA, TAUA, B, LDB, WORK, LWORK, INFO );
      LOPT = MAX( LOPT, INT( WORK( 1 ) ) );

      // RQ factorization of N-by-P matrix B: B = T*Z.

      zgerqf(N, P, B, LDB, TAUB, WORK, LWORK, INFO );
      WORK( 1 ) = MAX( LOPT, INT( WORK( 1 ) ) );

      return;

      // End of ZGGQRF

      }
