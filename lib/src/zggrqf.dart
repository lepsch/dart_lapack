      void zggrqf(M, P, N, final Matrix<double> A, final int LDA, TAUA, final Matrix<double> B, final int LDB, TAUB, final Array<double> WORK, final int LWORK, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDB, LWORK, M, N, P;
      Complex         A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ), WORK( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LOPT, LWKOPT, NB, NB1, NB2, NB3;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQRF, ZGERQF, ZUNMRQ
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN

      // Test the input parameters

      INFO = 0;
      NB1 = ilaenv( 1, 'ZGERQF', ' ', M, N, -1, -1 );
      NB2 = ilaenv( 1, 'ZGEQRF', ' ', P, N, -1, -1 );
      NB3 = ilaenv( 1, 'ZUNMRQ', ' ', M, N, P, -1 );
      NB = max( NB1, NB2, NB3 );
      LWKOPT = max( 1, max( N, M, P )*NB );
      WORK[1] = LWKOPT;
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
         xerbla('ZGGRQF', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // RQ factorization of M-by-N matrix A: A = R*Q

      zgerqf(M, N, A, LDA, TAUA, WORK, LWORK, INFO );
      LOPT = INT( WORK( 1 ) );

      // Update B := B*Q**H

      zunmrq('Right', 'Conjugate Transpose', P, N, min( M, N ), A( max( 1, M-N+1 ), 1 ), LDA, TAUA, B, LDB, WORK, LWORK, INFO );
      LOPT = max( LOPT, INT( WORK( 1 ) ) );

      // QR factorization of P-by-N matrix B: B = Z*T

      zgeqrf(P, N, B, LDB, TAUB, WORK, LWORK, INFO );
      WORK[1] = max( LOPT, INT( WORK( 1 ) ) );

      }
