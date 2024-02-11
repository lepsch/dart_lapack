      void zggqrf(final int N, final int M, final int P, final Matrix<double> A, final int LDA, final int TAUA, final Matrix<double> B, final int LDB, final int TAUB, final Array<double> WORK, final int LWORK, final Box<int> INFO,) {

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
      // EXTERNAL XERBLA, ZGEQRF, ZGERQF, ZUNMQR
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN

      // Test the input parameters

      INFO = 0;
      NB1 = ilaenv( 1, 'ZGEQRF', ' ', N, M, -1, -1 );
      NB2 = ilaenv( 1, 'ZGERQF', ' ', N, P, -1, -1 );
      NB3 = ilaenv( 1, 'ZUNMQR', ' ', N, M, P, -1 );
      NB = max( NB1, NB2, NB3 );
      LWKOPT = max( 1, max( N, M, P )*NB );
      WORK[1] = LWKOPT;
      LQUERY = ( LWORK == -1 );
      if ( N < 0 ) {
         INFO = -1;
      } else if ( M < 0 ) {
         INFO = -2;
      } else if ( P < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -8;
      } else if ( LWORK < max( 1, N, M, P ) && !LQUERY ) {
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

      zunmqr('Left', 'Conjugate Transpose', N, P, min( N, M ), A, LDA, TAUA, B, LDB, WORK, LWORK, INFO );
      LOPT = max( LOPT, INT( WORK( 1 ) ) );

      // RQ factorization of N-by-P matrix B: B = T*Z.

      zgerqf(N, P, B, LDB, TAUB, WORK, LWORK, INFO );
      WORK[1] = max( LOPT, INT( WORK( 1 ) ) );

      }
