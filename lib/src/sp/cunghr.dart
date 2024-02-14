      void cunghr(final int N, final int ILO, final int IHI, final Matrix<double> A_, final int LDA, final int TAU, final Array<double> WORK_, final int LWORK, final Box<int> INFO,) {
  final A = A_.dim();
  final WORK = WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                IHI, ILO, INFO, LDA, LWORK, N;
      Complex            A( LDA, * ), TAU( * ), WORK( * );
      // ..

      Complex            ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) ;
      bool               LQUERY;
      int                I, IINFO, J, LWKOPT, NB, NH;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CUNGQR, XERBLA
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN

      // Test the input arguments

      INFO = 0;
      NH = IHI - ILO;
      LQUERY = ( LWORK == -1 );
      if ( N < 0 ) {
         INFO = -1;
      } else if ( ILO < 1 || ILO > max( 1, N ) ) {
         INFO = -2;
      } else if ( IHI < min( ILO, N ) || IHI > N ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LWORK < max( 1, NH ) && !LQUERY ) {
         INFO = -8;
      }

      if ( INFO == 0 ) {
         NB = ilaenv( 1, 'CUNGQR', ' ', NH, NH, NH, -1 );
         LWKOPT = max( 1, NH )*NB;
         WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }

      if ( INFO != 0 ) {
         xerbla('CUNGHR', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         WORK[1] = 1;
         return;
      }

      // Shift the vectors which define the elementary reflectors one
      // column to the right, and set the first ilo and the last n-ihi
      // rows and columns to those of the unit matrix

      for (J = IHI; J >= ILO + 1; J--) { // 40
         for (I = 1; I <= J - 1; I++) { // 10
            A[I][J] = ZERO;
         } // 10
         for (I = J + 1; I <= IHI; I++) { // 20
            A[I][J] = A( I, J-1 );
         } // 20
         for (I = IHI + 1; I <= N; I++) { // 30
            A[I][J] = ZERO;
         } // 30
      } // 40
      for (J = 1; J <= ILO; J++) { // 60
         for (I = 1; I <= N; I++) { // 50
            A[I][J] = ZERO;
         } // 50
         A[J][J] = ONE;
      } // 60
      for (J = IHI + 1; J <= N; J++) { // 80
         for (I = 1; I <= N; I++) { // 70
            A[I][J] = ZERO;
         } // 70
         A[J][J] = ONE;
      } // 80

      if ( NH > 0 ) {

         // Generate Q(ilo+1:ihi,ilo+1:ihi)

         cungqr(NH, NH, NH, A( ILO+1, ILO+1 ), LDA, TAU( ILO ), WORK, LWORK, IINFO );
      }
      WORK[1] = SROUNDUP_LWORK(LWKOPT);
      }
