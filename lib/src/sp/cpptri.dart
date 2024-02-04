      void cpptri(UPLO, N, AP, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      Complex            AP( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, JC, JJ, JJN;
      REAL               AJJ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- COMPLEX            CDOTC;
      // EXTERNAL lsame, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPR, CSSCAL, CTPMV, CTPTRI, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('CPPTRI', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Invert the triangular Cholesky factor U or L.

      ctptri(UPLO, 'Non-unit', N, AP, INFO );
      if (INFO > 0) return;
      if ( UPPER ) {

         // Compute the product inv(U) * inv(U)**H.

         JJ = 0;
         for (J = 1; J <= N; J++) { // 10
            JC = JJ + 1;
            JJ = JJ + J;
            if (J > 1) chpr( 'Upper', J-1, ONE, AP( JC ), 1, AP );
            AJJ = REAL( AP( JJ ) );
            csscal(J, AJJ, AP( JC ), 1 );
         } // 10

      } else {

         // Compute the product inv(L)**H * inv(L).

         JJ = 1;
         for (J = 1; J <= N; J++) { // 20
            JJN = JJ + N - J + 1;
            AP[JJ] = REAL( CDOTC( N-J+1, AP( JJ ), 1, AP( JJ ), 1 ) );
            if (J < N) ctpmv( 'Lower', 'Conjugate transpose', 'Non-unit', N-J, AP( JJN ), AP( JJ+1 ), 1 );
            JJ = JJN;
         } // 20
      }

      return;
      }
