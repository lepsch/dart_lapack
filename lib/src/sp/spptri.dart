      void spptri(UPLO, N, AP, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, N;
      double               AP( * );
      // ..

      double               ONE;
      const              ONE = 1.0 ;
      bool               UPPER;
      int                J, JC, JJ, JJN;
      double               AJJ;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SDOT;
      // EXTERNAL lsame, SDOT
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSPR, STPMV, STPTRI, XERBLA

      // Test the input parameters.

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      }
      if ( INFO != 0 ) {
         xerbla('SPPTRI', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Invert the triangular Cholesky factor U or L.

      stptri(UPLO, 'Non-unit', N, AP, INFO );
      if (INFO > 0) return;

      if ( UPPER ) {

         // Compute the product inv(U) * inv(U)**T.

         JJ = 0;
         for (J = 1; J <= N; J++) { // 10
            JC = JJ + 1;
            JJ = JJ + J;
            if (J > 1) sspr( 'Upper', J-1, ONE, AP( JC ), 1, AP );
            AJJ = AP( JJ );
            sscal(J, AJJ, AP( JC ), 1 );
         } // 10

      } else {

         // Compute the product inv(L)**T * inv(L).

         JJ = 1;
         for (J = 1; J <= N; J++) { // 20
            JJN = JJ + N - J + 1;
            AP[JJ] = SDOT( N-J+1, AP( JJ ), 1, AP( JJ ), 1 );
            if (J < N) stpmv( 'Lower', 'Transpose', 'Non-unit', N-J, AP( JJN ), AP( JJ+1 ), 1 );
            JJ = JJN;
         } // 20
      }

      }
