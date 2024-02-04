      void zlaqhe(UPLO, N, A, LDA, S, SCOND, AMAX, EQUED ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, UPLO;
      int                LDA, N;
      double             AMAX, SCOND;
      // ..
      // .. Array Arguments ..
      double             S( * );
      Complex         A( LDA, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, THRESH;
      const              ONE = 1.0, THRESH = 0.1 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH;
      // EXTERNAL lsame, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N <= 0 ) {
         EQUED = 'N';
         return;
      }

      // Initialize LARGE and SMALL.

      SMALL = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' );
      LARGE = ONE / SMALL;

      if ( SCOND >= THRESH && AMAX >= SMALL && AMAX <= LARGE ) {

         // No equilibration

         EQUED = 'N';
      } else {

         // Replace A by diag(S) * A * diag(S).

         if ( lsame( UPLO, 'U' ) ) {

            // Upper triangle of A is stored.

            for (J = 1; J <= N; J++) { // 20
               CJ = S( J );
               for (I = 1; I <= J - 1; I++) { // 10
                  A[I, J] = CJ*S( I )*A( I, J );
               } // 10
               A[J, J] = CJ*CJ*(A( J, J )).toDouble();
            } // 20
         } else {

            // Lower triangle of A is stored.

            for (J = 1; J <= N; J++) { // 40
               CJ = S( J );
               A[J, J] = CJ*CJ*(A( J, J )).toDouble();
               for (I = J + 1; I <= N; I++) { // 30
                  A[I, J] = CJ*S( I )*A( I, J );
               } // 30
            } // 40
         }
         EQUED = 'Y';
      }

      return;
      }