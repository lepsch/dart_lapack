      void claqp2(M, N, OFFSET, A, LDA, JPVT, TAU, VN1, VN2, WORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, M, N, OFFSET;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      double               VN1( * ), VN2( * );
      Complex            A( LDA, * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      Complex            CONE;
      const              ZERO = 0.0, ONE = 1.0, CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, ITEMP, J, MN, OFFPI, PVT;
      double               TEMP, TEMP2, TOL3Z;
      Complex            AII;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLARF, CLARFG, CSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, MAX, MIN, SQRT
      // ..
      // .. External Functions ..
      //- int                ISAMAX;
      //- REAL               SCNRM2, SLAMCH;
      // EXTERNAL ISAMAX, SCNRM2, SLAMCH
      // ..
      // .. Executable Statements ..

      MN = min( M-OFFSET, N );
      TOL3Z = sqrt(SLAMCH('Epsilon'));

      // Compute factorization.

      for (I = 1; I <= MN; I++) { // 20

         OFFPI = OFFSET + I;

         // Determine ith pivot column and swap if necessary.

         PVT = ( I-1 ) + ISAMAX( N-I+1, VN1( I ), 1 );

         if ( PVT != I ) {
            cswap(M, A( 1, PVT ), 1, A( 1, I ), 1 );
            ITEMP = JPVT( PVT );
            JPVT[PVT] = JPVT( I );
            JPVT[I] = ITEMP;
            VN1[PVT] = VN1( I );
            VN2[PVT] = VN2( I );
         }

         // Generate elementary reflector H(i).

         if ( OFFPI < M ) {
            clarfg(M-OFFPI+1, A( OFFPI, I ), A( OFFPI+1, I ), 1, TAU( I ) );
         } else {
            clarfg(1, A( M, I ), A( M, I ), 1, TAU( I ) );
         }

         if ( I < N ) {

            // Apply H(i)**H to A(offset+i:m,i+1:n) from the left.

            AII = A( OFFPI, I );
            A[OFFPI][I] = CONE;
            clarf('Left', M-OFFPI+1, N-I, A( OFFPI, I ), 1, CONJG( TAU( I ) ), A( OFFPI, I+1 ), LDA, WORK( 1 ) );
            A[OFFPI][I] = AII;
         }

         // Update partial column norms.

         for (J = I + 1; J <= N; J++) { // 10
            if ( VN1( J ) != ZERO ) {

               // NOTE: The following 4 lines follow from the analysis in
               // Lapack Working Note 176.

               TEMP = ONE - ( ( A( OFFPI, J ) ).abs() / VN1( J ) )**2;
               TEMP = max( TEMP, ZERO );
               TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2;
               if ( TEMP2 <= TOL3Z ) {
                  if ( OFFPI < M ) {
                     VN1[J] = SCNRM2( M-OFFPI, A( OFFPI+1, J ), 1 );
                     VN2[J] = VN1( J );
                  } else {
                     VN1[J] = ZERO;
                     VN2[J] = ZERO;
                  }
               } else {
                  VN1[J] = VN1( J )*sqrt( TEMP );
               }
            }
         } // 10

      } // 20

      return;
      }
