      void ztrexc(COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ;
      int                IFST, ILST, INFO, LDQ, LDT, N;
      // ..
      // .. Array Arguments ..
      Complex         Q( LDQ, * ), T( LDT, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               WANTQ;
      int                K, M1, M2, M3;
      double             CS;
      Complex         SN, T11, T22, TEMP;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZLARTG, ZROT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCONJG, MAX
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters.

      INFO = 0;
      WANTQ = lsame( COMPQ, 'V' );
      if ( !lsame( COMPQ, 'N' ) && !WANTQ ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDT < max( 1, N ) ) {
         INFO = -4;
      } else if ( LDQ < 1 || ( WANTQ && LDQ < max( 1, N ) ) ) {
         INFO = -6;
      } else if (( IFST < 1 || IFST > N ) && ( N > 0 )) {
         INFO = -7;
      } else if (( ILST < 1 || ILST > N ) && ( N > 0 )) {
         INFO = -8;
      }
      if ( INFO != 0 ) {
         xerbla('ZTREXC', -INFO );
         return;
      }

      // Quick return if possible

      if (N <= 1 || IFST == ILST) return;

      if ( IFST < ILST ) {

         // Move the IFST-th diagonal element forward down the diagonal.

         M1 = 0;
         M2 = -1;
         M3 = 1;
      } else {

         // Move the IFST-th diagonal element backward up the diagonal.

         M1 = -1;
         M2 = 0;
         M3 = -1;
      }

      for (K = IFST + M1; M3 < 0 ? K >= ILST + M2 : K <= ILST + M2; K += M3) { // 10

         // Interchange the k-th and (k+1)-th diagonal elements.

         T11 = T( K, K );
         T22 = T( K+1, K+1 );

         // Determine the transformation to perform the interchange.

         zlartg(T( K, K+1 ), T22-T11, CS, SN, TEMP );

         // Apply transformation to the matrix T.

         if (K+2 <= N) zrot( N-K-1, T( K, K+2 ), LDT, T( K+1, K+2 ), LDT, CS, SN );
         zrot(K-1, T( 1, K ), 1, T( 1, K+1 ), 1, CS, DCONJG( SN ) );

         T[K, K] = T22;
         T[K+1, K+1] = T11;

         if ( WANTQ ) {

            // Accumulate transformation in the matrix Q.

            zrot(N, Q( 1, K ), 1, Q( 1, K+1 ), 1, CS, DCONJG( SN ) );
         }

      } // 10

      return;
      }