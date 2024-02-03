      void zhetrs_aa_2stage(UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, IPIV2, B, LDB, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // IMPLICIT NONE

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, NRHS, LDA, LTB, LDB, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IPIV2( * );
      Complex         A( LDA, * ), TB( * ), B( LDB, * );
      // ..

// =====================================================================

      Complex         ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                LDTB, NB;
      bool               UPPER;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGBTRS, ZLASWP, ZTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( NRHS < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LTB < ( 4*N ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('ZHETRS_AA_2STAGE', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      // Read NB and compute LDTB

      NB = INT( TB( 1 ) );
      LDTB = LTB/N;

      if ( UPPER ) {

         // Solve A*X = B, where A = U**H*T*U.

         if ( N > NB ) {

            // Pivot, P**T * B -> B

            zlaswp(NRHS, B, LDB, NB+1, N, IPIV, 1 );

            // Compute (U**H \ B) -> B    [ (U**H \P**T * B) ]

            ztrsm('L', 'U', 'C', 'U', N-NB, NRHS, ONE, A(1, NB+1), LDA, B(NB+1, 1), LDB);

         }

         // Compute T \ B -> B   [ T \ (U**H \P**T * B) ]

         zgbtrs('N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB, INFO);
         if ( N > NB ) {

            // Compute (U \ B) -> B   [ U \ (T \ (U**H \P**T * B) ) ]

            ztrsm('L', 'U', 'N', 'U', N-NB, NRHS, ONE, A(1, NB+1), LDA, B(NB+1, 1), LDB);

            // Pivot, P * B -> B  [ P * (U \ (T \ (U**H \P**T * B) )) ]

            zlaswp(NRHS, B, LDB, NB+1, N, IPIV, -1 );

         }

      } else {

         // Solve A*X = B, where A = L*T*L**H.

         if ( N > NB ) {

            // Pivot, P**T * B -> B

            zlaswp(NRHS, B, LDB, NB+1, N, IPIV, 1 );

            // Compute (L \ B) -> B    [ (L \P**T * B) ]

            ztrsm('L', 'L', 'N', 'U', N-NB, NRHS, ONE, A(NB+1, 1), LDA, B(NB+1, 1), LDB);

         }

         // Compute T \ B -> B   [ T \ (L \P**T * B) ]

         zgbtrs('N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB, INFO);
         if ( N > NB ) {

            // Compute (L**H \ B) -> B   [ L**H \ (T \ (L \P**T * B) ) ]

            ztrsm('L', 'L', 'C', 'U', N-NB, NRHS, ONE, A(NB+1, 1), LDA, B(NB+1, 1), LDB);

            // Pivot, P * B -> B  [ P * (L**H \ (T \ (L \P**T * B) )) ]

            zlaswp(NRHS, B, LDB, NB+1, N, IPIV, -1 );

         }
      }

      return;
      }
