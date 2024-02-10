      void ssytrs_aa_2stage(final int UPLO, final int N, final int NRHS, final Matrix<double> A, final int LDA, final int TB, final int LTB, final Array<int> IPIV, final int IPIV2, final Matrix<double> B, final int LDB, final Box<int> INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      String             UPLO;
      int                N, NRHS, LDA, LTB, LDB, INFO;
      int                IPIV( * ), IPIV2( * );
      double               A( LDA, * ), TB( * ), B( LDB, * );
      // ..

// =====================================================================

      double               ONE;
      const              ONE = 1.0 ;
      int                LDTB, NB;
      bool               UPPER;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGBTRS, SLASWP, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
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
         xerbla('SSYTRS_AA_2STAGE', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0 || NRHS == 0) return;

      // Read NB and compute LDTB

      NB = INT( TB( 1 ) );
      LDTB = LTB/N;

      if ( UPPER ) {

         // Solve A*X = B, where A = U**T*T*U.

         if ( N > NB ) {

            // Pivot, P**T * B -> B

            slaswp(NRHS, B, LDB, NB+1, N, IPIV, 1 );

            // Compute (U**T \ B) -> B    [ (U**T \P**T * B) ]

            strsm('L', 'U', 'T', 'U', N-NB, NRHS, ONE, A(1, NB+1), LDA, B(NB+1, 1), LDB);

         }

         // Compute T \ B -> B   [ T \ (U**T \P**T * B) ]

         sgbtrs('N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB, INFO);
         if ( N > NB ) {

            // Compute (U \ B) -> B   [ U \ (T \ (U**T \P**T * B) ) ]

            strsm('L', 'U', 'N', 'U', N-NB, NRHS, ONE, A(1, NB+1), LDA, B(NB+1, 1), LDB);

            // Pivot, P * B -> B  [ P * (U \ (T \ (U**T \P**T * B) )) ]

            slaswp(NRHS, B, LDB, NB+1, N, IPIV, -1 );

         }

      } else {

         // Solve A*X = B, where A = L*T*L**T.

         if ( N > NB ) {

            // Pivot, P**T * B -> B

            slaswp(NRHS, B, LDB, NB+1, N, IPIV, 1 );

            // Compute (L \ B) -> B    [ (L \P**T * B) ]

            strsm('L', 'L', 'N', 'U', N-NB, NRHS, ONE, A(NB+1, 1), LDA, B(NB+1, 1), LDB);

         }

         // Compute T \ B -> B   [ T \ (L \P**T * B) ]

         sgbtrs('N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB, INFO);
         if ( N > NB ) {

            // Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ]

            strsm('L', 'L', 'T', 'U', N-NB, NRHS, ONE, A(NB+1, 1), LDA, B(NB+1, 1), LDB);

            // Pivot, P * B -> B  [ P * (L**T \ (T \ (L \P**T * B) )) ]

            slaswp(NRHS, B, LDB, NB+1, N, IPIV, -1 );

         }
      }

      }
