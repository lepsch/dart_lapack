      SUBROUTINE CSYTRS_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, IPIV2, B, LDB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      String             UPLO;
      int                N, NRHS, LDA, LTB, LDB, INFO;
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IPIV2( * );
      COMPLEX            A( LDA, * ), TB( * ), B( LDB, * )
      // ..

*  =====================================================================

      COMPLEX            ONE
      const              ONE  = ( 1.0E+0, 0.0E+0 ) ;
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
      // EXTERNAL CGBTRS, CLASWP, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( NRHS.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LTB.LT.( 4*N ) ) {
         INFO = -7
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -11
      }
      if ( INFO.NE.0 ) {
         xerbla('CSYTRS_AA_2STAGE', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      // Read NB and compute LDTB

      NB = INT( TB( 1 ) )
      LDTB = LTB/N

      if ( UPPER ) {

         // Solve A*X = B, where A = U**T*T*U.

         if ( N.GT.NB ) {

            // Pivot, P**T * B -> B

            claswp(NRHS, B, LDB, NB+1, N, IPIV, 1 );

            // Compute (U**T \ B) -> B    [ (U**T \P**T * B) ]

            ctrsm('L', 'U', 'T', 'U', N-NB, NRHS, ONE, A(1, NB+1), LDA, B(NB+1, 1), LDB);

         }

         // Compute T \ B -> B   [ T \ (U**T \P**T * B) ]

         cgbtrs('N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB, INFO);
         if ( N.GT.NB ) {

            // Compute (U \ B) -> B   [ U \ (T \ (U**T \P**T * B) ) ]

            ctrsm('L', 'U', 'N', 'U', N-NB, NRHS, ONE, A(1, NB+1), LDA, B(NB+1, 1), LDB);

            // Pivot, P * B -> B  [ P * (U \ (T \ (U**T \P**T * B) )) ]

            claswp(NRHS, B, LDB, NB+1, N, IPIV, -1 );

         }

      } else {

         // Solve A*X = B, where A = L*T*L**T.

         if ( N.GT.NB ) {

            // Pivot, P**T * B -> B

            claswp(NRHS, B, LDB, NB+1, N, IPIV, 1 );

            // Compute (L \ B) -> B    [ (L \P**T * B) ]

            ctrsm('L', 'L', 'N', 'U', N-NB, NRHS, ONE, A(NB+1, 1), LDA, B(NB+1, 1), LDB);

         }

         // Compute T \ B -> B   [ T \ (L \P**T * B) ]

         cgbtrs('N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB, INFO);
         if ( N.GT.NB ) {

            // Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ]

            ctrsm('L', 'L', 'T', 'U', N-NB, NRHS, ONE, A(NB+1, 1), LDA, B(NB+1, 1), LDB);

            // Pivot, P * B -> B  [ P * (L**T \ (T \ (L \P**T * B) )) ]

            claswp(NRHS, B, LDB, NB+1, N, IPIV, -1 );

         }
      }

      RETURN

      // End of CSYTRS_AA_2STAGE

      }
