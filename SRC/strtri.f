      SUBROUTINE STRTRI( UPLO, DIAG, N, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIAG, UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               NOUNIT, UPPER;
      int                J, JB, NB, NN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      // EXTERNAL LSAME, ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL STRMM, STRSM, STRTI2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      }
      if ( INFO.NE.0 ) {
         xerbla('STRTRI', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Check for singularity if non-unit.

      if ( NOUNIT ) {
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO ) RETURN
   10    CONTINUE
         INFO = 0
      }

      // Determine the block size for this environment.

      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
      if ( NB.LE.1 .OR. NB.GE.N ) {

         // Use unblocked code

         strti2(UPLO, DIAG, N, A, LDA, INFO );
      } else {

         // Use blocked code

         if ( UPPER ) {

            // Compute inverse of upper triangular matrix

            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )

               // Compute rows 1:j-1 of current block column

               strmm('Left', 'Upper', 'No transpose', DIAG, J-1, JB, ONE, A, LDA, A( 1, J ), LDA )                CALL STRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1, JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA );

               // Compute inverse of current diagonal block

               strti2('Upper', DIAG, JB, A( J, J ), LDA, INFO );
   20       CONTINUE
         } else {

            // Compute inverse of lower triangular matrix

            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               if ( J+JB.LE.N ) {

                  // Compute rows j+jb:n of current block column

                  strmm('Left', 'Lower', 'No transpose', DIAG, N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA, A( J+JB, J ), LDA )                   CALL STRSM( 'Right', 'Lower', 'No transpose', DIAG, N-J-JB+1, JB, -ONE, A( J, J ), LDA, A( J+JB, J ), LDA );
               }

               // Compute inverse of current diagonal block

               strti2('Lower', DIAG, JB, A( J, J ), LDA, INFO );
   30       CONTINUE
         }
      }

      RETURN

      // End of STRTRI

      }
