      RECURSIVE SUBROUTINE DPOTRF2( UPLO, N, A, LDA, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                N1, N2, IINFO;
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSYRK, DTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRF2', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // N=1 case

      IF( N.EQ.1 ) THEN

         // Test for non-positive-definiteness

         IF( A( 1, 1 ).LE.ZERO.OR.DISNAN( A( 1, 1 ) ) ) THEN
            INFO = 1
            RETURN
         END IF

         // Factor

         A( 1, 1 ) = SQRT( A( 1, 1 ) )

      // Use recursive code

      } else {
         N1 = N/2
         N2 = N-N1

         // Factor A11

         CALL DPOTRF2( UPLO, N1, A( 1, 1 ), LDA, IINFO )
         IF ( IINFO.NE.0 ) THEN
            INFO = IINFO
            RETURN
         END IF

         // Compute the Cholesky factorization A = U**T*U

         IF( UPPER ) THEN

            // Update and scale A12

            CALL DTRSM( 'L', 'U', 'T', 'N', N1, N2, ONE, A( 1, 1 ), LDA, A( 1, N1+1 ), LDA )

            // Update and factor A22

            CALL DSYRK( UPLO, 'T', N2, N1, -ONE, A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA )
            CALL DPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
            IF ( IINFO.NE.0 ) THEN
               INFO = IINFO + N1
               RETURN
            END IF

         // Compute the Cholesky factorization A = L*L**T

         } else {

            // Update and scale A21

            CALL DTRSM( 'R', 'L', 'T', 'N', N2, N1, ONE, A( 1, 1 ), LDA, A( N1+1, 1 ), LDA )

            // Update and factor A22

            CALL DSYRK( UPLO, 'N', N2, N1, -ONE, A( N1+1, 1 ), LDA, ONE, A( N1+1, N1+1 ), LDA )
            CALL DPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
            IF ( IINFO.NE.0 ) THEN
               INFO = IINFO + N1
               RETURN
            END IF
         END IF
      END IF
      RETURN

      // End of DPOTRF2

      }
