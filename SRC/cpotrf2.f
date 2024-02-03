      RECURSIVE SUBROUTINE CPOTRF2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      COMPLEX            CONE
      PARAMETER          ( CONE = (1.0E+0, 0.0E+0) )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                N1, N2, IINFO;
      REAL               AJJ
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      // EXTERNAL LSAME, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHERK, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL, SQRT
      // ..
      // .. Executable Statements ..
*
      // Test the input parameters
*
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
         CALL XERBLA( 'CPOTRF2', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      // N=1 case
*
      IF( N.EQ.1 ) THEN
*
         // Test for non-positive-definiteness
*
         AJJ = REAL( A( 1, 1 ) )
         IF( AJJ.LE.ZERO.OR.SISNAN( AJJ ) ) THEN
            INFO = 1
            RETURN
         END IF
*
         // Factor
*
         A( 1, 1 ) = SQRT( AJJ )
*
      // Use recursive code
*
      ELSE
         N1 = N/2
         N2 = N-N1
*
         // Factor A11
*
         CALL CPOTRF2( UPLO, N1, A( 1, 1 ), LDA, IINFO )
         IF ( IINFO.NE.0 ) THEN
            INFO = IINFO
            RETURN
         END IF
*
         // Compute the Cholesky factorization A = U**H*U
*
         IF( UPPER ) THEN
*
            // Update and scale A12
*
            CALL CTRSM( 'L', 'U', 'C', 'N', N1, N2, CONE, A( 1, 1 ), LDA, A( 1, N1+1 ), LDA )
*
            // Update and factor A22
*
            CALL CHERK( UPLO, 'C', N2, N1, -ONE, A( 1, N1+1 ), LDA, ONE, A( N1+1, N1+1 ), LDA )
*
            CALL CPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
*
            IF ( IINFO.NE.0 ) THEN
               INFO = IINFO + N1
               RETURN
            END IF
*
         // Compute the Cholesky factorization A = L*L**H
*
         ELSE
*
            // Update and scale A21
*
            CALL CTRSM( 'R', 'L', 'C', 'N', N2, N1, CONE, A( 1, 1 ), LDA, A( N1+1, 1 ), LDA )
*
            // Update and factor A22
*
            CALL CHERK( UPLO, 'N', N2, N1, -ONE, A( N1+1, 1 ), LDA, ONE, A( N1+1, N1+1 ), LDA )
*
            CALL CPOTRF2( UPLO, N2, A( N1+1, N1+1 ), LDA, IINFO )
*
            IF ( IINFO.NE.0 ) THEN
               INFO = IINFO + N1
               RETURN
            END IF
*
         END IF
      END IF
      RETURN
*
      // End of CPOTRF2
*
      END
