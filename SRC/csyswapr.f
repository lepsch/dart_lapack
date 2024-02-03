      SUBROUTINE CSYSWAPR( UPLO, N, A, LDA, I1, I2)

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String           UPLO;
      int              I1, I2, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX          A( LDA, N )

*  =====================================================================

      // ..
      // .. Local Scalars ..
      bool               UPPER;
      COMPLEX            TMP

      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSWAP
      // ..
      // .. Executable Statements ..

      UPPER = LSAME( UPLO, 'U' )
      if (UPPER) {

          // UPPER
          // first swap
           // - swap column I1 and I2 from I1 to I1-1
         cswap(I1-1, A(1,I1), 1, A(1,I2), 1 );

           // second swap :
           // - swap A(I1,I1) and A(I2,I2)
           // - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1
         TMP=A(I1,I1)
         A(I1,I1)=A(I2,I2)
         A(I2,I2)=TMP

         cswap(I2-I1-1, A(I1,I1+1), LDA, A(I1+1,I2), 1 );

           // third swap
           // - swap row I1 and I2 from I2+1 to N
         IF ( I2.LT.N ) CALL CSWAP( N-I2, A(I1,I2+1), LDA, A(I2,I2+1), LDA )

        } else {

          // LOWER
          // first swap
           // - swap row I1 and I2 from I1 to I1-1
         cswap(I1-1, A(I1,1), LDA, A(I2,1), LDA );

          // second swap :
           // - swap A(I1,I1) and A(I2,I2)
           // - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1
          TMP=A(I1,I1)
          A(I1,I1)=A(I2,I2)
          A(I2,I2)=TMP

          cswap(I2-I1-1, A(I1+1,I1), 1, A(I2,I1+1), LDA );

          // third swap
           // - swap col I1 and I2 from I2+1 to N
          IF ( I2.LT.N ) CALL CSWAP( N-I2, A(I2+1,I1), 1, A(I2+1,I2), 1 )

      ENDIF
      END SUBROUTINE CSYSWAPR
