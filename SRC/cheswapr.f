      SUBROUTINE CHESWAPR( UPLO, N, A, LDA, I1, I2)

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
      int                I;
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
         CALL CSWAP( I1-1, A(1,I1), 1, A(1,I2), 1 )

           // second swap :
           // - swap A(I1,I1) and A(I2,I2)
           // - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1
           // - swap A(I2,I1) and A(I1,I2)

         TMP=A(I1,I1)
         A(I1,I1)=A(I2,I2)
         A(I2,I2)=TMP

         DO I=1,I2-I1-1
            TMP=A(I1,I1+I)
            A(I1,I1+I)=CONJG(A(I1+I,I2))
            A(I1+I,I2)=CONJG(TMP)
         END DO

          A(I1,I2)=CONJG(A(I1,I2))


          t // hird swap
           // - swap row I1 and I2 from I2+1 to N
         DO I=I2+1,N
            TMP=A(I1,I)
            A(I1,I)=A(I2,I)
            A(I2,I)=TMP
         END DO

        } else {

          // LOWER
          // first swap
           // - swap row I1 and I2 from 1 to I1-1
         CALL CSWAP ( I1-1, A(I1,1), LDA, A(I2,1), LDA )

          // second swap :
           // - swap A(I1,I1) and A(I2,I2)
           // - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1
           // - swap A(I2,I1) and A(I1,I2)

          TMP=A(I1,I1)
          A(I1,I1)=A(I2,I2)
          A(I2,I2)=TMP

          DO I=1,I2-I1-1
             TMP=A(I1+I,I1)
             A(I1+I,I1)=CONJG(A(I2,I1+I))
             A(I2,I1+I)=CONJG(TMP)
          END DO

          A(I2,I1)=CONJG(A(I2,I1))

         t // hird swap
           // - swap col I1 and I2 from I2+1 to N
          DO I=I2+1,N
             TMP=A(I,I1)
             A(I,I1)=A(I,I2)
             A(I,I2)=TMP
          END DO

      ENDIF

      END SUBROUTINE CHESWAPR
