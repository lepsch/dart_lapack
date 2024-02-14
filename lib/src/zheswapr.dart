      void zheswapr(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final int I1, final int I2,) {
  final A = A_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String           UPLO;
      int              I1, I2, LDA, N;
      Complex          A( LDA, N );

// =====================================================================

      bool               UPPER;
      int                I;
      Complex            TMP;

      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZSWAP

      UPPER = lsame( UPLO, 'U' );
      if (UPPER) {

          // UPPER
          // first swap
          //  - swap column I1 and I2 from I1 to I1-1
         zswap(I1-1, A(1,I1), 1, A(1,I2), 1 );

           // second swap :
           // - swap A(I1,I1) and A(I2,I2)
           // - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1
           // - swap A(I2,I1) and A(I1,I2)

         TMP=A(I1,I1);
         A(I1,I1)=A(I2,I2);
         A(I2,I2)=TMP;

         for (I = 1; I <= I2-I1-1; I++) {
            TMP=A(I1,I1+I);
            A(I1,I1+I)=DCONJG(A(I1+I,I2));
            A(I1+I,I2)=DCONJG(TMP);
         }

          A(I1,I2)=DCONJG(A(I1,I2));


           // third swap
           // - swap row I1 and I2 from I2+1 to N
         for (I = I2+1; I <= N; I++) {
            TMP=A(I1,I);
            A(I1,I)=A(I2,I);
            A(I2,I)=TMP;
         }

        } else {

          // LOWER
          // first swap
          //  - swap row I1 and I2 from 1 to I1-1
         zswap(I1-1, A(I1,1), LDA, A(I2,1), LDA );

          // second swap :
          //  - swap A(I1,I1) and A(I2,I2)
          //  - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1
          //  - swap A(I2,I1) and A(I1,I2)

          TMP=A(I1,I1);
          A(I1,I1)=A(I2,I2);
          A(I2,I2)=TMP;

          for (I = 1; I <= I2-I1-1; I++) {
             TMP=A(I1+I,I1);
             A(I1+I,I1)=DCONJG(A(I2,I1+I));
             A(I2,I1+I)=DCONJG(TMP);
          }

          A(I2,I1)=DCONJG(A(I2,I1));

          // third swap
          //  - swap col I1 and I2 from I2+1 to N
          for (I = I2+1; I <= N; I++) {
             TMP=A(I,I1);
             A(I,I1)=A(I,I2);
             A(I,I2)=TMP;
          }

      }

      END SUBROUTINE ZHESWAPR;
